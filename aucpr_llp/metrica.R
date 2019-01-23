library(randomForest)
library(SparseM)
library(grDevices)
library(foreign)
library(mlbench)
library(e1071)
library(PRROC)
library(zoo)

readlibsvm <- function(file){
  dd <- read.matrix.csr(file, fac=F)
  dtst <- as.data.frame(as.matrix(dd$x))
  dtst$Y <- dd$y+1
  dtst$Y[which(dtst$Y==0)] <- 1
  dtst$Y<-factor(dtst$Y)
  return(dtst[sample(1:nrow(dtst)),])
}

# Obtain the LPs from a full labeling by forming bags 
# of size bs with contiguous instances
obtain_lps <- function(classVar, bs=50) {
  nlabels <- length(levels(classVar))
  ninsts <- length(classVar)
  nbags <- ceiling(length(classVar)/bs)
  lps <- matrix(0,nbags,nlabels)
  
  for (i in 1:nbags) {
    first <- (i-1)*bs+1
    last <- min(i*bs,ninsts)
    lps[i,] <- as.vector(table(classVar[first:last]))
  }
  return(lps)
}


probMeanRank <- function(lps,probs) {
  res <- vector("numeric",nrow(lps))
  iAct <- 1
  for (i in 1:nrow(lps)) {
    bs <- sum(lps[i,])
    iNext <- iAct+bs
    
    vals <- c(0,as.vector(probs[iAct:(iNext-1),2]),1)
    res[i] <- mean(vals[(lps[i,1]+1):(lps[i,1]+2)])
    
    iAct <- iNext
  }
  
  return(list(m=mean(res),s=sd(res)))
}

probMinRank <- function(lps,probs) {
  res <- vector("numeric",nrow(lps))
  iAct <- 1
  for (i in 1:nrow(lps)) {
    bs <- sum(lps[i,])
    iNext <- iAct+bs
    
    vals <- c(0,as.vector(probs[iAct:(iNext-1),2]),1)
    res[i] <- vals[lps[i,1]+1]+0.001
    
    iAct <- iNext
  }
  
  return(list(m=mean(res),s=sd(res)))
}

probMaxRank <- function(lps,probs) {
  res <- vector("numeric",nrow(lps))
  iAct <- 1
  for (i in 1:nrow(lps)) {
    bs <- sum(lps[i,])
    iNext <- iAct+bs
    
    vals <- c(0,as.vector(probs[iAct:(iNext-1),2]),1)
    res[i] <- vals[lps[i,1]+2]-0.001
    
    iAct <- iNext
  }
  
  return(list(m=mean(res),s=sd(res)))
}

# PR in regular supervised setting
precisionRecall <- function(realClass,probs,draw=T) {
  realClass = realClass-1
  totalRealPos <- sum(realClass)

  thrs <- sort(unique(c(0,probs,1)),decreasing = T)

  res <- matrix(0,length(thrs),3)

  for (t in 1:length(thrs)) {
    pred <- vector("numeric",nrow(probs))
    pred[which(probs[,2]>=thrs[t])] <- 1

    TP <- sum(pred*realClass)
    res[t,]<-c(thrs[t],TP,sum(pred)-TP)
  }
  
  prPoints <- interpolaMacro(res, totalRealPos)
  
  if (draw) {
    rbPal <- colorRampPalette(c("green","red"))
    cool <- rbPal(nrow(prPoints))[as.numeric(cut(prPoints[,1],breaks = nrow(prPoints)))]
    
    par(fig=c(0,0.9,0,1))
    plot(prPoints[1,2], prPoints[1,3], pch=".", col = cool[1],
         xlab="Recall", ylab="Precision", xlim = c(0,1), ylim=c(0,1))
    for (i in 2:nrow(prPoints)) {
      lines(prPoints[(i-1):i,2], prPoints[(i-1):i,3], col=cool[i])
    }
    
    par(fig=c(0.6,1,0,1), new=TRUE)
    legend_image <- as.raster(matrix(cool, ncol=1))
    plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
    text(x=1.5, y = seq(0,1,l=5), cex=0.6, srt=270, labels=seq(0,1,l=5))
    rasterImage(legend_image, 0, 0, 1,1)
    
    mtext("PR Curve", side=3, outer=TRUE, line=-3)
  }
  return( list( AUC=sum(diff(prPoints[,2])*rollmean(prPoints[,3],2)) ) )
}


# Interpolation of the line for PR-curve using macro averaging
interpolaMacro <- function(res, totalRealPos) {
  prPoints <- matrix(0,1,3) # threshold, Recall, Precision
  prPoints[1,] <- c(res[1,1],res[1,2]/totalRealPos,res[1,2]/(res[1,2]+res[1,3]))
  for (i in 2:nrow(res)) { # interpolation
    leap <- res[i,2]-res[i-1,2]-1
    fact <- (res[i,3]-res[i-1,3])/(res[i,2]-res[i-1,2])
    factT <- (res[i,1]-res[i-1,1])/(res[i,2]-res[i-1,2])
    if (leap>0) {
      for (l in 1:leap) {
        nTP <- res[i-1,2] + l
        nFP <- res[i-1,3] + fact*l
        nT <- res[i-1,1] + factT*l
        prPoints <- rbind(prPoints,c(nT,nTP/totalRealPos,nTP/(nTP+nFP)))
      }
    } 
    prPoints <- rbind(prPoints,c(res[i,1],res[i,2]/totalRealPos,res[i,2]/(res[i,2]+res[i,3])))
  }
  prPoints[is.nan(prPoints)] <- 0
  return(prPoints)
}


# PR-curve using macro averaging
precisionRecallMacro <- function(lps,probs,draw=T) {
  totalRealPos <- sum(lps[,2])
  thrs <- sort(unique(c(0,probs,1)),decreasing = T)
  resMin <- matrix(0,length(thrs),3) # threshold, TP, FP
  resAle <- matrix(0,length(thrs),3) # threshold, TP, FP
  resMed <- matrix(0,length(thrs),3) # threshold, TP, FP
  resMax <- matrix(0,length(thrs),3) # threshold, TP, FP
  
  BSs <- rowSums(lps)
  recTria <- lps[,2]>mean(lps[,2])
  
  for (t in 1:length(thrs)) {
    pred <- vector("numeric",nrow(probs))
    pred[which(probs[,2]>=thrs[t])] <- 1
    totalPredPos <- sum(pred)

    PredPos <- vector("numeric",length=nrow(lps))
    iAct <- 1
    for (b in 1:nrow(lps)){
      iNext <- iAct+BSs[b]
      PredPos[b] <- sum(pred[iAct:(iNext-1)])
      iAct <- iNext
    }
    precTria <- PredPos>mean(PredPos)
    
    tpMin <- vector("numeric",length=nrow(lps))
    tpAle <- vector("numeric",length=nrow(lps))
    tpMed <- vector("numeric",length=nrow(lps))
    tpMax <- vector("numeric",length=nrow(lps))

    for (b in 1:nrow(lps)){
      tpMin[b] <- max( 0 , PredPos[b]-lps[b,1] )
      tpAle[b] <- PredPos[b]*lps[b,2]/BSs[b]
      tpMax[b] <- min( lps[b,2] , PredPos[b] ) 
    }
    pesos <- (tpMax-tpMin)/(2*max(tpMax-tpMin))
    pesos[is.nan(pesos)] <- 0
    pesos <- 1-pesos
    recprev <- sum(pesos[recTria]*tpMin[recTria])/sum(pesos[recTria]*lps[recTria,2])
    precprev <- sum(pesos[precTria]*tpMin[precTria])/sum(pesos[precTria]*PredPos[precTria])
    if (is.nan(precprev) | precprev < 0 | precprev > 1) { precprev <- 0}
    for (b in 1:nrow(lps)) {
      prevision<-min(max(tpAle[b],precprev*PredPos[b],recprev*lps[b,2]),tpMax[b])
      tpMed[b] <- round(mean(c(tpMax[b],prevision)))#round((tpMax[b]-tpMin[b])*wei[b]+tpMin[b])#
    }
    resMin[t,]<-c(thrs[t],sum(tpMin),totalPredPos-sum(tpMin))
    resAle[t,]<-c(thrs[t],sum(tpAle),totalPredPos-sum(tpAle))
    resMed[t,]<-c(thrs[t],sum(tpMed),totalPredPos-sum(tpMed))
    resMax[t,]<-c(thrs[t],sum(tpMax),totalPredPos-sum(tpMax))
  }
  prMinPoints <- interpolaMacro(resMin, totalRealPos)
  prAlePoints <- interpolaMacro(resAle, totalRealPos)
  prMaxPoints <- interpolaMacro(resMax, totalRealPos)
  prMedPoints <- interpolaMacro(resMed, totalRealPos)

  if (draw) {
    rbPal <- colorRampPalette(c("green","red"))
    cool <- rbPal(nrow(prMinPoints))[as.numeric(cut(prMinPoints[,1],breaks = nrow(prMinPoints)))]
    
    par(fig=c(0,0.9,0,1))
    plot(prMinPoints[1,2], prMinPoints[1,3], pch=".", col = cool[1],
         xlab="Recall", ylab="Precision", xlim = c(0,1), ylim=c(0,1))
    for (i in 2:nrow(prMinPoints)) {
      lines(prMinPoints[(i-1):i,2], prMinPoints[(i-1):i,3], col=cool[i])
    }
    polygon(c(prMaxPoints[,2],rev(prMinPoints[,2])), c(prMaxPoints[,3],rev(prMinPoints[,3])),col='gray90',border=NA)
    
    cool <- rbPal(nrow(prAlePoints))[as.numeric(cut(prAlePoints[,1],breaks = nrow(prAlePoints)))]
    
    par(fig=c(0,0.9,0,1), new=TRUE)
    plot(prAlePoints[1,2], prAlePoints[1,3], pch=".", col=cool[1],
         xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), axes = F)
    for (i in 2:nrow(prAlePoints)) {
      lines(prAlePoints[(i-1):i,2],prAlePoints[(i-1):i,3],col = cool[i])
    }
    
    cool <- rbPal(nrow(prMedPoints))[as.numeric(cut(prMedPoints[,1],breaks = nrow(prMedPoints)))]
    
    par(fig=c(0,0.9,0,1), new=TRUE)
    plot(prMedPoints[1,2], prMedPoints[1,3], pch=".", col=cool[1],
         xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), axes = F)
    for (i in 2:nrow(prMedPoints)) {
      lines(prMedPoints[(i-1):i,2],prMedPoints[(i-1):i,3],col = cool[i])
    }
    
    cool <- rbPal(nrow(prMaxPoints))[as.numeric(cut(prMaxPoints[,1],breaks = nrow(prMaxPoints)))]
    
    par(fig=c(0,0.9,0,1), new=TRUE)
    plot(prMaxPoints[1,2], prMaxPoints[1,3], pch=".", col=cool[1],
         xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), axes = F)
    for (i in 2:nrow(prMaxPoints)) {
      lines(prMaxPoints[(i-1):i,2],prMaxPoints[(i-1):i,3],col = cool[i])
    }
  
    par(fig=c(0.6,1,0,1), new=TRUE)
    legend_image <- as.raster(matrix(cool, ncol=1))
    plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
    text(x=1.5, y = seq(0,1,l=5), cex=0.6, srt=270, labels=seq(0,1,l=5))
    rasterImage(legend_image, 0, 0, 1,1)
    
    mtext("PR Curve", side=3, outer=TRUE, line=-3)
  }

  return( list(
    AUCmin = sum(diff(prMinPoints[,2])*rollmean(prMinPoints[,3],2)), 
    AUCale = sum(diff(prAlePoints[,2])*rollmean(prAlePoints[,3],2)), 
    AUCmed = sum(diff(prMedPoints[,2])*rollmean(prMedPoints[,3],2)), 
    AUCmax = sum(diff(prMaxPoints[,2])*rollmean(prMaxPoints[,3],2)) 
  ) )
}


# PR-curve using micro averaging
precisionRecallMicro <- function(lps,probs,draw=T) {
  totalRealPos <- sum(lps[,2])
  prodRealPos <- prod(lps[which(lps[,2]>0),2])
  numPRP <- sum(lps[,2]>0)
  thrs <- sort(unique(c(0,probs,1)),decreasing = T)
  resMin <- matrix(0,length(thrs),3) # threshold, TP, FP
  resAle <- matrix(0,length(thrs),3) # threshold, TP, FP
  resMed <- matrix(0,length(thrs),3) # threshold, TP, FP
  resMax <- matrix(0,length(thrs),3) # threshold, TP, FP
  
  BSs <- rowSums(lps)

  recTria <- lps[,2]>mean(lps[,2])
  
  for (t in 1:length(thrs)) {
    pred <- vector("numeric",nrow(probs))
    pred[which(probs[,2]>=thrs[t])] <- 1
    totalPredPos <- sum(pred)
    
    PredPos <- vector("numeric",length=nrow(lps))
    iAct <- 1
    for (b in 1:nrow(lps)){
      iNext <- iAct+BSs[b]
      PredPos[b] <- sum(pred[iAct:(iNext-1)])
      iAct <- iNext
    }
    
    tpMin <- vector("numeric",length=nrow(lps))
    tpAle <- vector("numeric",length=nrow(lps))
    tpMed <- vector("numeric",length=nrow(lps))
    tpMax <- vector("numeric",length=nrow(lps))
    iAct <- 1
    for (b in 1:nrow(lps)){
      iNext <- iAct+BSs[b]
      tpMin[b] <- max( 0 , PredPos[b]-lps[b,1] )
      tpAle[b] <- PredPos[b]*lps[b,2]/BSs[b]
      tpMax[b] <- min( lps[b,2], PredPos[b] )
      iAct <- iNext
    }

    precTria <- PredPos>mean(PredPos)
    pesos <- (tpMax-tpMin)/(2*max(tpMax-tpMin))
    pesos[is.nan(pesos)] <- 0
    pesos <- 1-pesos
    recprev <- sum(pesos[recTria]*tpMin[recTria])/sum(pesos[recTria]*lps[recTria,2])
    precprev <- sum(pesos[precTria]*tpMin[precTria])/sum(pesos[precTria]*PredPos[precTria])
    if (is.nan(precprev) | precprev < 0 | precprev > 1) { precprev <- 0}
    for (b in 1:nrow(lps)) {
      prevision<-min(max(tpAle[b],precprev*PredPos[b],recprev*lps[b,2]),tpMax[b])
      tpMed[b] <- round(mean(c(tpMax[b],prevision)))
    }
    
    
    meanPrec <- tpMin/PredPos
    meanPrec[is.na(meanPrec)] <- 0
    meanRec <- tpMin/lps[,2]
    meanRec[is.na(meanRec)] <- 0
    resMin[t,]<-c( thrs[t], 100*mean(meanPrec[which(PredPos>0)]) , 100*mean(meanRec[which(lps[,2]>0)]) )
    meanPrec <- tpAle/PredPos
    meanPrec[is.na(meanPrec)] <- 0
    meanRec <- tpAle/lps[,2]
    meanRec[is.na(meanRec)] <- 0
    resAle[t,]<-c( thrs[t], 100*mean(meanPrec[which(PredPos>0)]) , 100*mean(meanRec[which(lps[,2]>0)]) )
    meanPrec <- tpMed/PredPos
    meanPrec[is.na(meanPrec)] <- 0
    meanRec <- tpMed/lps[,2]
    meanRec[is.na(meanRec)] <- 0
    resMed[t,]<-c( thrs[t], 100*mean(meanPrec[which(PredPos>0)]) , 100*mean(meanRec[which(lps[,2]>0)]) )
    meanPrec <- tpMax/PredPos
    meanPrec[is.na(meanPrec)] <- 0
    meanRec <- tpMax/lps[,2]
    meanRec[is.na(meanRec)] <- 0
    resMax[t,]<-c( thrs[t], 100*mean(meanPrec[which(PredPos>0)]) , 100*mean(meanRec[which(lps[,2]>0)]) )
  }
  resMin[is.nan(resMin)]<-0
  resAle[is.nan(resAle)]<-0
  resMed[is.nan(resMed)]<-0
  resMax[is.nan(resMax)]<-0
  # print(resMin)
  
  prMinPoints <- interpolaMicro(resMin)
  prAlePoints <- interpolaMicro(resAle)
  prMedPoints <- interpolaMicro(resMed)
  prMaxPoints <- interpolaMicro(resMax)

  if (draw) {
    rbPal <- colorRampPalette(c("green","red"))
    cool <- rbPal(nrow(prMinPoints))[as.numeric(cut(prMinPoints[,1],breaks = nrow(prMinPoints)))]
    
    par(fig=c(0,0.9,0,1))
    plot(prMinPoints[1,2], prMinPoints[1,3], pch=".", col = cool[1],
         xlab="Recall", ylab="Precision", xlim = c(0,1), ylim=c(0,1))
    for (i in 2:nrow(prMinPoints)) {
      lines(prMinPoints[(i-1):i,2], prMinPoints[(i-1):i,3], col=cool[i])
    }
    polygon(c(prMaxPoints[,2],rev(prMinPoints[,2])), c(prMaxPoints[,3],rev(prMinPoints[,3])),col='gray90',border=NA)
    
    cool <- rbPal(nrow(prAlePoints))[as.numeric(cut(prAlePoints[,1],breaks = nrow(prAlePoints)))]
    
    par(fig=c(0,0.9,0,1), new=TRUE)
    plot(prAlePoints[1,2], prAlePoints[1,3], pch=".", col=cool[1],
         xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), axes = F)
    for (i in 2:nrow(prAlePoints)) {
      lines(prAlePoints[(i-1):i,2],prAlePoints[(i-1):i,3],col = cool[i])
    }
    
    cool <- rbPal(nrow(prMedPoints))[as.numeric(cut(prMedPoints[,1],breaks = nrow(prMedPoints)))]
    
    par(fig=c(0,0.9,0,1), new=TRUE)
    plot(prMedPoints[1,2], prMedPoints[1,3], pch=".", col=cool[1],
         xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), axes = F)
    for (i in 2:nrow(prMedPoints)) {
      lines(prMedPoints[(i-1):i,2],prMedPoints[(i-1):i,3],col = cool[i])
    }
    
    cool <- rbPal(nrow(prMaxPoints))[as.numeric(cut(prMaxPoints[,1],breaks = nrow(prMaxPoints)))]
    
    par(fig=c(0,0.9,0,1), new=TRUE)
    plot(prMaxPoints[1,2], prMaxPoints[1,3], pch=".", col=cool[1],
         xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), axes = F)
    for (i in 2:nrow(prMaxPoints)) {
      lines(prMaxPoints[(i-1):i,2],prMaxPoints[(i-1):i,3],col = cool[i])
    }
    
    
    
    # BAR LEGEND Threshold
    par(fig=c(0.6,1,0,1), new=TRUE)
    legend_image <- as.raster(matrix(cool, ncol=1))
    plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
    text(x=1.5, y = seq(0,1,l=5), cex=0.6, srt=270, labels=seq(0,1,l=5))
    rasterImage(legend_image, 0, 0, 1,1)
    
    mtext("PR Curve", side=3, outer=TRUE, line=-3)
  }
  return( list(
    AUCmin = sum(diff(prMinPoints[,2])*rollmean(prMinPoints[,3],2)), 
    AUCale = sum(diff(prAlePoints[,2])*rollmean(prAlePoints[,3],2)), 
    AUCmed = sum(diff(prMedPoints[,2])*rollmean(prMedPoints[,3],2)), 
    AUCmax = sum(diff(prMaxPoints[,2])*rollmean(prMaxPoints[,3],2)) 
    ) )
}

# Interpolation of the line for PR-curve using micro averaging
interpolaMicro <- function(res){
  prPoints <- matrix(0,1,3) # threshold, Recall, Precision
  prPoints[1,] <- c(res[1,1],res[1,3]/100,res[1,2]/100)
  for (i in 2:nrow(res)) { # interpolation
    leap <- res[i,2]-res[i-1,2]-1
    factF <- (100-res[i,2]-100+res[i-1,2])/(res[i,2]-res[i-1,2])
    factRec <- (res[i,3]-res[i-1,3])/(res[i,2]-res[i-1,2])
    factT <- (res[i,1]-res[i-1,1])/(res[i,2]-res[i-1,2])
    if (leap>=1) {
      for (l in 1:leap) {
        nTP <- res[i-1,2] + l
        nFP <- (100-res[i-1,2]) + l*factF
        nTPrec <- res[i-1,3] + factRec*l
        nT <- res[i-1,1] + factT*l
        prPoints <- rbind(prPoints,c(nT,nTPrec/100,nTP/(nTP+nFP)))
      }
    }
    prPoints <- rbind(prPoints,c(res[i,1],res[i,3]/100,res[i,2]/100))
  }
  prPoints[is.nan(prPoints)]<-0
  return (prPoints)
}



############################ DATA LOAD ############################




# dtst <- read.table("../_data_binary_class/biodeg.csv", sep=";", header=F)
# iclass <- 42
# # 1055 x 42 (66/34)
# 
# Pima Indians Diabetes Database
data(PimaIndiansDiabetes)
dtst <- PimaIndiansDiabetes
dtst <- dtst[sample(1:nrow(dtst)),]
iclass <- 9
dtst[,iclass]<-factor(as.numeric(dtst[,iclass]))
rm(PimaIndiansDiabetes)
# 768 x 9 (65/35)
# 
# # Titanic
# data(Titanic)
# Titanic_df=as.data.frame(Titanic)
# repeating_sequence <- rep.int(seq_len(nrow(Titanic_df)), Titanic_df$Freq)
# dtst <- Titanic_df[repeating_sequence,]
# dtst$Freq=NULL
# dtst<-dtst[sample(1:nrow(dtst)),]
# iclass <- 4
# dtst[,iclass]<-factor(as.numeric(dtst[,iclass]))
# rm(Titanic,Titanic_df,repeating_sequence)
# # 2201 x 4 (68/32)
# 
# dtst <- readlibsvm("../_data_binary_class/australian.libsvm")
# iclass <- 15
# # 690 x 15 (56/44)
# 
# dtst <- readlibsvm("../_data_binary_class/splice.libsvm")
# iclass <- 61
# # 3175 x 61 (48/52)
# 
# dtst <- readlibsvm("../_data_binary_class/svmguide1.libsvm")
# iclass <- 5
# # 7089 x 5 (44/56)
# 
# # Wisconsin Breast Cancer Database
# data(BreastCancer)
# dtst <- BreastCancer
# dtst <- dtst[sample(1:nrow(dtst)),-1]
# for (i in 1:9) {
#   dtst[,i] <- as.numeric(dtst[,i])
# }
# iclass <- 10
# dtst[,iclass]<-factor(as.numeric(dtst[,iclass]))
# rm(BreastCancer)
# # 699 x 10 (66/34)




###################### LEARNING ######################

frm <- as.formula(paste(names(dtst)[iclass]," ~ .",sep=""))

fit=naiveBayes(frm, data=dtst)
origClassProbs <- predict(fit, dtst, type="raw")+0.002
origClassProbs <- origClassProbs/rowSums(origClassProbs)


namefile <- "phishing"
lPrSw <- c(0.1,0.2,0.3,0.4)
lBS <- c(5,10,25,50)
dr <- F
nrep <- 10

resultados <- c()
for (bs in lBS) {
  for (propSwaps in lPrSw) {
    
    swaps <- sort(sample(1:nrow(origClassProbs),size = propSwaps*nrow(origClassProbs),replace = F))
    classProbs <- origClassProbs
    classProbs[swaps,] <- classProbs[swaps,c(2,1)]
    
    resMacro <- matrix(0,nrow=nrep,ncol=4)
    resMicro <- matrix(0,nrow=nrep,ncol=4)
    resLoss <- matrix(0,nrow=nrep,ncol=1)
    
    for (r in 1:nrep) {
      print(c(bs,propSwaps,r))
      orden <- sample(1:nrow(dtst))
      real_lps <- obtain_lps(dtst[orden,iclass],bs = bs)
      pred_lps <- obtain_lps(factor((classProbs[orden,2]>=0.5)+1),bs = bs)
      resMacro[r,] <- unlist(precisionRecallMacro(real_lps,classProbs[orden,], draw=dr))
      resMicro[r,] <- unlist(precisionRecallMicro(real_lps,classProbs[orden,], draw=dr))
      resLoss[r,1] <- 1-sum(abs(real_lps[,2]-pred_lps[,2]))/nrow(dtst)
      
    }
    
    resultados <- rbind(resultados,
                        c(colSums(resMacro)/nrow(resMacro),
                          colSums(resMicro)/nrow(resMicro),
                          mean(resLoss),
                          sum(diag(table(dtst[,iclass], (classProbs[,2]>=0.5)+1)))/nrow(dtst),
                          precisionRecall(as.numeric(dtst[,iclass]),classProbs, draw=dr))
    )
  }
}

write.csv(resultados,paste("results/",namefile,".resout",sep=""),row.names = F,quote=F)

