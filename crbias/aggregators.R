
library(gtools)


Amean <- function(tab) {
  vectrec <- diag(tab)/apply(tab,1,sum)
  vectrec[which(is.nan(vectrec)|is.na(vectrec))]<-0
  
  return ( mean(vectrec) )  
}

Acc <- function(tab) {
  return( sum(diag(tab))/(sum(tab)) )
}

f1 <- function(tab) {
  vectprec <- diag(tab)/apply(tab,2,sum)
  vectprec[which(is.nan(vectprec)|is.na(vectprec))]<-0
  prec <- mean(vectprec)
  
  vectrec <- diag(tab)/apply(tab,1,sum)
  vectrec[which(is.nan(vectrec)|is.na(vectrec))]<-0
  rec <- mean(vectrec) 
  
  if ((prec + rec) != 0) {
    return (2*prec*rec / (prec+rec) )
  } else {
    return ( 0 )
  }
}


#################################### K-Means ###################################

processDataForKmeans <- function(matriz,nLabels) {
  mat <- matrix(0,nrow(matriz),nLabels+1)
  for (i in 1:nrow(matriz)) {
    mat[i,]<-processDataForKmeansInner(matriz[i,],nLabels)
  }
  return (mat);
}
processDataForKmeansInner <- function(vector,nLabels) {
  tab <- table(factor(vector,1:nLabels))
  vct <- c( as.vector(tab/sum(tab)) , 0 )
  
  for (i in 2:nLabels) {
    vct[nLabels+1] <- vct[nLabels+1] + (vct[i]-vct[i-1])
  }
  vct[nLabels+1] <- vct[nLabels+1] / nLabels
  return (vct)
}

realMatchClusterClass <- function(tab,firstAssignment,nLabels) {
  lClusters <- vector("numeric",nLabels)
  reassignmnt <- matrix(0,nLabels,nLabels)
  for (l in 1:nLabels){
    indsOfClass <- which(firstAssignment==l)
    lClusters[l] <- length(indsOfClass)
    reassignmnt[l,]<-apply(tab[indsOfClass,1:nLabels,drop=F],2,sum)
  }
  indexes <- sort(lClusters, index.return = TRUE)$ix
  
  realKs <- vector("numeric",nLabels)
  for (l in nLabels:1) {
    realKs[ indexes[l] ] <- which.max(reassignmnt[ indexes[l], ] )
    reassignmnt[, realKs[ indexes[l] ] ] <- -1 # avoiding other l's taking the same label
  }
  
  resKmeans <- factor(rep(1,nrow(tab)), levels=1:nLabels)
  for (l in 1:nLabels){
    resKmeans[ which(firstAssignment==l) ] <- realKs[l]
  }
  
  return(resKmeans)
}

areEquals <- function(tab,centers) {
  prim <- which( centers == centers[anyDuplicated(centers)] )
  if (length(prim > 0)) return ( matrix(prim,1,length(prim)) )
  allCombs <- combinations(n = length(centers), r = 2, v = centers, repeats.allowed = F)
  return( allCombs[apply(allCombs,1,function(X,Y) return(all(Y[X[1],]==Y[X[2],])),tab), ,drop=F] )
}

pickUpCenters <- function(tab,nLabels) {
  listsForCenters<-apply(tab,2,function(X) return(sort(X,index.return=T,decreasing=T)$ix))[,1:nLabels]
  iCenters <- rep(1,nLabels)
  centers <- sapply(1:nLabels,function(X,Y,Z) return(Y[Z[X],X]),listsForCenters,iCenters)
  
  mEquals <- areEquals(tab,centers)
  if (nrow(mEquals)!=0) {
    repeat{
      inds <- mEquals[1,]
      inds <- inds[sample(1:length(inds), 1)]
      iCenters[inds]<-iCenters[inds]+1
      centers[inds] <- listsForCenters[iCenters[inds],inds]
      
      mEquals <- areEquals(tab,centers)
      if (length(unique(centers)) == length(centers) & (nrow(mEquals)==0)) {
        break
      }
    }
  }

  return(centers)
}

runKmeans <- function(tab,nLabels,grTruth=NULL,alg="Lloyd") {
  centers <- pickUpCenters(tab,nLabels)
  
  cl <- kmeans(tab, tab[centers,],algorithm = alg)
  preKmeans <- factor(cl$cluster, levels=1:nLabels)
  resKmeans <- realMatchClusterClass(tab,preKmeans,nLabels)
  
  if (is.null(grTruth)) {
    return( table(resKmeans) )
  } else {
    return( table(grTruth,resKmeans) )
  }
}


###################################### MD ######################################

individualAbsoluteDifferences <- function(myTab,nLabels,alg="add") { 
  # calculate mean proportions
  props <- apply(myTab,2,sum)/nrow(myTab)
  
  resMine <- factor( rep(1,nrow(myTab)) , levels=1:nLabels )
  
  if ( alg == "add") {
    for (i in 1:nrow(myTab)) {
      resMine[i] <- which.max(myTab[i,]-props)
    }
  } else if ( alg == "mult") {
    for (i in 1:nrow(myTab)) {
      resMine[i] <- which.max(myTab[i,]/props)
    }
  }
  
  return( resMine )
}
runIndAbsDiff <- function(tab,nLabels,grTruth=NULL,alg="add") {
  
  resMine <- individualAbsoluteDifferences(tab[,1:nLabels],nLabels,alg)
  
  if (is.null(grTruth)) {
    return( table(resMine) )
  } else {
    return( table(grTruth,resMine) )
  }
}

##################################################################################


###################################### MV ######################################


MV <- function(matriz,nLabels) {
  vct <- vector("numeric",nrow(matriz))
  for (i in 1:nrow(matriz)) {
    vct[i]<-MVinner(matriz[i,],nLabels)
  }
  return (vct);
}

MVinner <- function(vector,nLabels) {
  vct <- as.vector(table(factor(vector,1:nLabels)))
  
  vctMax<-which(vct==max(vct))
  return (vctMax[sample(1:length(vctMax),1)])
}

runMV <- function(labels,nLabels,grTruth=NULL,rep=10) {
  resMV <- factor(levels=1:nLabels)
  
  for (i in 1:rep) resMV <- factor( c(resMV,MV(labels,nLabels)), levels=1:nLabels)
  
  if (is.null(grTruth)) {
    return( table(resMV)/rep )
  } else {
    nInsts <- length(grTruth)
    tabMV <- matrix(0,nLabels,nLabels)
    for (i in 1:rep) {
      tabMV <- tabMV+table(grTruth,resMV[(nInsts*(i-1)+1):(nInsts*i)])
    }
    
    return( tabMV )
  }
}

################################################################################


########################### Weighted MV - confMatrix ###########################


weightsFromIndAbsDiff <- function(matriz, nLabels, global = FALSE, alg="add", weiType="confMatrix") {
  weis <- list()
  if (global) {
    myTab <- processDataForKmeans(matriz,nLabels)[,1:nLabels]
    vct <- individualAbsoluteDifferences(myTab,nLabels,alg)
    
    for (a in 1:ncol(matriz)) {
      weis[[a]]<- table(factor(matriz[,a],1:nLabels),factor(vct,1:nLabels))
      weis[[a]]<- weis[[a]]/apply(weis[[a]],1,sum)
      if (weiType == "perClass") {
        weis[[a]]<- diag(weis[[a]])
      }
    }
  } else {
    for (a in 1:ncol(matriz)) {
      myTab <- processDataForKmeans(matriz[,-a],nLabels)[,1:nLabels]
      vct <- individualAbsoluteDifferences(myTab,nLabels,alg)
      
      weis[[a]]<- table(factor(matriz[,a],1:nLabels),factor(vct,1:nLabels))
      weis[[a]]<- weis[[a]]/apply(weis[[a]],1,sum)
      if (weiType == "perClass") {
        weis[[a]] <- diag(weis[[a]])
      }
    }
  }
  return (weis);
}

weightsFromMV <- function(matriz, nLabels, global = FALSE, weiType="confMatrix") {
  weis <- list()
  if (global){ 
    vct <- MV(matriz,nLabels) #vector("numeric", nrow(matriz))
    
    for (a in 1:ncol(matriz)) {
      weis[[a]]<- table(factor(matriz[,a],1:nLabels),factor(vct,1:nLabels))
      weis[[a]]<- weis[[a]]/apply(weis[[a]],1,sum)
      if (weiType == "perClass") {
        weis[[a]]<- diag(weis[[a]])
      }
    }
  } else {
    for (a in 1:ncol(matriz)) {
      vct <- MV(matriz[,-a],nLabels)
      
      weis[[a]]<- table(factor(matriz[,a],1:nLabels),factor(vct,1:nLabels))
      weis[[a]]<- weis[[a]]/apply(weis[[a]],1,sum)
      if (weiType == "perClass") {
        weis[[a]]<- diag(weis[[a]])
      }
    }
  }
  return (weis);
}

weightsFromKmeans <- function(matriz, nLabels, global = FALSE, weiType="confMatrix") {
  weis <- list()
  if (global) {
    tab <- processDataForKmeans(matriz,nLabels)
    centers <- pickUpCenters(tab,nLabels)
    
    preKmeans<-kmeans(tab, tab[centers,])$cluster
    vct <- realMatchClusterClass(tab,preKmeans,nLabels)
    
    for (a in 1:ncol(matriz)) {
      weis[[a]]<- table(factor(matriz[,a],1:nLabels),factor(vct,1:nLabels))
      weis[[a]]<- weis[[a]]/apply(weis[[a]],1,sum)
      if (weiType == "perClass") {
        weis[[a]]<- diag(weis[[a]])
      }
    }
  } else {
    for (a in 1:ncol(matriz)) {
      tab <- processDataForKmeans(matriz[,-a],nLabels)
      centers <- pickUpCenters(tab,nLabels)
      
      preKmeans<-kmeans(tab, tab[centers,])$cluster
      vct <- realMatchClusterClass(tab,preKmeans,nLabels)
      
      weis[[a]]<- table(factor(matriz[,a],1:nLabels),factor(vct,1:nLabels))
      weis[[a]]<- weis[[a]]/apply(weis[[a]],1,sum)
      if (weiType == "perClass") {
        weis[[a]] <- diag(weis[[a]])
      }
    }
  }
  return (weis);
}

weightedMV <- function(matriz, weis,nLabels,weiType="confMatrix") {
  vct <- vector("numeric", nrow(matriz))
  if (weiType == "confMatrix") {
    for (i in 1:nrow(matriz)) {
      vct[i]<-weightedMVinnerMatrix(matriz[i,],weis,nLabels)
    }
  } else if (weiType == "perClass") {
    for (i in 1:nrow(matriz)) {
      vct[i]<-weightedMVinnerPerClass(matriz[i,],weis,nLabels)
    }
  }
  return (factor(vct,1:nLabels));
}
weightedMVinnerMatrix <- function(vctr,weis, nLabels) {
  pesos <- vector(mode = "numeric", length = nLabels);
  inds <- which(!is.na(vctr))
  for (a in 1:length(inds)) {
    pesos <- pesos + as.vector( weis[[ inds[a] ]][ vctr[ inds[a] ] , ] )
  }
  
  pesoMax<-which(pesos==max(pesos))
  return (pesoMax[sample(1:length(pesoMax),1)])
}
weightedMVinnerPerClass <- function(vctr,weis, nLabels) {
  pesos <- vector(mode = "numeric", length = nLabels);
  inds <- which(!is.na(vctr))
  for (a in 1:length(inds)) {
    lbl <- vctr[ inds[a] ]
    pesos[lbl] <- pesos[lbl] + weis[[ inds[a] ]][ lbl ]
  }
  
  pesoMax<-which(pesos==max(pesos))
  return (pesoMax[sample(1:length(pesoMax),1)])
}

runWeightedMV <- function(labels,nLabels,grTruth=NULL,alg="MV",weiType="confMatrix",rep=10, global = FALSE) {
  weis <- list()
  
  if (alg == "kmeans") {
    weis <- weightsFromKmeans(labels, nLabels, global, weiType)
  } else if (alg == "MV") {
    weis <- weightsFromMV(labels, nLabels, global, weiType)
  } else if (alg == "absDiffAdd") {
    weis <- weightsFromIndAbsDiff(labels, nLabels, global, alg="add", weiType)
  } else if (alg == "absDiffMult") {
    weis <- weightsFromIndAbsDiff(labels, nLabels, global, alg="mult", weiType)
  }
  
  resMV <- factor(levels=1:nLabels)
  for (i in 1:rep) resMV <- factor( c(resMV, weightedMV(labels , weis , nLabels, weiType)), levels=1:nLabels)
  
  if (is.null(grTruth)) {
    return( table(resMV)/rep )
  } else {
    nInsts <- length(grTruth)
    tabMV <- matrix(0,nLabels,nLabels)
    for (i in 1:rep) {
      tabMV <- tabMV+table(grTruth,resMV[(nInsts*(i-1)+1):(nInsts*i)])
    }

    return( tabMV )
  }
}
################################################################################

