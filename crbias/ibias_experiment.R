source("aggregators.R")

trickQ<-function(q, nLabels, b=4){
  newQ <- q
  
  for (c in 1:nLabels) {
    ib <- b-c+1
    if (ib == 2) {
      next;
    } else if (ib == 1) {
      aux <- mean(q[c,3:nLabels])
      if ((q[c,2]-aux) > 0) {
        q[c,1] <- q[c,1] + (q[c,2]-aux)
        q[c,2] <- q[c,2]-aux
      }
    } else if (ib < 1) {
      ib <- nLabels+ib
      newQ[c,ib] <- q[c,2]
      newQ[c,2] <- q[c,ib]
    } else {
      newQ[c,ib] <- q[c,2]
      newQ[c,2] <- q[c,ib]
    }
  }

  return(newQ)
}

run <- function(nLabels,nAnnots,relevance,degree) {
    
  alphs <- c(relevance, rep(1,nLabels-1))
  alphs <- 10*alphs/sum(alphs)
  altAlphs <- alphs
  altAlphs[2] <- alphs[1] - (alphs[1] - alphs[2])*0.25
  
  nAltAnnots <- round(nAnnots*degree)
  
  nSamples <- 1000
  nRep <- 10
  
  vals <- matrix(0,nRep*nRep,30*3)
  
  for (k1 in 1:nRep) {

    repeat{
      p<-as.vector(rdirichlet(1,alpha=rep(1,nLabels)))
      if (min(p)>0.01) break;
    }
    
    grTruth <- factor(apply(rmultinom(nSamples,1,p),2,which.max),1:nLabels)
    
    for (k2 in 1:nRep) {

      labels <- matrix(0,nSamples,nAnnots)
  
      if ((nAltAnnots > 0) && (nAltAnnots <= nAnnots)) {
        for (a in 1:nAltAnnots) {
          q <- rdirichlet(nLabels,alpha=altAlphs)
          q <- trucarQ(q,nLabels)
          for (c in 1:nLabels) {
            indsClass <- which(grTruth == c)
            modification <- apply(rmultinom(length(indsClass),1,q[c,]),2,which.max)-1
            labels[indsClass,a] <- ((c+modification-1)%%nLabels)+1
          }
        }
      }
      if ((nAltAnnots >= 0) && (nAltAnnots < nAnnots)) {
        for (a in (nAltAnnots+1):nAnnots) {
          q <- rdirichlet(nLabels,alpha=alphs)
          for (c in 1:nLabels) {
            indsClass <- which(grTruth == c)
            modification <- apply(rmultinom(length(indsClass),1,q[c,]),2,which.max)-1
            labels[indsClass,a] <- ((c+modification-1)%%nLabels)+1
          }
        }
      }

      exp <- 0
      
      tab <- processDataForKmeans(labels,nLabels)

      res<-matrix(NA,2,2)
      res <- runKmeans(tab,nLabels,grTruth,"Hartigan-Wong")
      exp <- exp+1
      vals[(k1-1)*nRep+k2,(-2:0)+exp*3] <- c(Acc(res),f1(res),Amean(res))
      
      res[]<-NA
      res <- runMV(labels,nLabels,grTruth)
      exp <- exp+1
      vals[(k1-1)*nRep+k2,(-2:0)+exp*3] <- c(Acc(res),f1(res),Amean(res))
  
      res[]<-NA
      res <- runIndAbsDiff(tab,nLabels,grTruth, alg="add")
      exp <- exp+1
      vals[(k1-1)*nRep+k2,(-2:0)+exp*3] <- c(Acc(res),f1(res),Amean(res))
      
      res[]<-NA
      res <- runIndAbsDiff(tab,nLabels,grTruth, alg="mult")
      exp <- exp+1
      vals[(k1-1)*nRep+k2,(-2:0)+exp*3] <- c(Acc(res),f1(res),Amean(res))
  
      res[]<-NA
      res <- runWeightedMV(labels,nLabels,grTruth,alg="kmeans",weiType="perClass", global = TRUE)
      exp <- exp+1
      vals[(k1-1)*nRep+k2,(-2:0)+exp*3] <- c(Acc(res),f1(res),Amean(res))
  
      res[]<-NA
      res <- runWeightedMV(labels,nLabels,grTruth,alg="MV",weiType="perClass", global = TRUE)
      exp <- exp+1
      vals[(k1-1)*nRep+k2,(-2:0)+exp*3] <- c(Acc(res),f1(res),Amean(res))
  
      res[]<-NA
      res <- runWeightedMV(labels,nLabels,grTruth,alg="absDiffAdd",weiType="perClass", global = TRUE)
      exp <- exp+1
      vals[(k1-1)*nRep+k2,(-2:0)+exp*3] <- c(Acc(res),f1(res),Amean(res))
  
      res[]<-NA
      res <- runWeightedMV(labels,nLabels,grTruth,alg="absDiffMult",weiType="perClass", global = TRUE)
      exp <- exp+1
      vals[(k1-1)*nRep+k2,(-2:0)+exp*3] <- c(Acc(res),f1(res),Amean(res))

    }
  }
  means <- apply(vals,2,mean,na.rm = T)
  stdev <- apply(vals,2,sd,na.rm = T)
  
  return (cbind(means,stdev))
}

################################################################## 

nLabels <- 5
nAnnots <- 6
relevance <- 3
degree <- 0


print( run(nLabels,nAnnots,relevance,degree) )
  