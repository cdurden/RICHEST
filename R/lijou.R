# ascending factorial function
fLogAscFact <- function(a,n) {
  if (n > 0) {
    return(sum(log(a+c(0:(n-1)))))
  } else {
    return(0)
  }
}

# logarithm of pitman sampling formula
fLogPitman <- function(x,vSpeciesCount,sampleRichness,sampleSize) {
  if (x[2] > 1 | x[1] < -1 | x[1] < -x[2] | x[1] < (1-sampleRichness)*x[2]) {
    return(-Inf)
  } else {
    return(sum(log(x[1] + c(1:(sampleRichness-1))*x[2]))-fLogAscFact(x[1]+1,sampleSize-1)+sum(mapply(fLogAscFact,rep(1-x[2],sampleRichness),vSpeciesCount-1)))
  }
}

fGradLogPitman <- function(x,vSpeciesCount,sampleRichness,sampleSize) {
  if (x[2] > 1 | x[1] < -1 | x[1] < -x[2] | x[1] < (1-sampleRichness)*x[2]) {
    return(NA)
  } else {
    d <- c(0,0)
    d[1] <- sum(1/(x[1] + c(1:(sampleRichness-1))*x[2])) - sum(1/(x[1]+c(1:(sampleSize-1))))
    d[2] <- sum(c(1:(sampleRichness-1))/(x[1] + c(1:(sampleRichness-1))*x[2]))
    for(j in which(vSpeciesCount>1)) {
      d[2] <- d[2] + sum(1/(x[2]-c(1:(vSpeciesCount[j]-1))))
    }
    return(d)
  }
}

BayesNPSampleRichness = function(SampleSize,sampleRichness,sampleSize,Theta,Sigma) {
  if(SampleSize<0) {
    return(NA)
  } else {
    return(max((sampleRichness+Theta/Sigma)*(exp(fLogAscFact(Theta+Sigma+sampleSize,(SampleSize))-fLogAscFact(Theta+sampleSize,(SampleSize)))-1),0))
  }
}
 
BayesNP = function(tCountsPerSpecies, sampleSizes) {
  tCountsPerSpecies=tCountsPerSpecies[tCountsPerSpecies[,2]>0,]
  sampleRichness = dim(tCountsPerSpecies)[1]
  sampleSize = sum(tCountsPerSpecies[,2])
  if (sampleSize > max(sampleSizes)) {
    return(NULL)
  }
  control = list()
  control$fnscale = -1
  x_0 = c(0,0)
  x_1 = c(100,0.5)
  it = 0
  Optim = optim(par=x_1,fn=fLogPitman,gr=fGradLogPitman,vSpeciesCount=tCountsPerSpecies[,2],sampleRichness=sampleRichness,sampleSize=sampleSize,control=control)
  richest = sampleRichness+mapply(BayesNPSampleRichness,sampleSizes-sampleSize,sampleRichness,sampleSize,Optim$par[1],Optim$par[2])
  return(richest)
}
