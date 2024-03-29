SampleRichnessVariance = function(p,vSampleSize) {
  S <- length(p)
  variances <- vector("numeric",length(vSampleSize))
  for (k in 1:length(vSampleSize)) {
    m <- vSampleSize[k]
    summand1 <- (1-p[1])^m*(1-(1-p[1])^m)
    summand2 <- 0
    for(i in 2:S) {
      summand1 <- summand1 + (1-p[i])^m*(1-(1-p[i])^m)
      summand2 <- summand2 + sum((1-p[i]-p[1:(i-1)])^m-((1-p[i])*(1-p[1:(i-1)]))^m)
    }
    variances[k] <- summand1 + summand2
  }
  return(cbind(vSampleSize,variances))
}

sampleRichnessApprox = function(p,sampleSizes) {
  sampleRichness = vector("numeric",length(sampleSizes))
  populationRichness = sum(p>0)
  for (i in seq(length(sampleSizes))) {
    sampleRichness[i] = populationRichness - sum((1-p)^(sampleSizes[i]))
  }
  return(sampleRichness)
}

sampleRichness = function(p,sampleSizes,populationSize) {
  sampleRichness = vector("numeric",length(sampleSizes))
  populationRichness = sum(p>0)
  for (i in seq(length(sampleSizes))) {
    if (sampleSizes[i] < populationSize) {
      sampleRichness[i] = sum(1-mapply(choose,populationSize-p*populationSize,sampleSizes[i])/mapply(choose,populationSize,sampleSizes[i]))
    } else {
      sampleRichness[i] = populationRichness
    }
  }
  return(as.matrix(cbind(sampleSizes,sampleRichness)))
}

RunEMSpeciesPerCount_C = function(tSpeciesPerCount,nSupportPoints,nSpecies) {
  tSpeciesPerCount = as.matrix(tSpeciesPerCount)
  nSpeciesObs = sum(tSpeciesPerCount[,2])
  ZeroIndex = which(tSpeciesPerCount[,1]==0)
  if(length(ZeroIndex)==0) {
    tSpeciesPerCount = rbind(tSpeciesPerCount,c(0,nSpecies-nSpeciesObs))
  } else {
    tSpeciesPerCount[ZeroIndex,2] = c(nSpecies-nSpeciesObs,rep(0,length(ZeroIndex)-1))
  }
  nCounts = dim(tSpeciesPerCount)[1]
  sSpeciesPerCount1 = paste("[",paste(tSpeciesPerCount[,1],collapse=", "),"]",sep="")
  sSpeciesPerCount2 = paste("[",paste(tSpeciesPerCount[,2],collapse=", "),"]",sep="")
  MLE = .Call("R_FindMLESpeciesPerCount",as.character(sSpeciesPerCount1),as.character(sSpeciesPerCount2),as.integer(nCounts),as.integer(nSupportPoints),as.integer(nSpecies),PACKAGE="richest")
  MixingF = data.frame(lambda=MLE[1:nSupportPoints],pi=MLE[(nSupportPoints+1):(2*nSupportPoints)])
  return(MixingF)
}

speciesProbabilities = function(estMixingDist,populationRichness) {
  p <- vector("numeric",populationRichness)
  quantiles <- c(0,cumsum(estMixingDist$pi))
  cat(quantiles)
  for(i in 1:populationRichness) {
    p[i] <- estMixingDist$lambda[min(max(which(quantiles<(i/populationRichness))),length(quantiles)-1)]
  }
  p <- p/sum(p)
  return(p)
}
