gtRichest = function(speciesPerCount, sampleSizes) {
  sampleSize = sum(apply(as(speciesPerCount,"matrix"),1,prod))
  sampleRichness = sum(speciesPerCount[,2])
  if (sampleSize > max(sampleSizes)) {
    return(NULL)
  }
  richest = rep(NA,length(sampleSizes))
  names(richest) = sampleSizes
  for(i in which(sampleSizes>sampleSize)) {
    richest[i] = sampleRichness + sum((-1)^(speciesPerCount[,1]+1)*((sampleSizes[i]-sampleSize)/sampleSize)^speciesPerCount[,1]*speciesPerCount[,2])
  }
  return(richest)
}
