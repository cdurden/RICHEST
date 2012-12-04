setMethod("richest",signature=c("richestInput","numeric","character"),
def = function(x,sampleSizes,method,...) {
  spc = slot(as(x,"richestSPC"),"x") 
  if(any(sampleSizes<sum(spc[,1]*spc[,2]))) {
    warning("At least one sample sizes given for richness estimation is less than the size of the given sample")
  }
  richest = new("richest",sampleSizes,...)
  if(method=="GT") {
    x = as(x,"richestSPC")
    richest@estimates = gtRichest(x@x,richest@sampleSizes)
  }
  if(method=="Lijou") {
    x = as(x,"richestCPS")
    richest@estimates = BayesNP(x@x,richest@sampleSizes)
  }
  if(method=="PNPMLE") {
    if(richest@pop==0) { richest@pop = chao(as(x,"richestCPS")@x) }
    x = as(x,"richestSPC")
    if(richest@supports==0 || richest@supports>nrow(x@x)) richest@supports = nrow(x@x)
    richest@mix = RunEMSpeciesPerCount_C(x@x,richest@supports,richest@pop)
    p = speciesProbabilities(richest@mix,richest@pop)
    richest@estimates = sampleRichnessApprox(p,richest@sampleSizes)
  }
  return(richest)
})
