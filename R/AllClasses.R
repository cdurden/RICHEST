setClass("richestInput",representation(x="data.frame"))
setClass("richestSPC",contains="richestInput")
setClass("richestCPS",contains="richestInput")

setMethod("initialize",signature="richestInput",
def=function(.Object,x) {
  .Object@x = x
  return(.Object)
})

setClass("richest", representation(estimates="numeric",sampleSizes="numeric",supports="numeric",pop="numeric",mix="data.frame",likelihood="numeric"))

setMethod("initialize",signature="richest",
def=function(.Object,sampleSizes,supports=0,pop=0) { 
  .Object@supports = supports
  .Object@pop = pop
  .Object@sampleSizes = sampleSizes
  return(.Object)
})

