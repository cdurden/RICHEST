setMethod("plot",signature="richest",where=".GlobalEnv",
  def=function(x,xlim=c(0,max(x@samplSizes)),ylim=c(0,max(x@estimates),na.rm=TRUE),...) {
    plot(x@sampleSizes,x@estimates,xlab='sample size',ylab='estimated sample richness',...)
  }
)
setMethod("points",signature="richest",where=".GlobalEnv",
  def=function(x,...) {
    points(x@sampleSizes,x@estimates,...)
  }
)
write.richest = function(x,...) {
  write.table(cbind(x@sampleSizes,x@estimates),...,sep="\t",quote=FALSE,row.names=FALSE)
}
