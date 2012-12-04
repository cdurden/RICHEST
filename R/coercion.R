setAs("richestCPS","richestSPC",
function(from) {
  spc = NULL
  for(i in 1:length(from@x[[1]])) {
    if(!any(spc[,1]==from@x[[2]][i])) {
      spc <- rbind(spc, c(from@x[[2]][i],1))
    } else {
      j <- which(spc[,1]==from@x[[2]][i])
      spc[j,2] <- spc[j,2]+1
    }
  }
  return(new("richestSPC",data.frame(spc)))
})
setAs("richestSPC","richestCPS",
function(from) {
  counts <- NULL
  for(i in 1:length(from@x[[1]])) {
    counts <- c(counts,rep(from@x[[1]][i],from@x[[2]][i]))
  }
  return(new("richestCPS",as.data.frame(cbind(c(1:length(counts)),counts))))
})
