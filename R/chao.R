chao <- function(cps) {
  chaoLB <- length(cps) + sum(cps==1)^2/(2*sum(cps==2))
  return(chaoLB)
}
