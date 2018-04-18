slopeHeuristic <- function(mod) {
  fin<-length(mod$all_results)
  K=c((mod$all_results[[1]]$K):(mod$all_results[[1]]$K+fin-1))
  nbparam=mod$all_results[[1]]$complexity
  loglik=mod$all_results[[1]]$loglik
  for (i in 2:fin){
  nbparam=cbind(nbparam,mod$all_results[[i]]$complexity)
  loglik=cbind(loglik,mod$all_results[[i]]$loglik)
  }
  #library(MASS)
  dd=data.frame(nbp=c(nbparam),ll=c(loglik),K=c(K))
  fit=MASS::rlm(ll ~ nbp,data=dd,method='MM')
  if(fit$coefficients[2]<0) fit$coefficients[2]=0
  dd$llpen = dd$ll- 2* fit$coefficients[2]*dd$nbp
  par(mfrow=c(1,2))
  plot(dd$nbp,dd$ll,type='b',xlab='Model dimension',ylab='Log-likelihood')
  abline(fit,col='red')
  plot(dd$K,dd$llpen,type='p',xlab='K',ylab='Penalized log-likelihood')
  points(dd$K[which.max(dd$llpen)],max(dd$llpen),pch=19,col='red')
  return(K[which.max(dd$llpen)])
}
