predict.funHDDC<-function(object,newdata,...){
  model=object
  if (class(newdata)!='list') {x = t(newdata$coefs)
  }else {x = t(newdata[[1]]$coefs); for (i in 2:length(newdata)) x = cbind(x,t(newdata[[i]]$coefs))}
  p <- ncol(x)
  N <- nrow(x)
  K <- model$K
  a <- model$a
  b <- model$b
  mu <- model$mu
  d <- model$d
  prop <- model$prop
  Q <- model$Q

  b[b<1e-6] <- 1e-6

  if (model$model%in%("AKJBKQKDK")){
    K_pen <- matrix(0,K,N)
    for (i in 1:K) {
      s <- sum(log(a[i,1:d[i]]))
      X <- x - matrix(mu[i,], N, p, byrow=TRUE)
      Qi = model$fpca[[i]]$W %*% Q[[i]]
      proj <- (X%*%Qi)%*%t(Qi)
      A <- (-proj)%*%Qi%*%sqrt(diag(1/a[i,1:d[i]],d[i]))
      B <- X-proj
      K_pen[i,] <- rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])+p*log(2*pi)
    }
  }else if (model$model%in%("AKJBQKDK")){
    K_pen <- matrix(0,K,N)
    for (i in 1:K) {
    s <- sum(log(a[i,1:d[i]]))
    X <- x - matrix(mu[i,], N, p, byrow=TRUE)
    Qi = model$fpca[[i]]$W %*% Q[[i]]
    proj <- (X%*%Qi)%*%t(Qi)
    A <- (-proj)%*%Qi%*%sqrt(diag(1/a[i,1:d[i]],d[i]))
    B <- X-proj
    K_pen[i,] <- rowSums(A^2)+1/b[1]*rowSums(B^2)+s+(p-d[i])*log(b[1])-2*log(prop[i])+p*log(2*pi)
    }
  }else if (model$model%in%("ABQKDK")){
    K_pen <- matrix(0,K,N)
    for (i in 1:K) {
      s <- sum(log(a[1]))
      X <- x - matrix(mu[i,], N, p, byrow=TRUE)
      Qi = model$fpca[[i]]$W %*% Q[[i]]
      proj <- (X%*%Qi)%*%t(Qi)
      A <- (-proj)%*%Qi%*%sqrt(diag(1/a[1],d[i]))
      B <- X-proj
      K_pen[i,] <- rowSums(A^2)+1/b[1]*rowSums(B^2)+s+(p-d[i])*log(b[1])-2*log(prop[i])+p*log(2*pi)
    }
  } else if (model$model%in%("AKBKQKDK")){
    K_pen <- matrix(0,K,N)
    for (i in 1:K) {
      s <- sum(log(a[1,i]))
      X <- x - matrix(mu[i,], N, p, byrow=TRUE)
      Qi = model$fpca[[i]]$W %*% Q[[i]]
      proj <- (X%*%Qi)%*%t(Qi)
      A <- (-proj)%*%Qi%*%sqrt(diag(1/a[1,i],d[i]))
      B <- X-proj
      K_pen[i,] <- rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])+p*log(2*pi)
    }
  }else if(model$model%in%("ABKQKDK")){
    K_pen <- matrix(0,K,N)
    for (i in 1:K) {
      s <- sum(log(a[1]))
      X <- x - matrix(mu[i,], N, p, byrow=TRUE)
      Qi = model$fpca[[i]]$W %*% Q[[i]]
      proj <- (X%*%Qi)%*%t(Qi)
      A <- (-proj)%*%Qi%*%sqrt(diag(1/a[1],d[i]))
      B <- X-proj
      K_pen[i,] <- rowSums(A^2)+1/b[i]*rowSums(B^2)+s+(p-d[i])*log(b[i])-2*log(prop[i])+p*log(2*pi)
    }
  }else if (model$model%in%("AKBQKDK")){
    K_pen <- matrix(0,K,N)
    for (i in 1:K) {
      s <- sum(log(a[1,i]))
      X <- x - matrix(mu[i,], N, p, byrow=TRUE)
      Qi = model$fpca[[i]]$W %*% Q[[i]]
      proj <- (X%*%Qi)%*%t(Qi)
      A <- (-proj)%*%Qi%*%sqrt(diag(1/a[1,i],d[i]))
      B <- X-proj
      K_pen[i,] <- rowSums(A^2)+1/b[1]*rowSums(B^2)+s+(p-d[i])*log(b[1])-2*log(prop[i])+p*log(2*pi)
    }
  }
  A <- -1/2*t(K_pen)
  L <- sum(log(rowSums(exp(A-apply(A,1,max))))+apply(A,1,max))

  t <- matrix(0,N,K)
  for (i in 1:K) t[,i] <- 1/rowSums(exp((K_pen[i,]-t(K_pen))/2))
  class<-matrix(nrow=N,ncol=1)
  for (l in 1:N) class[l,]<-which.max(t[l,])
  list(class=class,t=t,L=L)
}

