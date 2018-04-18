plot.mfpca<-function(x,nharm=3,threshold=0.05,...){
  mfpcaobj<-x
  if (nharm>dim(mfpcaobj$harmonics$coefs)[[2]]){
    warning("The number of harmonics you choose is superior of the one of your mfpca")
  }

  if(length(mfpcaobj$call[[2]])==1){
  #Data plot
  v<-as.character(mfpcaobj$call['fdobj'])
  data_obj<-get(v)
  m<-paste(v,"plot",sep=" ")
  plot(data_obj,main=m)

  #Scores plot
  for (j in 1:nharm){
    if (mfpcaobj$varprop[j]<threshold) break
  xax<-paste('Dimension',j,sep=" ")
  yaxis<-paste('Dimension',j+1,sep=" ")
  plot(mfpcaobj$scores[,j],mfpcaobj$scores[,j+1],main="Scores plot",xlab = xax,ylab=yaxis)
  text(x=mfpcaobj$scores[,j],y=mfpcaobj$scores[,j+1], labels = c(1:nrow(mfpcaobj$scores)), cex=0.7,pos=2)
  }

  #Definition of elements needed for the 2 graphics
  range<-mfpcaobj$harmonics$basis$rangeval
  nbas<-mfpcaobj$harmonics$basis$nbasis

  #Mean variations
  t=seq(range[1],range[2],length=nbas)
  mean_function_values<-cbind(eval.fd(t,mfpcaobj$meanfd),t)
  fdmat<-eval.fd(t,mfpcaobj$harmonic)
  for (i in 1:nharm){
  if (mfpcaobj$varprop[i]<threshold) break
  uplim<-cbind(mean_function_values[,1]+sqrt(mfpcaobj$eigval[i])*fdmat[,i],t)
  lowlim<-cbind(mean_function_values[,1]-sqrt(mfpcaobj$eigval[i])*fdmat[,i],t)

  if(min(uplim[,1])<min(lowlim[,1])){
   downy<-min(uplim[,1])-1
  }else{
   downy<-min(lowlim[,1])-1
  }
  if(max(uplim[,1])<max(lowlim[,1])){
  upy<-max(lowlim[,1])+1
  }else{
  upy<-max(uplim[,1])+1
  }

  xax<-paste('Variation of the mean curve, harmonic',i,sep=" ")
  plot(mfpcaobj$meanfd,ylim=c(downy,upy),main=xax)
  lines(x=uplim[,2],y=uplim[,1],col='red',cex=0.6,lty=2)
  lines(x=lowlim[,2],y=lowlim[,1],col='blue',cex=0.6,lty=2)
  }

  #Plot of eigenfunctions
  t=seq(range[1],range[2],length=nbas+600)
  fdmat<-eval.fd(t,mfpcaobj$harmonics)
  for (k in 1:nharm){
  if (mfpcaobj$varprop[k]<threshold) break
  xax<-paste('Dimension',k, ", Proportion of variance:",round(mfpcaobj$varprop[k],2) ,sep=" ")
  plot(t,fdmat[,k],main=xax,xlab="x",ylab="y",type="l")
  }

  }else if (length(mfpcaobj$call[[2]])>1){

    ##Data plot
    v<-as.character(mfpcaobj$call[[2]])
    d<-length(v)
    for (n in 2:d){
      data_obj<-get(v[[n]])
      m<-paste(v[[n]],"plot",sep=" ")
      plot(data_obj,main=m)

    }

    ##Scores plot
    for (j in 1:nharm){
      if (mfpcaobj$varprop[j]<threshold) break
      xax<-paste('Dimension',j,sep=" ")
      yaxis<-paste('Dimension',j+1,sep=" ")
      plot(mfpcaobj$scores[,j],mfpcaobj$scores[,j+1],main="Scores plot",xlab = xax,ylab=yaxis)
      text(x=mfpcaobj$scores[,j],y=mfpcaobj$scores[,j+1], labels = c(1:nrow(mfpcaobj$scores)), cex=0.7,pos=2)
    }

    #Definition of elements needed for the 2 graphics
    range<-mfpcaobj$harmonics$basis$rangeval
    nbas<-mfpcaobj$harmonics$basis$nbasis
    nvar<-length(as.character(mfpcaobj$call[[2]]))-1

    ##Mean variations
    t=seq(range[1],range[2],length=nbas)
    if (mfpcaobj$harmonics$basis[[2]]=="fourier"){
      new.basis<-create.fourier.basis(c(mfpcaobj$harmonics$basis$rangeval),nbasis=nbas)
    }else{
      new.basis<-create.bspline.basis(c(mfpcaobj$harmonics$basis$rangeval),nbasis=nbas)
    }

    #Cas de la première variable
    mean_function_values<-cbind(eval.fd(t,fd(mfpcaobj$meanfd[[1]]$coefs,mfpcaobj$meanfd[[1]]$basis)),t)
    harmo<-mfpcaobj$harmonics
    harmo$basis<-new.basis
    harmo$coefs<-mfpcaobj$harmonics$coefs[1:nbas,1:nbas]
    fdmat<-eval.fd(t,harmo)
    for (i in 1:nharm){
      if (mfpcaobj$varprop[i]<threshold) break
      uplim<-cbind(mean_function_values[,1]+sqrt(mfpcaobj$eigval[i])*fdmat[,i],t)
      lowlim<-cbind(mean_function_values[,1]-sqrt(mfpcaobj$eigval[i])*fdmat[,i],t)

      if(min(uplim[,1])<min(lowlim[,1])){
        downy<-min(uplim[,1])-1
      }else{
        downy<-min(lowlim[,1])-1
      }
      if(max(uplim[,1])<max(lowlim[,1])){
        upy<-max(lowlim[,1])+1
      }else{
        upy<-max(uplim[,1])+1
      }
      xax<-paste('Variation of the mean curve, Variable 1, Dimension',i,sep=" ")
      plot(mfpcaobj$meanfd[[1]],ylim=c(downy,upy),main=xax)
      lines(x=uplim[,2],y=uplim[,1],col='red',cex=0.6,lty=2)
      lines(x=lowlim[,2],y=lowlim[,1],col='blue',cex=0.6,lty=2)
    }

    #Cas des autres variables
    for (m in 1:(nvar-1)){
    mean_function_values<-cbind(eval.fd(t,fd(mfpcaobj$meanfd[[(m+1)]]$coefs,mfpcaobj$meanfd[[(m+1)]]$basis)),t)
    harmo2<-mfpcaobj$harmonics
    harmo2$basis<-new.basis
    harmo2$coefs<-mfpcaobj$harmonics$coefs[((m*nbas)+1):((m+1)*nbas),((m*nbas)+1):((m+1)*nbas)]
    fdmat<-eval.fd(t,harmo2)
    for (p in 1:nharm){
      if (mfpcaobj$varprop[p]<threshold) break
      uplim<-cbind(mean_function_values[,1]+sqrt(mfpcaobj$eigval[p])*fdmat[,p],t)
      lowlim<-cbind(mean_function_values[,1]-sqrt(mfpcaobj$eigval[p])*fdmat[,p],t)

      if(min(uplim[,1])<min(lowlim[,1])){
        downy<-min(uplim[,1])-1
      }else{
        downy<-min(lowlim[,1])-1
      }
      if(max(uplim[,1])<max(lowlim[,1])){
        upy<-max(lowlim[,1])+1
      }else{
        upy<-max(uplim[,1])+1
      }
      xax<-paste('Variation of the mean curve, Variable', (m+1), ' Dimension',p,sep=" ")
      plot(mfpcaobj$meanfd[[(m+1)]],ylim=c(downy,upy),main=xax)
      lines(x=uplim[,2],y=uplim[,1],col='red',cex=0.6,lty=2)
      lines(x=lowlim[,2],y=lowlim[,1],col='blue',cex=0.6,lty=2)
    }
    }

    ##Plot of eigenfunctions
    t=seq(range[1],range[2],length=nbas+600)
    if (mfpcaobj$harmonics$basis[[2]]=="fourier"){
      new.basis<-create.fourier.basis(c(mfpcaobj$harmonics$basis$rangeval),nbasis=nbas)
    }else{
      new.basis<-create.bspline.basis(c(mfpcaobj$harmonics$basis$rangeval),nbasis=nbas)
    }

    #Cas de la première variable
    harmo<-mfpcaobj$harmonics
    harmo$basis<-new.basis
    harmo$coefs<-mfpcaobj$harmonics$coefs[1:nbas,1:nbas]
    fdmat<-eval.fd(t,harmo)

    for (k in 1:nharm){
      if (mfpcaobj$varprop[k]<threshold) break
      xax<-paste('Eigenfunction plot, Dimension',k, ", Proportion of variance:",round(mfpcaobj$varprop[k],2),", Variable 1",sep=" ")
      plot(t,fdmat[,k],main=xax,xlab="x",ylab="y",type="l")
    }

    #Cas des variables suivantes
    for (l in 1:(nvar-1)){
    harmo<-mfpcaobj$harmonics
    harmo$basis<-new.basis
    harmo$coefs<-mfpcaobj$harmonics$coefs[((l*nbas)+1):((l+1)*nbas),((l*nbas)+1):((l+1)*nbas)]
    fdmat<-eval.fd(t,harmo)

    for (q in 1:nharm){
      if (mfpcaobj$varprop[q]<threshold) break
      xax<-paste('Eigenfunction plot, Dimension',q, ", Proportion of variance:",round(mfpcaobj$varprop[q],2),", Variable",(l+1) ,sep=" ")
      plot(t,fdmat[,q],main=xax,xlab="x",ylab="y",type="l")
    }

    }
  }
}

