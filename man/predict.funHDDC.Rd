\name{predict.funHDDC}
\alias{predict.funHDDC}
\title{Predict new results from a funHDDC model}
\usage{
  \method{predict}{funHDDC}(object,newdata,...)
}
\description{
  It provides predictions for a new functional dataset using the funHDDC model you train.
}
\arguments{
  \item{object}{An object produced by funHDDC function}
  \item{newdata}{In the univariate case: a functional data object produced by the fda package, in the multivariate     case: a list of functional data objects, for which you want to predict the group membership from the model.}
  \item{...}{Some additional parameters.}
  }
    \value{
    \item{class}{The clustering partition.}
    \item{t}{The probability of each individual to belong to each cluster.}
    \item{L}{The loglikelihood.}

    }
\examples{
##Univariate example
    data("trigo")
    basis<- create.bspline.basis(c(0,1), nbasis=25)
    suppressWarnings(RNGversion("3.5.0"))
    set.seed(1)
    z<-sample(1:100,0.8*100,replace=FALSE)
    var1_train<-smooth.basis(argvals=seq(0,1,length.out = 100),y=t(trigo[z,1:100]),
                             fdParobj=basis)$fd
    var1_test<-smooth.basis(argvals=seq(0,1,length.out = 100),y=t(trigo[-z,1:100]),
                            fdParobj=basis)$fd

    model<-funHDDC(var1_train,K=2)
    pred<-predict(model,var1_test)

##Multivariate example
    data("triangle")
    basis<- create.bspline.basis(c(1,21), nbasis=25)
    suppressWarnings(RNGversion("3.5.0"))
    set.seed(1)
    z<-sample(1:100,0.8*100,replace=FALSE)
    var1_train<-smooth.basis(argvals=seq(1,21,length.out = 101),y=t(triangle[z,1:101]),
                              fdParobj=basis)$fd
    var1_test<-smooth.basis(argvals=seq(1,21,length.out = 101),y=t(triangle[-z,1:101]),
                            fdParobj=basis)$fd
    var2_train<-smooth.basis(argvals=seq(1,21,length.out = 101),y=t(triangle[z,102:202]),
                            fdParobj=basis)$fd
    var2_test<-smooth.basis(argvals=seq(1,21,length.out = 101),y=t(triangle[-z,102:202]),
                            fdParobj=basis)$fd

    model<-funHDDC(list(var1_train,var2_train),K=3)
    pred<-predict(model,list(var1_test,var2_test))

##NOT RUN: another example on Canadian data
    ##Clustering the "Canadian temperature" data (Ramsey & Silverman): univariate case
    daybasis65 <- create.fourier.basis(c(0, 365), nbasis=65, period=365)
    #z<-sample(1:35,0.9*35,replace=FALSE)
    #daytempfd_train <- smooth.basis(day.5, CanadianWeather$dailyAv[,z,"Temperature.C"],
    #daybasis65,fdnames=list("Day", "Station", "Deg C"))$fd
    #daytempfd_test <- smooth.basis(day.5, CanadianWeather$dailyAv[,-z,"Temperature.C"],
    #daybasis65,fdnames=list("Day", "Station", "Deg C"))$fd

     #model<-funHDDC(daytempfd_train,K=3)
     #pred<-predict.funHDDC(model,daytempfd_test)


    ##Clustering the "Canadian temperature" data (Ramsey & Silverman): multivariate case
    #daybasis65 <- create.fourier.basis(c(0, 365), nbasis=65, period=365)
    #z<-sample(1:35,0.9*35,replace=FALSE)
    #daytempfd_train <- smooth.basis(day.5, CanadianWeather$dailyAv[,z,"Temperature.C"],
    #daybasis65,fdnames=list("Day", "Station", "Deg C"))$fd
    #daytempfd_test <- smooth.basis(day.5, CanadianWeather$dailyAv[,-z,"Temperature.C"],
    #daybasis65,fdnames=list("Day", "Station", "Deg C"))$fd
    #dayprecfd_train<-smooth.basis(day.5, CanadianWeather$dailyAv[,z,"Precipitation.mm"],
    #daybasis65,fdnames=list("Day", "Station", "Mm"))$fd
    #dayprecfd_test<-smooth.basis(day.5, CanadianWeather$dailyAv[,-z,"Precipitation.mm"],
    #daybasis65,fdnames=list("Day", "Station", "Mm"))$fd

    #model<-funHDDC(list(daytempfd_train,dayprecfd_train),K=3)
    #pred<-predict.funHDDC(model,list(daytempfd_test,dayprecfd_test))

    }
