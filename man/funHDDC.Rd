\name{funHDDC}
\alias{funHDDC}
\title{Clustering univariate and multivariate functional data in group-specific functional subspaces}
\usage{
funHDDC(data,K=1:10,model="AkjBkQkDk",threshold=0.2,itermax=200,eps=1e-6,init="kmeans",
criterion="bic",algo='EM', d_select="Cattell", init.vector=NULL, show=TRUE,
mini.nb=c(5, 10),min.individuals=2, mc.cores=1, nb.rep=1, keepAllRes=TRUE,
kmeans.control = list(), d_max=100)
}
\description{
The funHDDC algorithm  allows to cluster functional univariate (Bouveyron and Jacques, 2011, <doi:10.1007/s11634-011-0095-6>) or multivariate data (Schmutz et al., 2018) by modeling each group within a specific functional subspace.
}
\arguments{
    \item{data}{In the univariate case: a functional data object produced by the fda package, in the multivariate     case: a list of functional data objects.}
    \item{K}{The number of clusters, you can test one partition or multiple partitions at the same time, for example K=2:10.}
    \item{model}{The chosen model among 'AkjBkQkDk', 'AkjBQkDk', 'AkBkQkDk', 'ABkQkDk', 'AkBQkDk', 'ABQkDk'.          'AkjBkQkDk' is the default. You can test multiple models at the same time with the command c(). For example c("AkjBkQkDk","AkjBQkDk").}
    \item{threshold}{The threshold of the Cattell' scree-test used for selecting the group-specific intrinsic         dimensions.}
    \item{itermax}{The maximum number of iterations.}
    \item{eps}{The threshold of the convergence criterion.}
    \item{init}{A character string. It is the way to initialize the E-M algorithm. There are five ways of
    initialization: “kmeans” (default), “param”, “random”, “mini-em” or “vector”. See details for more information.
    It can also be directly initialized with a vector containing the prior classes of the observations.}
    \item{criterion}{The criterion used for model selection: bic (default) or icl. But we recommand using              slopeHeuristic function provided in this package.}
    \item{algo}{A character string indicating the algorithm to be used. The available algorithms are the
    Expectation-Maximisation ("EM"), the Classification E-M ("CEM") and the Stochastic E-M ("SEM"). The default        algorithm is the "EM".}
    \item{d_select}{Either “Cattell” (default) or “BIC”. This parameter selects which method to use to select the       intrinsic dimensions of subgroups.}
    \item{init.vector}{A vector of integers or factors. It is a user-given initialization. It should be of the same     length as of the data. Only used when init="vector".}
    \item{show}{Use show = FALSE to settle off the informations that may be printed.}
    \item{mini.nb}{A vector of integers of length two. This parameter is used in the “mini-em” initialization. The     first integer sets how many times the algorithm is repeated; the second sets the maximum number of iterations      the algorithm will do each time. For example, if init=“mini-em” and mini.nb=c(5,10), the algorithm wil be          launched 5 times, doing each time 10 iterations; finally the algorithm will begin with the initialization that     maximizes the log-likelihood.}
    \item{min.individuals}{This parameter is used to control for the minimum population of a class. If the             population of a class becomes stricly inferior to 'min.individuals' then the algorithm stops and gives the         message: 'pop<min.indiv.'. Here the meaning of "population of a class" is the sum of its posterior
    probabilities. The value of 'min.individuals' cannot be lower than 2.}
    \item{mc.cores}{Positive integer, default is 1. If mc.cores>1, then parallel computing is used, using mc.cores     cores. Warning for Windows users only: the parallel computing can sometimes be slower than using one single        core (due to how parLapply works).}
    \item{nb.rep}{A positive integer (default is 1 for kmeans initialization and 20 for random initialization).        Each estimation (i.e. combination of (model, K, threshold)) is repeated nb.rep times and only the estimation       with the highest log-likelihood is kept.}
    \item{keepAllRes}{Logical. Should the results of all runs be kept? If so, an argument all_results is created in     the results. Default is TRUE.}
    \item{kmeans.control}{A list. The elements of this list should match the parameters of the kmeans
    initialization (see kmeans help for details). The parameters are “iter.max”, “nstart” and “algorithm”.}
    \item{d_max}{A positive integer. The maximum number of dimensions to be computed. Default is 100. It means that     the instrinsic dimension of any cluster cannot be larger than d_max. It quickens a lot the algorithm for            datasets with a large number of variables (e.g. thousands).}

}
\value{
  \item{d}{The number of dimensions for each cluster.}
  \item{a}{Values of parameter a for each cluster.}
  \item{b}{Values of parameter b for each cluster.}
  \item{mu}{The mean of each cluster in the original space.}
  \item{prop}{The proportion of individuals in each cluster.}
  \item{loglik}{The maximum of log-likelihood.}
  \item{loglik_all}{The log-likelihood at each iteration.}
  \item{posterior}{The posterior probability for each individual to belong to each cluster.}
  \item{class}{The clustering partition.}
  \item{BIC}{The BIC value.}
  \item{ICL}{The ICL value.}
  \item{complexity}{the number of parameters that are estimated.}
  \item{all_results}{if multiple number of clusters are tested or multiple models, results of each models are stored here}

}
\details{
If you choose init="random", the algorithm is run 20 times with the same model options and the solution which maximises the log-likelihood is printed. This explains why sometimes with this initialization it runs a bit slower than with 'kmeans' initialization.

If the warning message: "In funHDDC(...) : All models diverged" is printed, it means that the algorithm found less classes that the number you choose (parameter K). Because of the EM algorithm, it should be because of a bad initialization of the EM algorithm. So you have to restart the algorithm multiple times in order to check if, with a new initialization of the EM algorithm the model converges, or if there is no solution with the number K you choose.

The different initializations are:
\describe{
\item{“param”:}{it is initialized with the parameters, the means being generated by a multivariate normal distribution and the covariance matrix being common to the whole sample.}
\item{“mini-em”:}{it is an initialization strategy, the classes are randomly initialized and the E-M algorithm makes several iterations, this action is repetead a few times (the default is 5 iterations and 10 times), at the end, the initialization choosen is the one which maximise the log-likelihood (see mini.nb for more information about its parametrization).}
\item{“random”:}{the classes are randomly given using a multinomial distribution}
\item{“kmeans”:}{the classes are initialized using the kmeans function (with: algorithm="Hartigan-Wong"; nstart=4; iter.max=50); note that the user can use his own arguments for kmeans using the dot-dot-dot argument}
\item{A prior class vector:}{It can also be directly initialized with a vector containing the prior classes of the observations. To do so use init="vector" and provide the vector in the argument init.vector.}
}

Note that the BIC criterion used in this function is to be maximized and is defined as 2*LL-k*log(n) where LL is the log-likelihood, k is the number of parameters and n is the number of observations.

}

\references{
- C.Bouveyron and J.Jacques, Model-based Clustering of Time Series in Group-specific Functional Subspaces, Advances in Data Analysis and Classification, vol. 5 (4), pp. 281-300, 2011 <doi:10.1007/s11634-011-0095-6>

- C. Bouveyron, L. Cheze, J. Jacques, P. Martin and A. Schmutz, Clustering multivariate functional data in group-specic functional subspaces, Preprint HAL 01652467, 2017.
}

\examples{
##Univariate example
library(fda)
data("trigo")
basis<- create.bspline.basis(c(0,1), nbasis=25)
var1<-smooth.basis(argvals=seq(0,1,length.out = 100),y=t(trigo[,1:100]),fdParobj=basis)$fd

res.uni<-funHDDC(var1,K=2,model="AkBkQkDk",init="kmeans",threshold=0.2)
plot(var1,col=res.uni$class)


##Multivariate example
library(fda)
data("triangle")
basis<- create.bspline.basis(c(1,21), nbasis=25)
var1<-smooth.basis(argvals=seq(1,21,length.out = 101),y=t(triangle[,1:101]),fdParobj=basis)$fd
var2<-smooth.basis(argvals=seq(1,21,length.out = 101),y=t(triangle[,102:202]),fdParobj=basis)$fd

res.multi<-funHDDC(list(var1,var2),K=3,model="AkjBkQkDk",init="kmeans",threshold=0.2)
par(mfrow=c(1,2))
plot(var1,col=res.multi$class)
plot(var2,col=res.multi$class)


##NOT RUN:
##You can test multiple number of groups at the same time
#data("trigo")
#basis<- create.bspline.basis(c(0,1), nbasis=25)
#var1<-smooth.basis(argvals=seq(0,1,length.out = 100),y=t(trigo[,1:100]),fdParobj=basis)$fd
#res.uni<-funHDDC(var1,K=2:4,model="AkBkQkDk",init="kmeans",threshold=0.2)

##You can test multiple models at the same time
#data("trigo")
#basis<- create.bspline.basis(c(0,1), nbasis=25)
#var1<-smooth.basis(argvals=seq(0,1,length.out = 100),y=t(trigo[,1:100]),fdParobj=basis)$fd
#res.uni<-funHDDC(var1,K=3,model=c("AkjBkQkDk","AkBkQkDk"),init="kmeans",threshold=0.2)


##another example on Canadian data
#library(fda)

###Clustering the "Canadian temperature" data (Ramsey & Silverman): univariate case
#daybasis65 <- create.fourier.basis(c(0, 365), nbasis=65, period=365)
#daytempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"], daybasis65,
#    fdnames=list("Day", "Station", "Deg C"))$fd

#res.uni<-funHDDC(daytempfd,K=3,model="AkjBkQkDk",init="random",threshold=0.2)

##Graphical representation of groups mean curves
#select1<-fd(daytempfd$coefs[,which(res.uni$class==1)],daytempfd$basis)
#select2<-fd(daytempfd$coefs[,which(res.uni$class==2)],daytempfd$basis)
#select3<-fd(daytempfd$coefs[,which(res.uni$class==3)],daytempfd$basis)

#plot(mean.fd(select1),col="lightblue",ylim=c(-20,20),lty=1,lwd=3)
#lines(mean.fd(select2),col="palegreen2",lty=1,lwd=3)
#lines(mean.fd(select3),col="navy",lty=1,lwd=3)


###Clustering the "Canadian temperature" data (Ramsey & Silverman): multivariate case
#daybasis65 <- create.fourier.basis(c(0, 365), nbasis=65, period=365)
#daytempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"], daybasis65,
#    fdnames=list("Day", "Station", "Deg C"))$fd
#dayprecfd<-smooth.basis(day.5, CanadianWeather$dailyAv[,,"Precipitation.mm"], daybasis65,
#    fdnames=list("Day", "Station", "Mm"))$fd

#res.multi<-funHDDC(list(daytempfd,dayprecfd),K=3,model="AkjBkQkDk",init="random",threshold=0.2)

##Graphical representation of groups mean curves
##Temperature
#select1<-fd(daytempfd$coefs[,which(res.multi$class==1)],daytempfd$basis)
#select2<-fd(daytempfd$coefs[,which(res.multi$class==2)],daytempfd$basis)
#select3<-fd(daytempfd$coefs[,which(res.multi$class==3)],daytempfd$basis)

#plot(mean.fd(select1),col="lightblue",ylim=c(-20,20),lty=1,lwd=3)
#lines(mean.fd(select2),col="palegreen2",lty=1,lwd=3)
#lines(mean.fd(select3),col="navy",lty=1,lwd=3)

##Precipitation
#select1b<-fd(dayprecfd$coefs[,which(res.multi$class==1)],dayprecfd$basis)
#select2b<-fd(dayprecfd$coefs[,which(res.multi$class==2)],dayprecfd$basis)
#select3b<-fd(dayprecfd$coefs[,which(res.multi$class==3)],dayprecfd$basis)

#plot(mean.fd(select1b),col="lightblue",ylim=c(0,5),lty=1,lwd=3)
#lines(mean.fd(select2b),col="palegreen2",lty=1,lwd=3)
#lines(mean.fd(select3b),col="navy",lty=1,lwd=3)


}

