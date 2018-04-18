\name{funHDDC-package}
\alias{funHDDC-package}
\docType{package}
\title{
Clustering Univariate and Multivariate Functional Data in Group-specific Functional Subspaces
}
\description{
The funHDDC algorithm allows to cluster functional univariate and multivariate
data by modeling each cluster within a specific functional subspace. See Bouveyron and Jacques (2011) <doi:10.1007/s11634-011-0095-6> and Schmutz et al. (2017) <hal:01652467> for details on the model parameters described below.
}
\details{
\tabular{ll}{
Package: \tab funHDDC\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2018-04-18\cr
License: \tab GPL-2\cr
}
}
\author{
A Schmutz, J. Jacques & C. Bouveyron
Maintainer: Charles Bouveyron <charles.bouveyron@gmail.com>
}
\references{
- C. Bouveyron and J. Jacques, Model-based Clustering of Time Series in Group-specific Functional Subspaces, Advances in Data Analysis and Classification, vol. 5 (4), pp. 281-300, 2011.
- A. Schmutz, C. Bouveyron, J. Jacques, L. Cheze and  P. Martin, Clustering multivariate functional data in group-specic functional subspaces, Preprint HAL 01652467, Université Côte d'Azur, 2017.
}
\keyword{ package }
%\seealso{}
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

res.multi<-funHDDC(list(var1,var2),K=3,model="AkBkQkDk",init="kmeans",threshold=0.2)
par(mfrow=c(1,2))
plot(var1,col=res.multi$class)
plot(var2,col=res.multi$class)
}
