\name{funHDDC-package}
\alias{funHDDC-package}
\docType{package}
\title{
Model-based clustering in group-specific functional subspaces with the funHDDC algorithm
}
\description{
It provides the funHDDC algorithm (Bouveyron & Jacques, 2011) which allows to cluster functional data by modeling each group within a specific functional subspace. 
}
\details{
\tabular{ll}{
Package: \tab funHDDC\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-09-06\cr
License: \tab GPL-2\cr
}
}
\author{
Charles Bouveyron & Julien Jacques
Maintainer: Charles Bouveyron <charles.bouveyron@parisdescartes.fr>
}
\references{
C. Bouveyron & J. Jacques, Model-based Clustering of Time Series in Group-specific Functional Subspaces, Advances in Data Analysis and Classification, vol. 5 (4), pp. 281-300, 2011.
}
\keyword{ package }
%\seealso{}
\examples{
# Clustering the well-known "Canadian temperature" data (Ramsay & Silverman)
basis <- create.bspline.basis(c(0, 365), nbasis=21, norder=4)
fdobj <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"],basis,
  fdnames=list("Day", "Station", "Deg C"))$fd
res = funHDDC(fdobj,4,model='AkBQkDk',init='hclust',thd=0.001)

# Visualization of the partition and the group means
par(mfrow=c(1,2))
plot(fdobj,col=res$cls,lwd=2,lty=1)
fdmeans = fdobj; fdmeans$coefs = t(res$prms$m)
plot(fdmeans,col=1:max(res$cls),lwd=2)

## DO NOT RUN
# # Map of the results
# par(mfrow=c(1,12))
# library(maps)
# map("world", "canada")
# text(-CanadianWeather$coordinates[,2],CanadianWeather$coordinates[,1],
#      labels=rownames(CanadianWeather$coordinates),col=res$cls,cex=0.75)
}
