\name{funHDDC}
\alias{funHDDC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model-based clustering in group-specific functional subspaces
}
\description{
It provides the funHDDC algorithm (Bouveyron & Jacques, 2011) which allows to cluster functional data by modeling each group within a specific functional subspace. 
}
\usage{
funHDDC(fd, K, init = "hclust", model = "AkBkQkDk", thd = 0.05, maxit = 50,
	 eps = 1e-6, ...)
}
\arguments{
  \item{fd}{
a functional data object produced by the fda package.
}
  \item{K}{
the number of clusters.
}
  \item{init}{
the initialization type ('random', 'kmeans' of 'hclust'). 'hclust' is the default.
}
  \item{model}{
the chosen model among 'AkjBkQkDk', 'AkjBQkDk','AkBkQkDk','AkBQkDk','ABkQkDk','ABQkDk'. See (Bouveyron & Jacques, 2011) for details.
}
  \item{thd}{
the threshold of the Cattell' scree-test used for selecting the group-specific initrisic dimensions. See (Bouveyron & Jacques, 2011) for details.
}
  \item{maxit}{
the maximum number of iterations.
}
  \item{eps}{
the threshold of the convergence criterion.
}
  \item{\dots}{
additional options for internal functions.
}
}
\value{
\item{cls}{the clustering partition}
\item{P}{the posterior probablities}
\item{prms}{the model parameters}
\item{bic}{the BIC value}
\item{aic}{the AIC value}
\item{icl}{the ICL value}
\item{loglik}{vector of the log-likelhood values}
}
\references{
C. Bouveyron & J. Jacques, Model-based Clustering of Time Series in Group-specific Functional Subspaces, Advances in Data Analysis and Classification, vol. 5 (4), pp. 281-300, 2011.
}
\author{
Charles Bouveyron \& Julien Jacques
}
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
\keyword{clustering}
\keyword{functional data}% __ONLY ONE__ keyword per line
