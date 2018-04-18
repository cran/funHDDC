\name{mfpca}
\alias{mfpca}
\title{Functional Principal Component Analysis for Univariate or Multivariate Functional Data}
\usage{
  mfpca(fdobj,center)
}
\description{
  Functional principal component analysis allows to build a low-dimensional representation for univariate or multivariate functional data.
}
\arguments{
  \item{fdobj}{For univariate FPCA: a functional data object produced by fd() function of fda package, for
   multivariate FPCA: a list of functional data objects.}
  \item{center}{if TRUE (default), it centers each lines of data coefficients by the mean before calculating the       FPCA. }

}
\value{
  \item{eigval}{a list of eigen values.}
  \item{harmonics}{a functional data object for the harmonics or eigenfunctions.}
  \item{scores}{a matrix of scores on the harmonics.}
  \item{varprop}{a vector giving the proportion of variance explained by each harmonic.}
  \item{meanfd}{a functional data object giving the mean function after centering (default) or the mean function of     raw data.}

}
\references{
- C. Bouveyron and J. Jacques, Model-based Clustering of Time Series in Group-specific Functional Subspaces, Advances in Data Analysis and Classification, vol. 5 (4), pp. 281-300, 2011.
- A. Schmutz, C. Bouveyron, J. Jacques, L. Cheze and  P. Martin, Clustering multivariate functional data in group-specic functional subspaces, Preprint HAL 01652467, Université Côte d'Azur, 2017.
}
\examples{
  library(fda)

  ####Univariate case: "Canadian temperature" data (Ramsey & Silverman)
  daybasis65 <- create.fourier.basis(c(0, 365), nbasis=65, period=365)
  daytempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"], daybasis65,
                                  fdnames=list("Day", "Station", "Deg C"))$fd

  res.pca<-mfpca(daytempfd)
  plot.mfpca(res.pca)


  ####Multivariate case: "Canadian temperature" data (Ramsey & Silverman)
  daybasis65 <- create.fourier.basis(c(0, 365), nbasis=65, period=365)
  daytempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"], daybasis65,
                                  fdnames=list("Day", "Station", "Deg C"))$fd
  dayprecfd<-smooth.basis(day.5, CanadianWeather$dailyAv[,,"Precipitation.mm"], daybasis65,
                                fdnames=list("Day", "Station", "Mm"))$fd

  res.pca<-mfpca(list(daytempfd,dayprecfd))
  plot.mfpca(res.pca)

}
