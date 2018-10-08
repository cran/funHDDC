\name{mfpca}
\alias{mfpca}
\title{Functional principal component analysis for univariate or multivariate functional data}
\usage{
  mfpca(fdobj,center)
}
\description{
  It provides functional principal component analysis for univariate or multivariate functional data.
}
\arguments{
  \item{fdobj}{For univariate FPCA: a functional data object produced by fd() function of fda package, for
   multivariate FPCA: a list of functional data objects.}
  \item{center}{If TRUE (default), it centers each lines of data coefficients by the mean before calculating the       FPCA. }

}
\value{
  \item{eigval}{A list of eigen values.}
  \item{harmonics}{A functional data object for the harmonics or eigenfunctions.}
  \item{scores}{A matrix of scores on the harmonics.}
  \item{varprop}{A vector giving the proportion of variance explained by each harmonic.}
  \item{meanfd}{A functional data object giving the mean function after centering (default) or the mean function of     raw data.}

}
\examples{

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
