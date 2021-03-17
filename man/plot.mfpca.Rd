\name{plot.mfpca}
\alias{plot.mfpca}
\title{Graphical representation for MFPCA}
\usage{
  \method{plot}{mfpca}(x,nharm,threshold,...)
}
\description{
It provides graphical representations for MFPCA: smoothed data plots, scores plots, variation of the mean curve      and eigenfunction plots.
}
\arguments{
  \item{x}{An object produced by mfpca function.}
  \item{nharm}{The number of harmonics for which you want graphical representations, default value is 3.}
  \item{threshold}{The threshold of proportion of variance that stop plotting, default value is 0.05.}
  \item{...}{Some additional parameters.}
}
\value{
\item{Data plot}{Plot of all smooth curves for each functional variable.}
\item{Scores plots}{Plot of curves coordinates on the number of eigenfunctions/harmonics selected and depending on the threshold choosen.}
\item{Mean curve plots}{Plot of variation of the mean curve. Variations are estimated based on the mean values with addition (red) and substraction (blue) of each eigenfunction.}
\item{Eigenfunction plots}{ visualization of each eigenfunction/harmonic selected and depending on the threshold choosen.}
}
\examples{
  ####Univariate case: "Canadian temperature" data (Ramsey & Silverman)
  daybasis65 <- create.fourier.basis(c(0, 365), nbasis=65, period=365)
  daytempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"], daybasis65,
                            fdnames=list("Day", "Station", "Deg C"))$fd

  res.pca<-mfpca(daytempfd)
  plot(res.pca,nharm=4)


  ####Multivariate case: "Canadian temperature" data (Ramsey & Silverman)
  daybasis65 <- create.fourier.basis(c(0, 365), nbasis=65, period=365)
  daytempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"], daybasis65,
                            fdnames=list("Day", "Station", "Deg C"))$fd
  dayprecfd<-smooth.basis(day.5, CanadianWeather$dailyAv[,,"Precipitation.mm"], daybasis65,
                          fdnames=list("Day", "Station", "Mm"))$fd

  res.pca<-mfpca(list(daytempfd,dayprecfd))
  plot(res.pca,nharm=4)

}
