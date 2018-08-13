\name{slopeHeuristic}
\alias{slopeHeuristic}
\title{Calculate slope heuristic}
\usage{
  slopeHeuristic(mod)
}
\description{
  It provides slope heuristic for number of groups selection or model selection (see Schmutz et al, 2017).
}
\arguments{
  \item{mod}{A funHDDC object with different number of groups tested or different models.}
  }
    \value{
    The index of the best model, among all estimated models.

    }
    \details{
    This function works only if all tested model stored in 'mod' have converged to a solution. If at least one         model is written ' -Inf pop<min.individuals' in mod output, an error message will be printed by slopeHeuristic     function.
    }
\examples{
    ##Multivariate example
    data("triangle")
    basis <- create.bspline.basis(c(1,21), nbasis=25)
    var1<-smooth.basis(argvals=seq(from=1,to=21,length.out = 101),y=t(triangle[,1:101]),
                        fdParobj=basis)$fd
    var2<-smooth.basis(argvals=seq(from=1,to=21,length.out = 101),y=t(triangle[,102:202]),
                        fdParobj=basis)$fd
    res.multi<-funHDDC(list(var1,var2),K=2:10,model="AkjBQkDk")
    slopeHeuristic(res.multi)

  \donttest{
    ##Univariate example
    data("trigo")
    basis<- create.bspline.basis(c(0,1), nbasis=25)
    var1<-smooth.basis(argvals=seq(0,1,length.out = 100),y=t(trigo[,1:100]),fdParobj=basis)$fd
    res.uni<-funHDDC(var1,K=2:10,model="AkBkQkDk")
    slopeHeuristic(res.uni)
    
    ##Clustering the "Canadian temperature" data (Ramsey & Silverman): univariate case
    daybasis65 <- create.fourier.basis(c(0, 365), nbasis=65, period=365)
    daytempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"], daybasis65,
    fdnames=list("Day", "Station", "Deg C"))$fd

    res.uni<-funHDDC(daytempfd,K=2:10,model="AkjBkQkDk",threshold=0.2)

    slopeHeuristic(res.uni)

    ##Clustering the "Canadian temperature" data (Ramsey & Silverman): multivariate case
    daybasis65 <- create.fourier.basis(c(0, 365), nbasis=65, period=365)
    daytempfd <- smooth.basis(day.5, CanadianWeather$dailyAv[,,"Temperature.C"], daybasis65,
    fdnames=list("Day", "Station", "Deg C"))$fd
    dayprecfd<-smooth.basis(day.5, CanadianWeather$dailyAv[,,"Precipitation.mm"], daybasis65,
    fdnames=list("Day", "Station", "Mm"))$fd

    res.multi<-funHDDC(list(daytempfd,dayprecfd),K=2:8,model="AkjBkQkDk",
    init="random",threshold=0.2)

    slopeHeuristic(res.multi)
    }
}
