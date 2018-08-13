slopeHeuristic <- function(mod){
  # Setting up the data for estimation => selection of valid models
  main_data = mod$allCriteria
  who_notNA_norInfinite = !is.na(main_data$complexity) & is.finite(main_data$LL)

  n_valid = sum(who_notNA_norInfinite)
  if(n_valid == 0){
    stop("There is not any valid model to be selected.")
  } else if(n_valid <= 2){
    stop("At least 3 valid models are necessary to perform the slope heuristic. Otherwise, use another criterion.")
  }

  # Robust estimation
  main_data = main_data[who_notNA_norInfinite, ]
  fit = MASS::rlm(LL ~ complexity, data=main_data, method='MM')

  # To avoid problems (no negative slope allowed => weird anyway)
  fit_coef = fit$coefficients
  if(fit_coef[2]<0){
    fit_coef[2]=0
  }

  # the new penalized likelihood
  llpen = main_data$LL- 2* fit_coef[2]*main_data$complexity
  SH = 2* llpen
  plot(main_data$complexity,main_data$LL,type='p',xlab='Model dimension',ylab='Log-likelihood')
  abline(fit,col='red')
  plot(main_data$K,llpen,type='p',xlab='K',ylab='Penalized log-likelihood')
  points(main_data$K[which.max(llpen)],max(llpen),pch=19,col='red')
  return(main_data$K[which.max(llpen)])
}

