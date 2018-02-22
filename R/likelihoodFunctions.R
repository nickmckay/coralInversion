


#set up likelihood functions
#' Likelihood function for Fsw
#' @param simFsw
#' @param obsFsw
#' @param uncFsw
#'
#' @return
#' @export
#'
#' @examples
Fsw.likelihood = function(simFsw,obsFsw,uncFsw = 0.1){
  #likelihood function for Fsw
  #assume normal distribution of Fsw for now
  like = dnorm(simFsw,mean = obsFsw, sd = uncFsw)
  return(like)
}

#' likelihood function for matching instrumental SST
#'
#' @param simSST
#' @param obsSST
#' @param uncSST
#'
#' @return
#' @export
#'
#' @examples
sst.likelihood = function(simSST,obsSST,uncSST = 2){
  scaledDiffs = (simSST-mean(simSST)) - (obsSST - mean(obsSST))
  RMSE = sqrt(t(scaledDiffs)%*%scaledDiffs/length(simSST))
  #assume normal disribution for SSE for now
  like = dnorm(RMSE,sd = uncSST)
  return(like)
}

#covariance likelihood
#' Likelihood function for matching SSS/SST covariance
#'
#' @param simCov
#' @param obsCov
#' @param unc
#'
#' @return
#' @export
#'
#' @examples
cov.likelihood = function(simCov,obsCov,unc = 1){
  #assume normal disribution for SSE for now
  like = dnorm(simCov,mean = obsCov, sd = unc)
  return(like)
}
