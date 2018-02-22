#' Calculate fraction seawater variability
#'
#' @param d18O.coral
#' @param SSS
#' @param SST
#' @param SST.coef
#' @param SSS.coef
#' @references Russon et al., 2013
#' @return Fsw The fraction of d18O coral variance explained by salinity variability
#' @export
#'
#' @examples
Fsw = function(d18O.coral, SSS, SST, SST.coef = -0.23, SSS.coef = 0.60){
  d18O.sw = SSS*SSS.coef
  return(  (var(d18O.sw) + 2*SST.coef*cov(SST,d18O.sw)) /
             var(d18O.coral) )
}
