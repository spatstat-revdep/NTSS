#' Test of independence between marks and points of a stationary marked point process
#'
#' @description Test of independence between marks and points of a stationary marked point process
#' based on the comparison of the mark K-function and the K-function of the underlying point process,
#' see Guan (2006).
#'
#' @details The test was proposed in Guan (2006) and used in Dvořák et al. (2022) as one of the tools 
#' to determine the dependence structure between points, marks and a covariate in a marked point process setting.
#' The marked point process is assumed to be stationary here.
#'
#' The observed marked point pattern should be supplied using the argument \code{X}.
#' The upper limit of the domain for the (marked) K-function used for constructing the test statistic
#' is given by the argument \code{R}. Furthermore, the variance used for constructing 
#' the test statistic is determined using a fitted covariance model for the mark field.
#' In the given implementation, the exponential model for the covariance is assumed.
#' The model parameters are estimated by fitting the corresponding parametric model
#' to the empirical variogram by the least-squares approach, with the upper bound on
#' the variogram argument given by \code{R.variogram}.
#'
#' @import spatstat
#' @import geoR
#'
#' @param X marked point pattern dataset (object of class \code{ppp})
#' @param R positive real number determining the domain for the (marked) K-function used for computing the test statistic
#' @param R.variogram positive real number determining the domain for variogram estimation
#' 
#' @return The p-value of the test of independence between marks and points of a stationary marked point process.
#'
#' @references Y. Guan (2006): Tests for independence between marks and points of a marked point Process. Biometrics 62, 126-134.
#' @references J. Dvořák, T. Mrkvička, J. Mateu, J.A. González (2022): Nonparametric testing of the dependence structure among points-marks-covariates in spatial point patterns. International Statistical Review 90(3), 592-621.
#'
#' @examples
#'
#' library(spatstat)
#' library(geoR)
#'
#' set.seed(123)
#'
#' X <- rpoispp(100)
#' X <- X %mark% runif(n=X$n)
#' 
#' out <- PM.test(X, R=0.1, R.variog=1)
#' out
#'
#' @export
#'
PM.test <- function(X, R, R.variogram){
  
  # Estimate the constant intensity
  lambda <- intensity(X)
  
  # Shift the mark values simultanously to avoid negative mark values, 
  # this does not affect the results (the respective terms cancel out)
  mX <- marks(X) - min(marks(X)) + 1
  
  # matrix of pairwise distances between points
  pd <- pairdist(X)
  
  
  ### Compute Int X(r) dr ###
  ###########################
  
  # weights for the inhomogeneous K-function (used since it allows using weights)
  wt <- pmax((R-pd),0)/lambda^2
  
  K <- Kinhom(X, r = c(0,R), reciplambda2 = wt, correction = "translate", renormalise = FALSE)
  intX <- K$trans[2]
  
  
  ### Compute Int Y(r) dr ###
  ###########################
  
  # weights for the inhomogeneous K-function (used since it allows using weights)
  wt <- pmax((R-pd),0) * matrix(mX,nrow=X$n,ncol=X$n) / lambda^2
  
  K <- Kinhom(X, r = c(0,R), reciplambda2 = wt, correction = "translate", renormalise = FALSE)
  intY <- K$trans[2]
  
  
  ### Compute estimate of mu ###
  ##############################
  
  mu.hat <- intY/intX
  mu.hat.std <- mu.hat - mean(mX)
  
  
  ### Compute the weights for variance ###
  ########################################
  
  aux <- (R-pd)*(pd <= R)*edge.Trans(X)/(area(X$window)*lambda^2*intX)
  diag(aux) <- 0
  wt <- rowSums(aux) - 1/X$n
  
  
  ### Estimate covariance structure of marks ###
  ##############################################
  
  # Estimate the variogram, exponential model
  geo.cov <- data.frame(x=X$x, y=X$y, z=mX)
  geo.cov <- as.geodata(geo.cov)
  v <- variog(geo.cov, uvec=seq(from=0, to=R.variogram, length.out=21), messages=FALSE)
  obj <- function(par){# par[1] = sigma^2, par[2] = s
    sum( (v$v - par[1]*(1-exp(-v$u/par[2])) )^2)
  }
  
  # estimated parameters
  est.par <- optim(par=c(1,1), fn=obj)$par
  
  # covariance function with estimated parameters
  Cfun <- function(r){
    est.par[1]*exp(-r/est.par[2])
  }
  
  # matrix of covariance values
  Cmat <- Cfun(pd)
  
  
  ### Estimate the variance ###
  #############################
  
  var.hat <- as.numeric(wt %*% Cmat %*% wt)
  
  
  ### Compute the test statistic and perform the test ###
  #######################################################
  
  test.stat <- mu.hat.std^2 / var.hat
  pval <- pchisq(test.stat, df=1, lower.tail=FALSE)
  
  
  names(test.stat) <- "test statistic"
  
  testname <- "Test of independence between marks and points of a stationary marked point process"
  alternative <- "one-sided"
  
  result <- structure(list(statistic = test.stat, parameter = list(R=R, R.variogram=R.variogram),
                           p.value = pval, method = testname, data.name = substitute(X),
                           alternative = alternative),
                      class = "htest")
  
  return(result)
}
