#' Random shift test of independence in a bivariate point process
#'
#' @description Nonparametric test of independence between a pair of point processes based on random shifts.
#' Either the torus correction or the variance correction can be used (note that the
#' variance correction is not yet implemented in this package but the corresponding source
#' codes can be obtained from the authors upon request).
#'
#' @details The test statistic is the expectation of the cross nearest-neighbor distance, see the paper
#' Mrkvička et al. (2021). The test statistic based on the cross K-function is
#' not provided here due to its poor performance in the simulation experiments
#' described in the paper.
#'
#' The two observed point patterns to be tested should be supplied using the
#' arguments \code{X} and \code{Y}. The shift vectors are generated from the
#' uniform distribution on the disk with radius given by the argument \code{radius}
#' and centered in the origin. The argument \code{verbose} determines if
#' auxiliary information and plots should be provided.
#'
#' @import spatstat
#' @import spatstat.explore
#' @import stats
#'
#' @param X first point pattern dataset (object of class \code{ppp})
#' @param Y second point pattern dataset (object of class \code{ppp})
#' @param N.shifts integer, how many random shifts should be performed in the random shift test
#' @param radius positive real number determining the radius of the disk on which the shift vectors are uniformly distributed
#' @param correction which correction should be applied in the random shift test (possible choices are "torus" and "variance")
#' @param verbose logical value indicating whether auxiliary information should be printed and auxiliary figures plotted during the computation
#'
#' @return The p-value of the random shift test of independence between a pair of point processes.
#'
#' @references T. Mrkvička, J. Dvořák, J.A. González, J. Mateu (2021): Revisiting the random shift approach for testing in spatial statistics. Spatial Statistics 42, 100430.
#'
#' @examples
#'
#' library(spatstat)
#'
#' set.seed(123)
#'
#' X <- rpoispp(150)
#' Y <- rpoispp(150)
#'
#' out <- PP.test(X=X, Y=Y, N.shifts=499, radius=0.5, correction="torus", verbose=TRUE)
#' out
#'
#' @export
#'
PP.test <- function(X, Y, N.shifts=999, radius, correction, verbose=FALSE){

  if (verbose){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  # Z ... bivariate point pattern to be tested, with marks "X" and "Y"
  Z <- superimpose(X = X, Y = Y)


  ### P-P test with torus correction

  if (correction=="torus"){

    values.simulated <- rep(NA, times=N.shifts+1)
    fun <- Gcross(Z, correction	=c("km"))
    values.simulated[1] <- stieltjes(function(x){x}, fun)$km

    test.stat <- values.simulated[1]
    names(test.stat) <- "EG_12"

    correct <- "torus"
    names(correct) <- "correction"

    for (k in 1:N.shifts){

      Z.shift <- rshift.ppp(Z, which='X', radius=radius, edge="torus")
      fun <- Gcross(Z.shift, correction	=c("km"))
      values.simulated[k+1] <- stieltjes(function(x){x}, fun)$km
    }

    test.rank <- rank(values.simulated)[1]
    pval <- 2*min(test.rank, N.shifts+1-test.rank)/(N.shifts+1)

    if (verbose){
      cat("Observed value of the test statistic: ")
      cat(values.simulated[1])
      cat("\n")
      cat("Random shift with torus correction, p-value: ")
      cat(pval)
      cat("\n")
      aux <- hist(values.simulated, plot=FALSE)
      par(mar=c(2,2,2,2))
      plot(aux, main="Test statistic values after random shifts")
      abline(v=values.simulated[1], col="red")
    }
  }


  ### P-P test with variance correction

  if (correction == "variance"){
    print("Variance correction for this test is not yet implemented but the corresponding source
           codes can be obtained from the authors upon request.")
    return(NULL)
  }

  res.N.shifts <- N.shifts
  names(res.N.shifts) <- "N.shifts"
  testname <- "Random shift test of independence between a pair of point processes"
  alternative <- "two-sided"

  result <- structure(list(statistic = test.stat, parameter = list(N.shifts=res.N.shifts,correction=correct),
                           p.value = pval, method = testname, data.name = paste(substitute(X), "and", substitute(Y)),
                           alternative = alternative),
                      class = "htest")

  return(result)
}
