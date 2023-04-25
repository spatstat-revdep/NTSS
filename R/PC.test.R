#' Random shift test of independence between a point process and a covariate
#'
#' @description Nonparametric test of independence between a point process and a random field (covariate) based on random shifts, see Dvořák et al. (2022).
#' Either the torus correction or the variance correction can be used, see Mrkvička et al. (2021).
#'
#' @details The test statistic is the mean covariate value observed at the points of the process, see the paper
#' Dvořák et al. (2022).
#'
#' The torus correction can be applied for rectangular windows. On the other hand,
#' the variance correction is applicable both for rectangular and for irregular windows.
#' The choice of the correction is given by the argument \code{correction}.
#' Based on the simulation studies in Dvořák et al. (2022),
#' the variance correction is recommended since it does not exhibit the liberality of the torus correction.
#'
#' The observed point pattern should be supplied using the
#' argument \code{X}, the realization of the covariate should be supplied using
#' the argument \code{covariate}.
#'
#' The shift vectors are generated from the
#' uniform distribution on the disk with radius given by the argument \code{radius}
#' and centered in the origin. The argument \code{verbose} determines if
#' auxiliary information and plots should be provided.
#'
#' @import spatstat
#' @import stats
#'
#' @param X point pattern dataset (object of class \code{ppp})
#' @param covariate random field (object of class \code{im})
#' @param N.shifts integer, how many random shifts should be performed in the random shift test
#' @param radius positive real number determining the radius of the disk on which the shift vectors are uniformly distributed
#' @param correction which correction should be applied in the random shift test (possible choices are "torus" and "variance")
#' @param verbose logical value indicating whether auxiliary information should be printed and auxiliary figures plotted during the computation
#'
#' @return The p-value of the random shift test of independence between a point process and a covariate.
#'
#' @references J. Dvořák, T. Mrkvička, J. Mateu, J.A. González (2022): Nonparametric testing of the dependence structure among points-marks-covariates in spatial point patterns. International Statistical Review 90(3), 592-621.
#' @references T. Mrkvička, J. Dvořák, J.A. González, J. Mateu (2021): Revisiting the random shift approach for testing in spatial statistics. Spatial Statistics 42, 100430.
#'
#' @examples
#'
#' library(spatstat)
#'
#' set.seed(123)
#'
#' elevation <- bei.extra$elev
#' plot(elevation)
#'
#' # tests run with only 99 shifts to speed up the computation
#' out1 <- PC.test(X=bei, covariate=elevation, N.shifts = 99, radius=250,
#'                 correction="torus", verbose=TRUE)
#' out1
#'
#' out2 <- PC.test(X=bei, covariate=elevation, N.shifts = 99, radius=250,
#'                 correction="variance", verbose=TRUE)
#' out2
#'
#' @export
#'
PC.test <- function(X, covariate, N.shifts=999, radius, correction, verbose=FALSE){

  if (verbose){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  ### P-C test with torus correction

  if (correction=="torus"){

    values.simulated <- rep(NA, times=N.shifts+1)
    values.simulated[1] <-  mean(covariate[X])

    test.stat <- values.simulated[1]
    names(test.stat) <- "mean"

    correct <- "torus"
    names(correct) <- "correction"

    for (k in 1:N.shifts){
      X.shift <- rshift(X, edge="torus", radius=radius)
      values.simulated[k+1] <- mean(covariate[X.shift])
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


  ### P-C test with variance correction

  if (correction == "variance"){

    values.simulated <- rep(NA, times=N.shifts+1)
    values.simulated[1] <- mean(covariate[X])
    n.simulated <- rep(NA, times=N.shifts+1)
    n.simulated[1] <- X$n

    # Random shifts
    for (k in 1:N.shifts){
      jump <- runifdisc(1, radius = radius)
      X.shifted <- shift(X, c(jump$x,jump$y))
      W.reduced <- intersect.owin(X$window, X.shifted$window)
      X.reduced <- X.shifted[W.reduced]
      values.simulated[k+1] <- mean(covariate[X.reduced])
      n.simulated[k+1] <- X.reduced$n
    }

    values.std <- (values.simulated - mean(values.simulated))*sqrt(n.simulated)
    test.rank <- rank(values.std)[1]
    pval <- 2*min(test.rank, N.shifts+1-test.rank)/(N.shifts+1)

    test.stat <- values.std[1]
    names(test.stat) <- "mean_std"

    correct <- "variance"
    names(correct) <- "correction"

    if (verbose){
      cat("Observed value of the test statistic (before variance standardization): ")
      cat(values.simulated[1])
      cat("\n")
      cat("Observed value of the test statistic (after variance standardization): ")
      cat(values.std[1])
      cat("\n")
      if (verbose){cat("Random shift with variance correction, p-value: ")}
      cat(pval)
      cat("\n")

      aux <- hist(values.std, plot=FALSE)
      par(mar=c(2,2,2,2))
      plot(aux, main="Test statistic values after random shifts")
      abline(v=values.std[1], col="red")
    }
  }

  res.N.shifts <- N.shifts
  names(res.N.shifts) <- "N.shifts"
  testname <- "Random shift test of independence between a point process and a covariate"
  alternative <- "two-sided"

  result <- structure(list(statistic = test.stat, parameter = list(N.shifts=res.N.shifts,correction=correct),
                           p.value = pval, method = testname, data.name = paste(substitute(X), "and", substitute(covariate)),
                           alternative = alternative),
                      class = "htest")

  return(result)
}
