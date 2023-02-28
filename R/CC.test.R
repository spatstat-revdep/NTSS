#' Random shift test of independence in a bivariate random field
#'
#' @description Nonparametric test of independence between a pair of random fields based on random shifts.
#' Either the torus correction or the variance correction can be used, see Mrkvička et al. (2021).
#' The variance correction is recommended for random fields with nontrivial autocorrelations.
#'
#' @details The test statistic can be either the sample covariance or the sample Kendall's or
#' Pearson's correlation coefficient. The choice of the test statistic is given
#' by the argument \code{type}.
#'
#' The torus correction can be applied for rectangular windows. On the other hand,
#' the variance correction is applicable both for rectangular and for irregular windows.
#' The choice of the correction is given
#' by the argument \code{correction}. Based on the simulation studies in Mrkvička et al. (2021),
#' the variance correction is recommended for random fields with nontrivial autocorrelations.
#'
#' The two realizations of the random fields (defined on the same domain) to be tested
#' should be supplied in the \code{covariateA, covariateB} arguments
#' as objects of the class \code{im} as used in the \code{spatstat} package.
#' The pattern of sampling points in which the test statistic is evaluated
#' is given in the argument \code{test.points} as an object of the class \code{ppp}.
#'
#' The shift vectors are generated from the
#' uniform distribution on the disk with radius given by the argument \code{radius}
#' and centered in the origin. The argument \code{verbose} determines if
#' auxiliary information and plots should be provided.
#'
#' @import spatstat
#' @import spatstat.random
#' @import spatstat.geom
#' @import stats
#' @import graphics
#'
#' @param covariateA first random field (object of class \code{im})
#' @param covariateB second random field (object of class \code{im})
#' @param test.points point pattern providing the set of sampling points (object of class \code{ppp})
#' @param N.shifts integer, how many random shifts should be performed in the random shift test
#' @param radius positive real number determining the radius of the disk on which the shift vectors are uniformly distributed
#' @param correction which correction should be applied in the random shift test (possible choices are "torus" and "variance")
#' @param type which test statistic should be used (possible choices are "Kendall", "Pearson" and "covariance")
#' @param verbose logical value indicating whether auxiliary information should be printed and auxiliary figures plotted during the computation
#'
#' @return The p-value of the random shift test of independence between a pair of random fields.
#'
#' @references T. Mrkvička, J. Dvořák, J.A. González, J. Mateu (2021): Revisiting the random shift approach for testing in spatial statistics. Spatial Statistics 42, 100430.
#'
#' @examples
#'
#' library(spatstat)
#'
#' set.seed(123)
#'
#' elevation <- bei.extra$elev
#' slope <- bei.extra$grad
#' plot(elevation)
#' plot(slope)
#'
#' test.points <- runifpoint(100, win=bei$window)
#'
#' out1 <- CC.test(covariateA=elevation, covariateB=slope, test.points=test.points, N.shifts=999,
#'                 radius=250, type="Kendall", correction="torus", verbose=TRUE)
#' out1
#'
#' out2 <- CC.test(covariateA=elevation, covariateB=slope, test.points=test.points, N.shifts=999,
#'                 radius=250, type="Kendall", correction="variance", verbose=TRUE)
#' out2
#'
#' @export
CC.test <- function(covariateA, covariateB, test.points, N.shifts=999, radius, correction, type="Kendall", verbose=FALSE){

  if (verbose){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  ### C-C test with torus correction

  if (correction=="torus"){

    values.simulated <- rep(NA, times=N.shifts+1)
    values.simulated[1] <- switch(type,
                                  covariance = cov(covariateA[test.points],covariateB[test.points]),
                                  Pearson    = cor(covariateA[test.points],covariateB[test.points],method="pearson"),
                                  Kendall    = cor(covariateA[test.points],covariateB[test.points],method="kendall"))

    test.stat <- values.simulated[1]
    names(test.stat) <- "T"

    correct <- "torus"
    names(correct) <- "correction"

    for (k in 1:N.shifts){
      test.points.shift <- rshift(test.points, edge="torus", radius=radius)
      values.simulated[k+1] <- switch(type,
                                      covariance = cov(covariateA[test.points.shift],covariateB[test.points]),
                                      Pearson    = cor(covariateA[test.points.shift],covariateB[test.points],method="pearson"),
                                      Kendall    = cor(covariateA[test.points.shift],covariateB[test.points],method="kendall"))
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


  ### C-C test with variance correction

  if (correction == "variance"){

    values.simulated <- rep(NA, times=N.shifts+1)
    values.simulated[1] <- switch(type,
                                  covariance = cov(covariateA[test.points],covariateB[test.points]),
                                  Pearson    = cor(covariateA[test.points],covariateB[test.points],method="pearson"),
                                  Kendall    = cor(covariateA[test.points],covariateB[test.points],method="kendall"))
    n.simulated <- rep(NA, times=N.shifts+1)
    n.simulated[1] <- test.points$n

    for (k in 1:N.shifts){
      jump <- runifdisc(1, radius = radius)
      test.points.shifted <- shift(test.points, c(jump$x,jump$y))
      W.reduced <- intersect.owin(test.points$window, test.points.shifted$window)
      test.points.reduced <- test.points[W.reduced]
      test.points.reduced.backshifted <- shift(test.points.reduced, c(-jump$x,-jump$y))
      covA <- covariateA[test.points.reduced.backshifted]
      covB <- covariateB[test.points.reduced]
      n.simulated[k+1] <- test.points.reduced$n
      values.simulated[k+1] <- switch(type,
                                      covariance = cov(covA,covB),
                                      Pearson    = cor(covA,covB,method="pearson"),
                                      Kendall    = cor(covA,covB,method="kendall"))
    }

    values.std <- (values.simulated - mean(values.simulated))*sqrt(n.simulated)
    test.rank <- rank(values.std)[1]
    pval <- 2*min(test.rank, N.shifts+1-test.rank)/(N.shifts+1)

    test.stat <- values.std[1]
    names(test.stat) <- "T_std"

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
  testname <- "Random shift test of independence between a pair of random fields (covariates)"
  alternative <- "two-sided"
  stat.type <- switch(type,
                      covariance = "sample covariance",
                      Pearson    = "Pearson's correlation coefficient",
                      Kendall    = "Kendall's correlation coefficient")

  result <- structure(list(statistic = test.stat, parameter = list(N.shifts=res.N.shifts,correction=correct,statistic=stat.type),
                           p.value = pval, method = testname, data.name = paste(substitute(covariateA), "and", substitute(covariateB)),
                           alternative = alternative),
                      class = "htest")

  return(result)
}
