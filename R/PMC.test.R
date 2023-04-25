#' Random shift test of independence between marks of a point process and a covariate
#'
#' @description Nonparametric test of independence between marks of a point process and a random field
#' (covariate) based on random shifts, see Dvořák et al. (2022).
#' Either the torus correction or the variance correction can be used, see Mrkvička et al. (2021).
#'
#' @details The test statistic can be either the sample covariance or the sample Kendall's or
#' Pearson's correlation coefficient, computed from the vector of marks and
#' the vector of covariate values observed at the points of the process, see the paper
#' Dvořák et al. (2022). These test statistics make the test sensitive to violation
#' of independence between the marks and the covariate. On the other hand, the
#' test can also be viewed as a test of independence between a marked point process
#' and a covariate.
#'
#' The torus correction can be applied for rectangular windows. On the other hand,
#' the variance correction is applicable both for rectangular and for irregular windows.
#' The choice of the correction is given by the argument \code{correction}.
#' Based on the simulation studies in Dvořák et al. (2022),
#' the variance correction is recommended since it does not exhibit the liberality of the torus correction.
#'
#' The choice of the test statistic is given
#' by the argument \code{type}. The observed marked point pattern should be supplied using the
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
#' @param X marked point pattern dataset (object of class \code{ppp})
#' @param covariate random field (object of class \code{im})
#' @param N.shifts integer, how many random shifts should be performed in the random shift test
#' @param radius positive real number determining the radius of the disk on which the shift vectors are uniformly distributed
#' @param correction which correction should be applied in the random shift test (possible choices are "torus" and "variance")
#' @param type ... which test statistic should be used (possible choices are "Kendall", "Pearson" and "covariance")
#' @param verbose logical value indicating whether auxiliary information should be printed and auxiliary figures plotted during the computation
#'
#' @return The p-value of the random shift test of independence between marks of a point process and a covariate.
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
#' slope <- bei.extra$grad
#' plot(elevation)
#' plot(slope)
#'
#' # use only part of the point pattern for this example
#' Xun <- rthin(bei,0.05)
#'
#' # use elevation values as marks
#' X <- Xun %mark% elevation[Xun]
#'
#' # use terrain gradient as covariate
#' covariate <- slope
#'
#' # tests run with only 99 shifts to speed up the computation
#' out1 <- PMC.test(X=X, covariate=covariate, N.shifts = 99, radius=250, correction="torus",
#'                  type="Kendall", verbose=TRUE)
#' out1
#'
#' out2 <- PMC.test(X=X, covariate=covariate, N.shifts = 99, radius=250, correction="variance",
#'                  type="Kendall", verbose=TRUE)
#' out2
#'
#' @export
#'
PMC.test <- function(X, covariate, N.shifts=999, radius, correction, type="Kendall", verbose=FALSE){

  if (verbose){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  ### PM-C test with torus correction

  if (correction=="torus"){

    values.simulated <- rep(NA, times=N.shifts+1)
    values.simulated[1] <- switch(type,
                                  covariance = cov(covariate[X],X$marks),
                                  Pearson    = cor(covariate[X],X$marks,method="pearson"),
                                  Kendall    = cor(covariate[X],X$marks,method="kendall"))

    test.stat <- values.simulated[1]
    names(test.stat) <- "T"

    correct <- "torus"
    names(correct) <- "correction"

    for (k in 1:N.shifts){
      X.shift <- rshift(X, edge="torus", radius=radius)
      values.simulated[k+1] <- switch(type,
                                      covariance = cov(covariate[X.shift],X.shift$marks),
                                      Pearson    = cor(covariate[X.shift],X.shift$marks,method="pearson"),
                                      Kendall    = cor(covariate[X.shift],X.shift$marks,method="kendall"))
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


  ### PM-C test with variance correction

  if (correction == "variance"){

    values.simulated <- rep(NA, times=N.shifts+1)
    values.simulated[1] <- switch(type,
                                  covariance = cov(covariate[X],X$marks),
                                  Pearson    = cor(covariate[X],X$marks,method="pearson"),
                                  Kendall    = cor(covariate[X],X$marks,method="kendall"))
    n.simulated <- rep(NA, times=N.shifts+1)
    n.simulated[1] <- X$n

    for (k in 1:N.shifts){
      jump <- runifdisc(1, radius = radius)
      X.shifted <- shift(X, c(jump$x,jump$y))
      W.reduced <- intersect.owin(X$window, X.shifted$window)
      X.reduced <- X.shifted[W.reduced]
      n.simulated[k+1] <- X.reduced$n
      values.simulated[k+1] <- switch(type,
                                      covariance = cov(covariate[X.reduced],X.reduced$marks),
                                      Pearson    = cor(covariate[X.reduced],X.reduced$marks,method="pearson"),
                                      Kendall    = cor(covariate[X.reduced],X.reduced$marks,method="kendall"))
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
  testname <- "Random shift test of independence between marks and a covariate"
  alternative <- "two-sided"
  stat.type <- switch(type,
                      covariance = "sample covariance",
                      Pearson    = "Pearson's correlation coefficient",
                      Kendall    = "Kendall's correlation coefficient")

  result <- structure(list(statistic = test.stat, parameter = list(N.shifts=res.N.shifts,correction=correct,statistic=stat.type),
                           p.value = pval, method = testname, data.name = paste(substitute(X), "and", substitute(covariate)),
                           alternative = alternative),
                      class = "htest")

  return(result)
}
