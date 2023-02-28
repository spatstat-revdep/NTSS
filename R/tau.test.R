#' Random shift test of independence between a point process and a covariate
#'
#' @description Nonparametric test of independence between a point process and a random field
#' (covariate of interest), taking into account the possible effect of nuisance covariates,
#' see Dvořák and Mrkvička (2022).
#' The test is based on random shifts. Either the torus correction or the variance
#' correction can be used, see Mrkvička et al. (2021).
#' This test has lower power than the test based on the covariate-weighted residuals
#' (see the function \code{CWR.test}), but it is still included for the sake of completeness.
#' Also, the test statistic of \code{tau.test} can be used to quantify
#' the partial correlation between the point process and the covariate of interest, taking into
#' account the possible effect of nuisance covariates, see also the function \code{tau.est}.
#'
#' @details The test statistic is the Kendall's correlation coefficient between the covariate of interest
#' and the smoothed residual field, sampled at a given number of test points
#' scattered independently in the observation window, see the paper Dvořák and Mrkvička (2022).
#' If no nuisance covariates are given, the null model
#' assumes a constant intensity function of the point process. If one or more nuisance covariates are
#' provided, the null model assumes an intensity function depending on the nuisance
#' covariates (but not on the covariate of interest) and the residuals are constructed
#' using this intensity function.
#'
#' For constructing the smoothed residual field an adaptive bandwidth selection
#' method can be used, see Dvořák and Mrkvička (2022). A vector of candidate
#' bandwidth values can be provided using the argument \code{bws}.
#'
#' The residuals can be constructed in a nonparametric way (see Baddeley et al. (2012))
#' or in a parametric way (using the \code{ppm} function from the \code{spatstat} package,
#' see Baddeley et al. (2015)). This choice is given by the argument \code{nonparametric}.
#' The raw residuals are considered here.
#'
#' The torus correction can be applied for rectangular windows. On the other hand,
#' the variance correction is applicable both for rectangular and for irregular windows.
#' The choice of the correction is given by the argument \code{correction}.
#' Based on the simulation studies in Dvořák and Mrkvička (2022),
#' the variance correction is recommended since it does not exhibit the liberality of the torus correction.
#'
#' The observed point pattern should be supplied using the argument \code{X}.
#' The realization of the covariate of interest should be supplied using
#' the argument \code{covariate.interest}. The set of nuisance covariates should
#' be supplied as a list using the argument \code{covariates.nuisance}. This list
#' can be empty if no nuisance covariates are considered.
#'
#' The shift vectors are generated from the
#' uniform distribution on the disk with radius given by the argument \code{radius}
#' and centered in the origin. The argument \code{verbose} determines if
#' auxiliary information and plots should be provided.
#'
#' In case the observation window accompanying the point pattern is irregular,
#' it must be specified in the form of a binary mask due to the specific implementation of the test.
#' For details on binary masks see the help for the \code{spatstat} function \code{owin}.
#'
#' @import spatstat
#' @import spatstat.explore
#' @import spatstat.model
#' @import stats
#' @import ks
#'
#' @param X point pattern dataset (object of class \code{ppp})
#' @param covariate.interest random field (object of class \code{im})
#' @param covariates.nuisance list of covariates (objects of class \code{im}) determining the nuisance covariates
#' @param N.shifts integer, how many random shifts should be performed in the random shift test
#' @param radius positive real number determining the radius of the disk on which the shift vectors are uniformly distributed
#' @param correction which correction should be applied in the random shift test (possible choices are "torus" and "variance")
#' @param bws vector of positive real values from which the bandwidth is adaptively chosen if at least one nuisance covariate is present; if no nuisance covariates are present, the first value is used
#' @param n.test.points the number of independent test points used in computing the test statistic value
#' @param nonparametric logical value indicating whether nonparametric residuals should be used when computing the test statistic
#' @param verbose logical value indicating whether auxiliary information should be printed and auxiliary figures plotted during the computation
#' @param bw.factor.rhonhat multiplicative factor used when determining the bandwidth in the nonparametric estimation of the intensity function depending on the nuisance covariates (defaults to 1)
#'
#' @return The p-value of the random shift test of independence between a point process and a covariate, taking into account possible effects of nuisance covariates.
#'
#' @references J. Dvořák, T. Mrkvička (2022): Nonparametric testing of the covariate significance for spatial point patterns under the presence of nuisance covariates. https://arxiv.org/abs/2210.05424
#' @references T. Mrkvička, J. Dvořák, J.A. González, J. Mateu (2021): Revisiting the random shift approach for testing in spatial statistics. Spatial Statistics 42, 100430.
#' @references A. Baddeley, E. Rubak, R. Turner (2015) Spatial Point Patterns: Methodology and Applications with R. Chapman & Hall Interdisciplinary Statistics Series. CRC Press, Boca Raton, Florida.
#' @references A. Baddeley, Y.-M. Chang, Y. Song, R. Turner (2012) Nonparametric estimation of the dependence of a point process on spatial covariates. Statistics and Its Interface 5(2), 221?236.
#'
#' @examples
#'
#' library(spatstat)
#' library(ks)
#'
#' # the point pattern
#' X <- bei
#' plot(X)
#'
#' # two covariates are available
#' elevation <- bei.extra$elev
#' slope <- bei.extra$grad
#' plot(elevation)
#' plot(slope)
#'
#' # candidate values for adaptive bandwidth selection
#' bws <- seq(from=12.5, to=100, by=12.5)
#'
#' # test with no nuisance covariates, with only 99 shifts to speed up the computation
#' out1 <- tau.test(X, covariate.interest=elevation, covariates.nuisance=NULL,
#'                  bws=bws, N.shifts = 99, verbose=TRUE, correction="torus", radius=250)
#' out1
#'
#' # test with one nuisance covariate, with only 99 shifts to speed up the computation
#' out2 <- tau.test(X, covariate.interest=elevation, covariates.nuisance=list(slope=slope),
#'                  bws=bws, N.shifts = 99, verbose=TRUE, correction="torus", radius=250)
#' out2
#'
#' @export
#'
tau.test <- function(X, covariate.interest, covariates.nuisance, N.shifts=999, radius, correction,
                     bws, n.test.points=1000, nonparametric=TRUE, verbose=FALSE, bw.factor.rhonhat=1){

  if ((correction=="torus") & (!is.rectangle(X$window))){stop("Torus correction only applicable for rectangular windows!")}

  if (verbose){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  ### Generate independent grid of test points

  test.points <- runifpoint(n.test.points, win=as.mask(density.ppp(X)))


  ### Computing the smoothed residual field

  ncovs <- length(covariates.nuisance)
  nbws <- length(bws)

  if (ncovs == 0){ # no covariates, use homogenenous model
    intX <- intensity(X)

    # bw is chosen as the first value in the vector "bws"
    bw <- bws[1]
    srf <- density.ppp(X, sigma=bw)
    srf$v <- srf$v - intX

    if (verbose){
      cat("Using model with constant intensity function for construction of residuals.\n")
      cat("User-specified value of bandwidth used: ")
      cat(bw)
      cat("\n")
    }

  } else { # at least one nuisance covariate present

    if (nonparametric){ # nonparametric estimate of intensity function

      est.rho <- rhonhat(X, covariates.nuisance, prediction.grid.from=covariate.interest, bw.factor=bw.factor.rhonhat)

      aux.bw.cor <- matrix(NA, nrow=nbws, ncol=ncovs)
      for (k.bw in 1:nbws){
        aux.bw <- bws[k.bw]
        dens <- density.ppp(X, sigma=aux.bw, W=as.owin(covariate.interest))
        sm.est.rho <- Smooth(est.rho, normalise=TRUE, sigma=aux.bw)
        res <- suppressWarnings(eval.im(dens-sm.est.rho))
        for (cov.rep in 1: ncovs){
          aux.bw.cor[k.bw, cov.rep] <- cor(res[test.points],covariates.nuisance[[cov.rep]][test.points],method="kendall")
        }
      }
      cor.sums <- rowSums(aux.bw.cor*aux.bw.cor)
      k.which <- which.min(cor.sums)
      bw <- bws[k.which]
      dens <- density.ppp(X, sigma=bw, W=as.owin(covariate.interest))
      sm.est.rho <- Smooth(est.rho, normalise=TRUE, sigma=bw)
      srf <- eval.im(dens-sm.est.rho)

      if (verbose){
        cat("Using nonparametric estimate of intensity function for construction of residuals.\n")
        cat("Data-driven value of bandwidth used: ")
        cat(bw)
        cat("\n")
      }

    } else { # parametric (log-linear) estimate of intensity function
      aux.bw.cor <- matrix(NA, nrow=nbws, ncol=ncovs)
      fit <- ppm(X ~ ., covariates=covariates.nuisance)
      res.fit <- residuals(fit)
      for (k.bw in 1:nbws){
        aux.bw <- bws[k.bw]
        res <- Smooth(res.fit, sigma=aux.bw)
        for (cov.rep in 1: ncovs){
          aux.bw.cor[k.bw, cov.rep] <- cor(res[test.points],covariates.nuisance[[cov.rep]][test.points],method="kendall")
        }
      }
      cor.sums <- rowSums(aux.bw.cor*aux.bw.cor)
      k.which <- which.min(cor.sums)
      bw <- bws[k.which]
      srf <- Smooth(res.fit, sigma=bw)
      if (verbose){
        cat("Using parametric (log-linear) estimate of intensity function for construction of residuals.\n")
        cat("Data-driven value of bandwidth used: ")
        cat(bw)
        cat("\n")
      }
    }

    if (verbose){
      par(mar=c(4,4,2,1))
      if (ncovs==1){
        plot(x=bws, y=cor.sums, xlab="bw", ylab="square of Kendall's tau", main="Adaptive procedure for bandwidth selection",
             ylim=c(-max(cor.sums),max(cor.sums)))
      } else {
        plot(x=bws, y=cor.sums, xlab="bw", ylab="sum of squares of Kendall's tau", main="Adaptive procedure for bandwidth selection",
             ylim=c(-max(cor.sums),max(cor.sums)))
      }
      points(x=bws[k.which], y=cor.sums[k.which], pch=19, col="red")
      abline(h=0, col="red")
    }
  }


  ### tau test with torus correction

  if (correction=="torus"){

    if (verbose){
      # Estimated value of the Kendall's correlation coefficient for the covariate of interest
      res <- cor(srf[test.points],covariate.interest[test.points],method="kendall")

      cat("Observed value of the test statistic: ")
      cat(res)
      cat("\n")
      cat("Random shift with torus correction, p-value: ")
    }

    values.simulated <- rep(NA, times=N.shifts+1)
    values.simulated[1] <- cor(srf[test.points],covariate.interest[test.points],method="kendall")

    test.stat <- values.simulated[1]
    names(test.stat) <- "tau"

    correct <- "torus"
    names(correct) <- "correction"

    for (k in 1:N.shifts){
      test.points.shift <- rshift(test.points, edge="torus", radius=radius)
      values.simulated[k+1] <- cor(srf[test.points.shift],covariate.interest[test.points.shift],method="kendall")
    }
    test.rank <- rank(values.simulated)[1]
    pval <- 2*min(test.rank, N.shifts+1-test.rank)/(N.shifts+1)

    if (verbose){
      cat(pval)
      cat("\n")
      aux <- hist(values.simulated, plot=FALSE)
      par(mar=c(2,2,2,2))
      plot(aux, main="Test statistic values after random shifts")
      abline(v=values.simulated[1], col="red")
    }
  }


  ### tau test with variance correction

  if (correction == "variance"){

    if (verbose){
      # Estimated value of the Kendall's correlation coefficient for the covariate of interest
      res <- cor(srf[test.points],covariate.interest[test.points],method="kendall")

      cat("Observed value of the test statistic (before variance standardization): ")
      cat(res)
      cat("\n")
    }

    values.simulated <- rep(NA, times=N.shifts+1)
    values.simulated[1] <- cor(srf[test.points],covariate.interest[test.points],method="kendall")
    n.simulated <- rep(NA, times=N.shifts+1)
    n.simulated[1] <- X$n

    for (k in 1:N.shifts){
      jump <- runifdisc(1, radius=radius)
      test.points.shifted <- shift(test.points, c(jump$x,jump$y))
      W.reduced <- intersect.owin(test.points$window, test.points.shifted$window)
      test.points.reduced <- test.points[W.reduced]
      test.points.reduced.backshifted <- shift(test.points.reduced, c(-jump$x,-jump$y))
      vals1 <- srf[test.points.reduced.backshifted, drop=FALSE]
      vals2 <- covariate.interest[test.points.reduced, drop=FALSE]
      values.simulated[k+1] <- cor(vals1, vals2, method="kendall", use="complete.obs")
      n.simulated[k+1] <- sum((!is.na(vals1))*(!is.na(vals2)))
    }
    values.std <- (values.simulated - mean(values.simulated))*sqrt(n.simulated)
    test.rank <- rank(values.std)[1]
    pval <- 2*min(test.rank, N.shifts+1-test.rank)/(N.shifts+1)

    test.stat <- values.std[1]
    names(test.stat) <- "tau_std"

    correct <- "variance"
    names(correct) <- "correction"

    if (verbose){
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
  bw.value <- bw
  names(bw.value) <- "bandwidth"
  testname <- "Random shift test of independence between a point process and a covariate, with nuisance covariates"
  alternative <- "two-sided"

  # Note: is there a reasonable way how to include the name of the list of nuisance covariates, too?
  dn <- paste(substitute(X), "and", substitute(covariate.interest))

  result <- structure(list(statistic = test.stat, parameter = list(N.shifts=res.N.shifts,correction=correct,bandwidth=bw.value),
                           p.value = pval, method = testname, data.name = dn,
                           alternative = alternative),
                      class = "htest")

  return(result)
}
