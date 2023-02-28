#' (Partial) correlation coefficient between a point process and a covariate
#'
#' @description Correlation coefficient between a point process and a random field
#' (covariate of interest), taking into account the possible effect of nuisance covariates,
#' see Dvořák and Mrkvička (2022). The random shift test based on this test statistics
#' is given in the function \code{tau.test}.
#'
#' @details This function computes the Kendall's correlation coefficient between the covariate of interest
#' and the smoothed residual field, sampled at a given number of test points
#' scattered independently in the observation window, see the paper Dvořák and Mrkvička (2022).
#' If no nuisance covariates are given, a constant intensity function of the point process is assumed
#' when constructing the residuals. If one or more nuisance covariates are
#' provided, an intensity function depending on the nuisance
#' covariates (but not on the covariate of interest) is assumed and the residuals are constructed
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
#' The observed point pattern should be supplied using the argument \code{X}.
#' The realization of the covariate of interest should be supplied using
#' the argument \code{covariate.interest}. The set of nuisance covariates should
#' be supplied as a list using the argument \code{covariates.nuisance}. This list
#' can be empty if no nuisance covariates are considered.
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
#' @param bws vector of positive real values from which the bandwidth is adaptively chosen if at least one nuisance covariate is present; if no nuisance covariates are present, the first value is used
#' @param n.test.points the number of independent test points used in computing the test statistic value
#' @param nonparametric logical value indicating whether nonparametric residuals should be used when computing the test statistic
#' @param verbose logical value indicating whether auxiliary information should be printed and auxiliary figures plotted during the computation
#' @param bw.factor.rhonhat multiplicative factor used when determining the bandwidth in the nonparametric estimation of the intensity function depending on the nuisance covariates (defaults to 1)
#'
#' @return Value of the (partial) Kendall's correlation coefficient.
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
#' # no nuisance covariates
#' out1 <- tau.est(X, covariate.interest=elevation, covariates.nuisance=NULL,
#'                 bws=bws, verbose=TRUE)
#' out1
#'
#' # one nuisance covariate
#' out2 <- tau.est(X, covariate.interest=elevation, covariates.nuisance=list(slope=slope),
#'                 bws=bws, verbose=TRUE)
#' out2
#'
#' @export
#'
tau.est <- function(X, covariate.interest, covariates.nuisance, bws, n.test.points=1000,
                    nonparametric=TRUE, verbose=FALSE, bw.factor.rhonhat=1){

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
    nm <- "tau"
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
    nm <- "tau_p"

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

  ### Estimated value of the Kendall's correlation coefficient for the covariate of interest
  res <- cor(srf[test.points],covariate.interest[test.points],method="kendall")
  names(res) <- nm

  if (verbose){
    cat("Observed value of the Kendall's correlation coefficient: ")
    cat(res)
    cat("\n")
  }

  return(res)
}
