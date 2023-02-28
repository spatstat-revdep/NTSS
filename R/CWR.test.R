#' Random shift test of independence between a point process and a covariate
#'
#' @description Nonparametric test of independence between a point process and a random field
#' (covariate of interest), taking into account the possible effect of nuisance covariates,
#' see Dvořák and Mrkvička (2022).
#' The test is based on random shifts. Either the torus correction or the variance
#' correction can be used, see Mrkvička et al. (2021).
#'
#' @details The test statistic is the covariate-weighted residual measure of the observation window, see the paper
#' Dvořák and Mrkvička (2022). If no nuisance covariates are given, the null model
#' assumes a constant intensity function of the point process. If one or more nuisance covariates are
#' provided, the null model assumes an intensity function depending on the nuisance
#' covariates (but not on the covariate of interest) and the residuals are constructed
#' using this intensity function.
#'
#' The residuals can be constructed in a nonparametric way (see Baddeley et al. (2012))
#' or in a parametric way (assuming a log-linear form of the intensity function
#' and using the \code{ppm} function from the \code{spatstat} package,
#' see Baddeley et al. (2015)). This choice is given by the argument \code{nonparametric}.
#' The nonparametric residuals are recommended if the log-linear form of the intensity function
#' is not clearly justified.
#' Also, different types of residuals can be considered (raw, Pearson or inverse,
#' see Baddeley et al. (2015)). This choice is given by the argument \code{type}.
#'
#' The torus correction can be applied for rectangular windows. On the other hand,
#' the variance correction is applicable both for rectangular and for irregular windows.
#' The choice of the correction is given by the argument \code{correction}.
#' Based on the simulation studies in Dvořák and Mrkvička (2022),
#' the torus and variance corrections perform almost equally well for the CWR test.
#' Hence, the torus correction is recommended for its smaller computational cost.
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
#' @param type which type of residuals should be used when computing the test statistic (possible choices are "raw", "Pearson" and "inverse")
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
#' # test with no nuisance covariates, with only 99 shifts to speed up the computation
#' out1 <- CWR.test(X, covariate.interest=elevation, covariates.nuisance=NULL, N.shifts = 99,
#'                  verbose=TRUE, correction="torus", radius=250)
#' out1
#'
#' # test with one nuisance covariate, with only 99 shifts to speed up the computation
#' out2 <- CWR.test(X, covariate.interest=elevation, covariates.nuisance=list(slope=slope),
#'                  N.shifts = 99, verbose=TRUE, correction="torus", radius=250)
#' out2
#'
#' @export
#'
CWR.test <- function(X, covariate.interest, covariates.nuisance, N.shifts=999, radius, correction,
                     type="raw", nonparametric=TRUE, verbose=FALSE, bw.factor.rhonhat=1){

  if (verbose){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  ### Estimation of intensity function

  ncovs <- length(covariates.nuisance)

  if (ncovs == 0){ # no covariates, use homogenenous model

    fit <- ppm(X ~ 1)
    est.rho <- predict(fit, dimyx=covariate.interest$dim, window=as.owin(covariate.interest))
    if (verbose){cat("Using model with constant intensity function for construction of residuals.\n")}

  } else { # at least one nuisance covariate present

    if (nonparametric){ # nonparametric estimate of intensity function
      if (verbose){cat("Using nonparametric estimate of intensity function for construction of residuals.\n")}
      est.rho <- rhonhat(X, covariates.nuisance, prediction.grid.from=covariate.interest, bw.factor=bw.factor.rhonhat)
    } else { # parametric (log-linear) estimate of intensity function
      if (verbose){cat("Using parametric (log-linear) estimate of intensity function for construction of residuals.\n")}
      fit <- ppm(X ~ ., covariates=covariates.nuisance)
      est.rho <- predict(fit, dimyx=covariate.interest$dim, window=as.owin(covariate.interest))
    }
  }

  if (verbose & (ncovs >= 1)){
    par(mar=c(1/2,1/2,2,2), mfrow=c(1,1))
    plot(est.rho, main="Estimated intensity function")
  }


  ### CWR test with torus correction

  if (correction=="torus"){
    values.simulated <- rep(NA, times=N.shifts+1)
    values.simulated[1] <- switch(type,
                                  raw     = sum(covariate.interest[X]) - summary(covariate.interest*est.rho)$integral,
                                  Pearson = sum(covariate.interest[X]/sqrt(est.rho[X])) - summary(covariate.interest*sqrt(est.rho))$integral,
                                  inverse = sum(covariate.interest[X]/est.rho[X]) - summary(covariate.interest)$integral)

    test.stat <- values.simulated[1]
    names(test.stat) <- "CWR"

    correct <- "torus"
    names(correct) <- "correction"

    for (k in 1:N.shifts){
      jump <- runifdisc(1, radius = radius)
      n.shift.x <- sum(covariate.interest$xcol < covariate.interest$xrange[1] + (jump$x %% diff(covariate.interest$xrange))) %% covariate.interest$dim[2]
      n.shift.y <- sum(covariate.interest$yrow < covariate.interest$yrange[1] + (jump$y %% diff(covariate.interest$yrange))) %% covariate.interest$dim[1]
      covariate.interest.shift <- covariate.interest
      if (n.shift.x > 0){
        covariate.interest.shift$v <- covariate.interest.shift$v[,c((covariate.interest$dim[2]-n.shift.x+1):covariate.interest$dim[2],1:(covariate.interest$dim[2]-n.shift.x))]
      }
      if (n.shift.y > 0){
        covariate.interest.shift$v <- covariate.interest.shift$v[c((covariate.interest$dim[1]-n.shift.y+1):covariate.interest$dim[1],1:(covariate.interest$dim[1]-n.shift.y)),]
      }
      suppressWarnings(values.simulated[k+1] <- switch(type,
                                                       raw     = sum(covariate.interest.shift[X]) - summary(covariate.interest.shift*est.rho)$integral,
                                                       Pearson = sum(covariate.interest.shift[X]/sqrt(est.rho[X])) - summary(covariate.interest.shift*sqrt(est.rho))$integral,
                                                       inverse = sum(covariate.interest.shift[X]/est.rho[X]) - summary(covariate.interest.shift)$integral))
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


  ### CWR test with variance correction

  if (correction == "variance"){

    values.simulated <- rep(NA, times=N.shifts+1)
    values.simulated[1] <- switch(type,
                                  raw     = sum(covariate.interest[X]) - summary(covariate.interest*est.rho)$integral,
                                  Pearson = sum(covariate.interest[X]/sqrt(est.rho[X])) - summary(covariate.interest*sqrt(est.rho))$integral,
                                  inverse = sum(covariate.interest[X]/est.rho[X]) - summary(covariate.interest)$integral)
    n.simulated <- rep(NA, times=N.shifts+1)
    n.simulated[1] <- X$n

    for (k in 1:N.shifts){
      jump <- runifdisc(1, radius = radius)
      jump <- c(jump$x,jump$y)
      jump[1] <- jump[1] - jump[1] %% covariate.interest$xstep
      jump[2] <- jump[2] - jump[2] %% covariate.interest$ystep

      covariate.interest.shift <- shift.im(covariate.interest, jump)
      W.shift <- shift.owin(X$window, jump)
      W.new <- intersect.owin(X$window, W.shift)
      X.new <- X[W.new]
      est.rho.new <- est.rho[W.new]
      covariate.interest.new <- covariate.interest.shift[W.new, drop=FALSE]

      suppressWarnings(values.simulated[k+1] <- switch(type,
                                                       raw     = sum(covariate.interest.new[X.new]) - summary(covariate.interest.new*est.rho.new)$integral,
                                                       Pearson = sum(covariate.interest.new[X.new]/sqrt(est.rho.new[X.new])) - summary(covariate.interest.new*sqrt(est.rho.new))$integral,
                                                       inverse = sum(covariate.interest.new[X.new]/est.rho.new[X.new]) - summary(covariate.interest.new)$integral))
      n.simulated[k+1] <- X.new$n
    }
    values.std <- (values.simulated - mean(values.simulated))/sqrt(n.simulated)
    test.rank <- rank(values.std)[1]
    pval <- 2*min(test.rank, N.shifts+1-test.rank)/(N.shifts+1)

    test.stat <- values.std[1]
    names(test.stat) <- "CWR_std"

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
  testname <- "Random shift test of independence between a point process and a covariate, with nuisance covariates"
  alternative <- "two-sided"
  ress <- switch(type,
                      raw = "raw",
                      Pearson    = "Pearson's",
                      inverse    = "inverse")

  # Note: is there a reasonable way how to include the name of the list of nuisance covariates, too?
  dn <- paste(substitute(X), "and", substitute(covariate.interest))

  result <- structure(list(statistic = test.stat, parameter = list(N.shifts=res.N.shifts,correction=correct,residuals=ress),
                           p.value = pval, method = testname, data.name = dn,
                           alternative = alternative),
                      class = "htest")

  return(result)
}



rhonhat <- function(X, covariates, bw.factor=1, plot=FALSE, prediction.grid.from=NULL){
  # Nonparametric estimation of the intensity depending on covariates
  # based on the paper Baddeley et al. (2012): Nonparametric estimation
  # of the dependence of a spatial point process on spatial covariates.

  # The implemented version is the 'ratio' estimator (8) from the paper.

  # This function provides predicted intensity function in the pixel
  # grid determined from the argument "prediction.grid.from" or
  # from the first covariate in the list "covariates".

  # X ... a point pattern
  # covariates ... a list of covariates (of the "im" format from the "spatstat" package)
  # bw.factor ... a positive real number, multiplicative factor in the bandwidth matrix
  # plot ... logical, should the estimated intensity function be plotted?
  # prediction.grid.from ... image (of the "im" format from the "spatstat" package) determining
  #                          the pixel grid on which the intensity function will be estimated;
  #                          if NULL, the first covariate from the "covariate" list will be used

  # Output value: pixel image of predicted intensity function

  # Number of covariates
  ncovs <- length(covariates)
  if (ncovs <= 4){ # for speeding up computation, implemented up to 4 dimensions
    binned <- TRUE
  } else {
    binned <- FALSE
  }

  # Auxiliary values
  if (is.null(prediction.grid.from)){prediction.grid.from <- covariates[[1]]}
  dim1 <- prediction.grid.from$dim[1]
  dim2 <- prediction.grid.from$dim[2]
  grid.pts.Wspace <- suppressWarnings(ppp(x=rep(prediction.grid.from$xcol, times=dim1), y=rep(prediction.grid.from$yrow, each=dim2), window=X$window))
  # grid.pts.Cspace <- matrix(NA, ncol=ncovs, nrow=dim1*dim2)
  grid.pts.Cspace <- matrix(NA, ncol=ncovs, nrow=grid.pts.Wspace$n)
  for (i in 1:ncovs){
    grid.pts.Cspace[,i] <- safelookup(covariates[[i]], grid.pts.Wspace)
  }
  x <- matrix(NA, ncol=ncovs, nrow=X$n)
  for (i in 1:ncovs){
    x[,i] <- safelookup(covariates[[i]], X)
  }

  # Scott's rule for choosing bandwidth
  if (ncovs > 1){
    H.aux <- diag(ncovs)
    for (i in 1:ncovs){
      H.aux[i,i] <- sd(x[,i])
    }
    H <- bw.factor*(H.aux*(X$n)^(-1/(ncovs+4)))^2
  } else {
    h <- bw.factor*(sd(x[,1])*(X$n)^(-1/(ncovs+4)))^2
  }

  # Kernel estimation
  if (ncovs > 1){
    fstar <- kde(x=x, eval.points=grid.pts.Cspace, H=H, binned=binned)$estimate
    gstar <- kde(x=grid.pts.Cspace, eval.points=grid.pts.Cspace, H=H, binned=binned)$estimate
  } else {
    fstar <- kde(x=x, eval.points=grid.pts.Cspace, h=h, binned=binned)$estimate
    gstar <- kde(x=grid.pts.Cspace, eval.points=grid.pts.Cspace, h=h, binned=binned)$estimate
  }
  pred.rho <- (X$n/area(X$window))*fstar/gstar
  est.rho <- prediction.grid.from

  if (grid.pts.Wspace$n == dim1*dim2){
    est.rho$v <- matrix(pred.rho, nrow=dim1, byrow=TRUE)
  } else {
    est.rho$v <- matrix(NA, nrow=dim1, ncol=dim2)
    aux.at <- 1
    aux.count <- 1
    aux.inside <- inside.owin(x=rep(prediction.grid.from$xcol, times=dim1), y=rep(prediction.grid.from$yrow, each=dim2), w=X$window)

    for (aux.i in 1:dim1){
      for (aux.j in 1:dim2){
        if (aux.inside[aux.at]){
          est.rho$v[aux.i,aux.j] <- pred.rho[aux.count]
          aux.count <- aux.count + 1
        }
        aux.at <- aux.at + 1
      }
    }
  }

  # Visualisation
  if (plot){plot(est.rho, main="Predicted intensity function")}

  # Output (pixel image of predicted intensity function)
  return(est.rho)
}
