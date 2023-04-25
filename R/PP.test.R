#' Random shift test of independence in a bivariate point process
#'
#' @description Nonparametric test of independence between a pair of point processes based on random shifts.
#' Either the torus correction or the variance correction can be used (note that the
#' variance correction is not yet implemented in this package but the corresponding source
#' codes can be obtained from the authors upon request).
#'
#' @details The test statistic is either the cross K-function K12 or the expectation of the
#' cross nearest-neighbor distance ED12, see the paper Mrkvička et al. (2021). It is recommended to use the K12 statistic
#' for regular and Poisson processes, while it is recommended to use ED12 for cluster processes due to the potential
#' liberality of K12 in this case.
#'
#' The two observed point patterns to be tested should be supplied using the
#' arguments \code{X} and \code{Y}. The shift vectors are generated from the
#' uniform distribution on the disk with radius given by the argument \code{radius}
#' and centered in the origin.
#'
#' The test statistic is determined by the argument \code{statistic}. For "K12", the test statistic
#' is functional, the stationary cross K-function. The range of arguments is from 0 to \code{rmax.K} (there
#' is a sensible default). The outcome of the test is determined by the global envelope test in this case.
#' Alternatively, the \code{statistic} can be set to "ED12", meaning that the test statistic is scalar,
#' the expectation of the cross nearest-neighbor distance. The outcome of the test is determined by the
#' classical univariate Monte Carlo test in this case.
#'
#' The argument \code{verbose} determines if auxiliary information and plots should be provided.
#'
#' @import spatstat
#' @import spatstat.explore
#' @import stats
#' @import GET
#'
#' @param X first point pattern dataset (object of class \code{ppp})
#' @param Y second point pattern dataset (object of class \code{ppp})
#' @param N.shifts integer, how many random shifts should be performed in the random shift test
#' @param radius positive real number determining the radius of the disk on which the shift vectors are uniformly distributed
#' @param correction which correction should be applied in the random shift test (possible choices are "torus" and "variance")
#' @param statistic which test statistic should be used (possible choices are "K12" and "ED12")
#' @param rmax.K positive real number, for the cross K-function determines the maximum argument to be considered
#' @param verbose logical value indicating whether auxiliary information should be printed and auxiliary figures plotted during the computation
#'
#' @return The p-value of the random shift test of independence between a pair of point processes.
#'
#' @references T. Mrkvička, J. Dvořák, J.A. González, J. Mateu (2021): Revisiting the random shift approach for testing in spatial statistics. Spatial Statistics 42, 100430.
#'
#' @examples
#'
#' library(spatstat)
#' library(GET)
#'
#' set.seed(123)
#'
#' X <- rpoispp(150)
#' Y <- rpoispp(150)
#'
#' # tests run with only 99 shifts to speed up the computation
#' out1 <- PP.test(X=X, Y=Y, N.shifts = 99, radius=0.5, statistic="K12",
#'                 correction="torus", verbose=TRUE)
#' out1
#' plot(out1$GET.outcome)
#'
#' out2 <- PP.test(X=X, Y=Y, N.shifts = 99, radius=0.5, statistic="ED12",
#'                 correction="torus", verbose=TRUE)
#' out2
#'
#' @export
#'
PP.test <- function(X, Y, N.shifts=999, radius, correction, statistic, rmax.K=NULL, verbose=FALSE){

  if (verbose){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }

  # Z ... bivariate point pattern to be tested, with marks "X" and "Y"
  Z <- superimpose(X = X, Y = Y)

  if ((statistic != "ED12") & (statistic != "K12")){
    print("The argument statistic must be either ED12 or K12, determining the choice of the test statistic.")
    return(NULL)
  }


  if (statistic=="ED12"){

    ### P-P test with torus correction

    if (correction=="torus"){

      values.simulated <- rep(NA, times=N.shifts+1)
      fun <- Gcross(Z, correction	=c("km"))
      values.simulated[1] <- stieltjes(function(x){x}, fun)$km

      test.stat.value <- values.simulated[1]
      names(test.stat.value) <- "ED_12"

      correct <- "torus"
      names(correct) <- "correction"
      GET.outcome <- NULL

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
      print("Variance correction for this test is not yet implemented but the corresponding source codes can be obtained from the authors upon request.")
      return(NULL)
    }
  }

  if (statistic=="K12"){

    ### P-P test with torus correction

    if (correction=="torus"){

      correct <- "torus"
      names(correct) <- "correction"

      if (is.null(rmax.K)){rmax.K <- rmax.rule("K", X$window, intensity(X))} # default from the spatstat package

      K.trans <- function(Z, rr = FALSE){
        # Estimate cross K-function with translation edge-correction
        # Z .... marked point pattern with marks "X" and "Y" (bivariate point pattern)
        # rr ... logical, should also the vector of argument values be returned?

        K.hat <- Kcross(Z, "X", "Y", r = seq(0, rmax.K, length.out = 513), correction = "translation")
        if (rr) return(list(r = K.hat$r, vals = K.hat$trans)) else return(K.hat$trans)
      }

      K.observed <- K.trans(Z, rr = TRUE)
      test.stat.value <- "K_12"
      names(test.stat.value) <- "K_12"

      r <- K.observed$r
      K.observed <- K.observed$vals

      RS.torus <- function() {
        Z.shift <- rshift.ppp(Z, which='X', radius=radius, edge="torus")
        Kreduced <- K.trans(Z.shift)
        return(Kreduced)
      }

      K.simulated <- replicate(N.shifts, RS.torus())

      # Global envelope test
      CS <- create_curve_set(list(r = r, obs = K.observed, sim_m = K.simulated))
      GET.outcome <- global_envelope_test(CS, type = "erl")
      pval <-  attr(GET.outcome, "p")


      if (verbose){
        cat("Random shift with torus correction, p-value: ")
        cat(pval)
        cat("\n")
        print(GET.outcome)
      }
    }


    ### P-P test with variance correction

    if (correction == "variance"){
      print("Variance correction for this test is not implemented here due to the necessity of high-dimensional numerical integration.")
      return(NULL)
    }
  }

  res.N.shifts <- N.shifts
  names(res.N.shifts) <- "N.shifts"
  testname <- "Random shift test of independence between a pair of point processes"
  alternative <- "two-sided"

  result <- structure(list(statistic = test.stat.value, parameter = list(N.shifts=res.N.shifts,correction=correct),
                           p.value = pval, method = testname, data.name = paste(substitute(X), "and", substitute(Y)),
                           alternative = alternative, GET.outcome=GET.outcome),
                      class = "htest")

  return(result)
}
