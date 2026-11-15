#' @include MRCCCFit-methods.R credibleIntervals.R
NULL

#' Receptor level at which the ligand effect reverses sign
#'
#' Computes the posterior distribution of the sign-reversal threshold
#' \eqn{\tau = -\beta_X / \beta_{XZ}}, the receptor expression level at
#' which the receptor-modulated ligand effect \eqn{\beta_X + \beta_{XZ} z}
#' changes sign, from paired posterior draws.
#'
#' @param x An [MRCCCFit-class] object returned by [mrccc()].
#' @param level Credible level in \eqn{(0, 1)} for the reported interval.
#' @param ... Not used.
#'
#' @details
#' \eqn{\tau} is a ratio of two correlated parameters, so its posterior is
#' computed from paired draws of \eqn{(\beta_X, \beta_{XZ})}. Draws whose
#' standardised interaction is negligible,
#' \eqn{|\beta_{XZ}|\,\mathrm{sd}(X)\,\mathrm{sd}(Z)/\mathrm{sd}(Y) \le
#' 10^{-8}}, are dropped, since \eqn{\tau} diverges there. All remaining
#' draws are used, including those from the spike, so for a fit with a PIP
#' well below one the threshold is best read together with the PIP.
#'
#' **Scale.** Receptor expression is centred before fitting, so \eqn{\tau}
#' is reported in standard deviations of receptor expression relative to
#' its mean over donors: \eqn{\tau = 0} is a donor with average receptor
#' expression and \eqn{\tau = 1} a donor one standard deviation above
#' average. The observed receptor range, on the same scale, is returned as
#' `zMin` and `zMax`, and `pInRange` is the posterior probability that
#' \eqn{\tau} lies inside it, which determines whether the reversal is
#' realised among the donors. The posterior of a ratio is heavy-tailed, so
#' the median and quantiles are reported rather than a mean. Note that the
#' interaction acts on the genetically predicted receptor level, whereas
#' the range refers to observed receptor expression.
#'
#' @return A one-row data frame with columns `tauMedian`, `tauLower`,
#'   `tauUpper` (equal-tailed credible bounds at `level`), `level`,
#'   `nDraws` (paired draws used), `zMin`, `zMax` and `pInRange`. The
#'   numeric summaries are `NA` when no usable draw is available.
#'
#' @seealso [mrccc()], [plotEffectCurve()], [signProbabilities()]
#'
#' @examples
#' sim <- simulateMRCCC(n = 200, betaX = 0.3, betaXZ = -0.3, seed = 1)
#' dat <- buildTriplet(sim$sender, sim$receiver, sim$genotypes,
#'                     ligand = "LIG1", receptor = "REC1",
#'                     pathway = sim$pathwayGenes,
#'                     covariates = c("sex", "age", "PC1", "PC2"))
#' fit <- mrccc(dat, nIter = 4000, burnIn = 1000, seed = 1)
#' signReversal(fit)
#'
#' @rdname signReversal
#' @export
setMethod("signReversal", "MRCCCFit", function(x, level = 0.95, ...) {
    .checkLevel(level)
    zr <- x@scale$zRangeStd
    tau <- .tauDraws(x)
    if (length(tau) == 0L) {
        return(data.frame(tauMedian = NA_real_, tauLower = NA_real_,
                          tauUpper = NA_real_, level = level, nDraws = 0L,
                          zMin = zr[1L], zMax = zr[2L], pInRange = NA_real_))
    }
    probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
    q <- stats::quantile(tau, probs = probs, names = FALSE)
    data.frame(tauMedian = q[2L], tauLower = q[1L], tauUpper = q[3L],
               level = level, nDraws = length(tau), zMin = zr[1L],
               zMax = zr[2L], pInRange = mean(tau >= zr[1L] & tau <= zr[2L]))
})

# Paired draws of tau on the standardised receptor scale.
.tauDraws <- function(x) {
    bx <- x@draws$beta_X
    bxz <- x@draws$beta_XZ
    bxzStd <- bxz * std_factors(x)[["beta_XZ"]]
    ok <- is.finite(bx) & is.finite(bxzStd) & abs(bxzStd) > 1e-8
    (-bx[ok] / bxz[ok]) / x@scale$sdZ
}
