#' @include MRCCCFit-methods.R
NULL

#' Credible intervals for the MR-CCC effect sizes
#'
#' Computes equal-tailed posterior credible intervals for \eqn{\beta_X},
#' \eqn{\beta_{XZ}} and \eqn{\beta_Z} from the pooled retained draws of an
#' MR-CCC fit.
#'
#' @param x An [MRCCCFit-class] object returned by [mrccc()].
#' @param level Credible level in \eqn{(0, 1)}. The default 0.95 gives the
#'   2.5 and 97.5 percent posterior quantiles.
#' @param standardized Logical. If `TRUE` (the default) the intervals are on
#'   the standardised scale (per standard deviation of ligand and receptor
#'   expression, in standard deviations of pathway activity). If `FALSE`,
#'   they are on the centred input scale.
#' @param ... Not used.
#'
#' @details
#' The prior on \eqn{(\beta_X, \beta_{XZ})} is a spike-and-slab mixture, so
#' the marginal posterior of each of these coefficients mixes draws
#' concentrated near zero (\eqn{\gamma = 0}) with slab draws
#' (\eqn{\gamma = 1}). An interval that covers zero therefore reflects the
#' posterior weight of the spike as much as the uncertainty within the slab.
#' See [signProbabilities()] for a decomposition.
#'
#' @return A data frame with columns `parameter` (`beta_X`, `beta_XZ`,
#'   `beta_Z`), `lower`, `median` and `upper`.
#'
#' @seealso [mrccc()], [effectSizes()], [signProbabilities()]
#'
#' @examples
#' sim <- simulateMRCCC(n = 200, seed = 1)
#' dat <- buildTriplet(sim$sender, sim$receiver, sim$genotypes,
#'                     ligand = "LIG1", receptor = "REC1",
#'                     pathway = sim$pathwayGenes,
#'                     covariates = c("sex", "age", "PC1", "PC2"))
#' fit <- mrccc(dat, nIter = 4000, burnIn = 1000, seed = 1)
#' credibleIntervals(fit)
#' credibleIntervals(fit, level = 0.9, standardized = FALSE)
#'
#' @rdname credibleIntervals
#' @export
setMethod("credibleIntervals", "MRCCCFit",
          function(x, level = 0.95, standardized = TRUE, ...) {
    .checkLevel(level)
    .checkFlag(standardized, "standardized")
    pars <- c("beta_X", "beta_XZ", "beta_Z")
    sf <- if (standardized) std_factors(x) else
        c(beta_X = 1, beta_XZ = 1, beta_Z = 1)
    probs <- c((1 - level) / 2, 0.5, 1 - (1 - level) / 2)
    q <- vapply(pars, function(p) {
        stats::quantile(x@draws[[p]] * sf[[p]], probs = probs, names = FALSE)
    }, numeric(3))
    data.frame(parameter = pars, lower = unname(q[1L, ]),
               median = unname(q[2L, ]), upper = unname(q[3L, ]),
               row.names = NULL, stringsAsFactors = FALSE)
})

.checkLevel <- function(level) {
    check_scalar(level, "level")
    if (level <= 0 || level >= 1) {
        stop("'level' must lie strictly between 0 and 1.", call. = FALSE)
    }
    invisible(level)
}
