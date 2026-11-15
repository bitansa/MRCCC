#' @include MRCCCFit-methods.R
NULL

#' Posterior sign probabilities of the interaction effect
#'
#' Computes the posterior probability that the ligand-by-receptor
#' interaction \eqn{\beta_{XZ}} is positive and that it is negative, both
#' marginally and conditionally on the causal block being included
#' (\eqn{\gamma = 1}).
#'
#' @param x An [MRCCCFit-class] object returned by [mrccc()].
#' @param ... Not used.
#'
#' @details
#' Under the spike-and-slab prior, the marginal posterior of
#' \eqn{\beta_{XZ}} mixes draws taken while \eqn{\gamma = 0}, which are
#' concentrated near zero and carry no directional information, with slab
#' draws taken while \eqn{\gamma = 1}.
#'
#' * The **marginal** probabilities `pPositive` and `pNegative` use every
#'   retained draw. They combine the evidence that there is an effect with
#'   the evidence about its direction:
#'   \eqn{P(\beta_{XZ} > 0) = \mathrm{PIP} \cdot P(\beta_{XZ} > 0 \mid
#'   \gamma = 1) + (1 - \mathrm{PIP}) \cdot P(\beta_{XZ} > 0 \mid
#'   \gamma = 0)}, and the last factor is close to one half.
#' * The **conditional** probabilities `pPositiveIncluded` and
#'   `pNegativeIncluded` use only the draws with \eqn{\gamma = 1}. They
#'   describe the direction of the interaction given that the ligand acts
#'   on the pathway. They are `NA` when no draw has \eqn{\gamma = 1}.
#'
#' @return A one-row data frame with columns `pip`, `pPositive`,
#'   `pNegative`, `pPositiveIncluded`, `pNegativeIncluded` and `nIncluded`
#'   (the number of pooled draws with \eqn{\gamma = 1}).
#'
#' @seealso [mrccc()], [credibleIntervals()], [signReversal()]
#'
#' @examples
#' sim <- simulateMRCCC(n = 200, seed = 1)
#' dat <- buildTriplet(sim$sender, sim$receiver, sim$genotypes,
#'                     ligand = "LIG1", receptor = "REC1",
#'                     pathway = sim$pathwayGenes,
#'                     covariates = c("sex", "age", "PC1", "PC2"))
#' fit <- mrccc(dat, nIter = 4000, burnIn = 1000, seed = 1)
#' signProbabilities(fit)
#'
#' @rdname signProbabilities
#' @export
setMethod("signProbabilities", "MRCCCFit", function(x, ...) {
    bxz <- x@draws$beta_XZ
    inc <- x@draws$gamma == 1
    nInc <- sum(inc)
    data.frame(
        pip = x@pip,
        pPositive = mean(bxz > 0),
        pNegative = mean(bxz < 0),
        pPositiveIncluded = if (nInc > 0L) mean(bxz[inc] > 0) else NA_real_,
        pNegativeIncluded = if (nInc > 0L) mean(bxz[inc] < 0) else NA_real_,
        nIncluded = as.integer(nInc),
        stringsAsFactors = FALSE)
})
