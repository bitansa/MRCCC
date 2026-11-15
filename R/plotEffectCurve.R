#' @include MRCCCFit-methods.R credibleIntervals.R
NULL

#' Plot the receptor-modulated ligand effect with a credible band
#'
#' Draws the posterior mean of the receptor-modulated ligand effect
#' \eqn{\beta_X + \beta_{XZ} z} over the observed range of receptor
#' expression, with a pointwise credible band from paired posterior draws.
#'
#' @param x An [MRCCCFit-class] object returned by [mrccc()].
#' @param level Credible level in \eqn{(0, 1)} of the pointwise band.
#' @param nDraws Maximum number of paired draws used for the band. A
#'   systematic subsample of the pooled draws of this size is taken.
#' @param standardized Logical. If `TRUE` (the default) the effect is on the
#'   standardised scale against receptor expression in standard deviations
#'   from its mean. If `FALSE`, it is on the centred input scale.
#' @param ... Not used.
#'
#' @details
#' At each receptor level \eqn{z} the band is the equal-tailed interval of
#' \eqn{\beta_X + \beta_{XZ} z} across paired draws, which accounts for the
#' posterior correlation of the two coefficients. The draws include both
#' spike and slab iterations. A dashed line marks an effect of zero. Where
#' the band crosses it, the direction of the ligand effect is uncertain.
#'
#' @return A `ggplot` object.
#'
#' @seealso [mrccc()], [signReversal()]
#'
#' @examples
#' sim <- simulateMRCCC(n = 200, betaX = 0.3, betaXZ = -0.3, seed = 1)
#' dat <- buildTriplet(sim$sender, sim$receiver, sim$genotypes,
#'                     ligand = "LIG1", receptor = "REC1",
#'                     pathway = sim$pathwayGenes,
#'                     covariates = c("sex", "age", "PC1", "PC2"))
#' fit <- mrccc(dat, nIter = 4000, burnIn = 1000, seed = 1)
#' plotEffectCurve(fit)
#'
#' @rdname plotEffectCurve
#' @export
setMethod("plotEffectCurve", "MRCCCFit",
          function(x, level = 0.95, nDraws = 2000L, standardized = TRUE,
                   ...) {
    .checkLevel(level)
    check_scalar(nDraws, "nDraws")
    if (nDraws < 2 || nDraws != round(nDraws)) {
        stop("'nDraws' must be an integer of at least 2.", call. = FALSE)
    }
    .checkFlag(standardized, "standardized")

    zRng <- x@scale$zRangeStd
    if (!standardized) zRng <- zRng * x@scale$sdZ
    sf <- if (standardized) std_factors(x) else
        c(beta_X = 1, beta_XZ = 1, beta_Z = 1)
    bx <- x@draws$beta_X * sf[["beta_X"]]
    bxz <- x@draws$beta_XZ * sf[["beta_XZ"]]
    if (length(bx) > nDraws) {
        idx <- round(seq(1, length(bx), length.out = nDraws))
        meanBx <- mean(bx)
        meanBxz <- mean(bxz)
        bx <- bx[idx]
        bxz <- bxz[idx]
    } else {
        meanBx <- mean(bx)
        meanBxz <- mean(bxz)
    }

    zg <- seq(zRng[1L], zRng[2L], length.out = 101L)
    curves <- outer(bx, rep(1, length(zg))) + outer(bxz, zg)
    probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
    band <- apply(curves, 2L, stats::quantile, probs = probs, names = FALSE)
    df <- data.frame(z = zg, effect = meanBx + meanBxz * zg,
                     lower = band[1L, ], upper = band[2L, ])

    lab <- x@labels
    recName <- if (is.na(lab[["receptor"]])) "Receptor" else lab[["receptor"]]
    ligName <- if (is.na(lab[["ligand"]])) "ligand" else lab[["ligand"]]
    xLab <- if (standardized) {
        paste0(recName, " expression (SD from mean)")
    } else {
        paste0(recName, " expression (centred)")
    }
    yLab <- if (standardized) {
        paste0("Effect of ", ligName, " (SD of pathway per SD)")
    } else {
        paste0("Effect of ", ligName, " (centred scale)")
    }
    ttl <- .tripletTitle(lab)

    ggplot2::ggplot(df, ggplot2::aes(x = .data$z)) +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
            fill = "grey70", alpha = 0.5) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                            colour = "grey40") +
        ggplot2::geom_line(ggplot2::aes(y = .data$effect), colour = "black",
                           linewidth = 0.8) +
        ggplot2::labs(
            x = xLab, y = yLab,
            title = if (is.null(ttl)) {
                "Receptor-modulated ligand effect"
            } else {
                ttl
            },
            subtitle = sprintf(
                "Posterior mean with %d%% pointwise credible band",
                as.integer(round(100 * level)))) +
        ggplot2::theme_bw()
})
