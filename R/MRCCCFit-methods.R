#' @include AllClasses.R AllGenerics.R mrccc.R
NULL

#' Accessors for MRCCCFit objects
#'
#' Extract results from a fitted MR-CCC model.
#'
#' @param x An [MRCCCFit-class] object returned by [mrccc()].
#' @param standardized Logical. If `TRUE` (the default for `effectSizes()`)
#'   coefficients are on the standardised scale (per standard deviation of
#'   ligand and receptor expression, in standard deviations of pathway
#'   activity). If `FALSE` they are on the centred input scale.
#' @param ... Not used.
#'
#' @return
#' * `pip()`: the posterior inclusion probability, a number in \[0, 1\].
#' * `pipMCSE()`: a named numeric vector with the Monte Carlo standard error
#'   (`mcse`) and effective sample size (`ess`) of the PIP.
#' * `effectSizes()`: a data frame with columns `parameter` and `mean`, the
#'   posterior mean. The parameters are `beta_X`, the ligand effect at the
#'   mean receptor level, `beta_XZ`, the ligand-by-receptor interaction, and
#'   `beta_Z`, the receptor main effect.
#' * `convergence()`: a data frame of Gelman-Rubin statistics (`rhat`) for
#'   `beta_X`, `beta_XZ`, `gamma` and the joint log-likelihood. The values
#'   are `NA` for a single chain.
#' * `instrumentStrength()`: a data frame with one row per exposure
#'   (`ligand`, `receptor`) and columns `nInstruments`, `F` (first-stage
#'   partial F statistic after the covariates) and `weak` (`F < 10`).
#' * `posteriorDraws()`: a data frame of the pooled retained draws with
#'   columns `chain`, `beta_X`, `beta_XZ`, `beta_Z`, `gamma`, `mu` and
#'   `loglik`. With `standardized = TRUE` the three coefficients are
#'   rescaled.
#' * `firstStageEffects()`: a data frame with one row per instrument and
#'   columns `exposure` (`ligand` or `receptor`), `snp`, `mean`, `lower` and
#'   `upper`: the posterior mean and 95% credible interval of the effect of
#'   one copy of the alternative allele on expression.
#' * `posteriorMeans()`: a list with the posterior means of the remaining
#'   parameters: `covariateEffects` (data frame with columns `covariate`,
#'   `ligand`, `receptor` and `pathway`), `residualVariances`,
#'   `inclusionProbabilityPrior` and `intercept`.
#' * `inputScale()`: a list with the standard deviations (`sdX`, `sdZ`,
#'   `sdY`) and means (`meanX`, `meanZ`, `meanY`) of ligand expression,
#'   receptor expression and pathway activity, and the observed receptor
#'   range in standard deviations from its mean (`zRangeStd`).
#' * `fitSettings()`: a list of the run settings (`nIter`, `burnIn`,
#'   `thin`, `nChains`, `nKeep`, `pipThreshold`, `hyper`, `g`, `initScale`
#'   and `priorScale`).
#' * `tripletLabels()`: the labels of the triplet.
#'
#' @seealso [mrccc()], [credibleIntervals()], [signProbabilities()],
#'   [signReversal()], [plotEffectCurve()]
#'
#' @examples
#' sim <- simulateMRCCC(n = 200, seed = 1)
#' dat <- buildTriplet(sim$sender, sim$receiver, sim$genotypes,
#'                     ligand = "LIG1", receptor = "REC1",
#'                     pathway = sim$pathwayGenes,
#'                     covariates = c("sex", "age", "PC1", "PC2"))
#' fit <- mrccc(dat, nIter = 4000, burnIn = 1000, seed = 1)
#' pip(fit)
#' pipMCSE(fit)
#' effectSizes(fit)
#' convergence(fit)
#' instrumentStrength(fit)
#' head(posteriorDraws(fit))
#' firstStageEffects(fit)
#' posteriorMeans(fit)$residualVariances
#' fitSettings(fit)$nKeep
#'
#' @name MRCCCFit-accessors
#' @rdname MRCCCFit-accessors
NULL

#' @rdname MRCCCFit-accessors
#' @export
setMethod("pip", "MRCCCFit", function(x) x@pip)

#' @rdname MRCCCFit-accessors
#' @export
setMethod("pipMCSE", "MRCCCFit",
          function(x) c(mcse = x@pipMCSE, ess = x@pipESS))

#' @rdname MRCCCFit-accessors
#' @export
setMethod("effectSizes", "MRCCCFit", function(x, standardized = TRUE, ...) {
    .checkFlag(standardized, "standardized")
    pars <- c("beta_X", "beta_XZ", "beta_Z")
    sf <- if (standardized) std_factors(x) else
        c(beta_X = 1, beta_XZ = 1, beta_Z = 1)
    data.frame(parameter = pars,
               mean = vapply(pars, function(p) mean(x@draws[[p]]) * sf[[p]],
                             numeric(1), USE.NAMES = FALSE),
               stringsAsFactors = FALSE)
})

#' @rdname MRCCCFit-accessors
#' @export
setMethod("convergence", "MRCCCFit", function(x) x@diagnostics)

#' @rdname MRCCCFit-accessors
#' @export
setMethod("instrumentStrength", "MRCCCFit", function(x) x@instruments)

#' @rdname MRCCCFit-accessors
#' @export
setMethod("posteriorDraws", "MRCCCFit",
          function(x, standardized = FALSE, ...) {
    .checkFlag(standardized, "standardized")
    d <- x@draws
    if (standardized) {
        sf <- std_factors(x)
        for (p in names(sf)) d[[p]] <- d[[p]] * sf[[p]]
    }
    d
})

#' @rdname MRCCCFit-accessors
#' @export
setMethod("tripletLabels", "MRCCCFit", function(x) x@labels)

#' @rdname MRCCCFit-accessors
#' @export
setMethod("firstStageEffects", "MRCCCFit", function(x) x@firstStage)

#' @rdname MRCCCFit-accessors
#' @export
setMethod("posteriorMeans", "MRCCCFit", function(x) x@posteriorMeans)

#' @rdname MRCCCFit-accessors
#' @export
setMethod("inputScale", "MRCCCFit", function(x) x@scale)

#' @rdname MRCCCFit-accessors
#' @export
setMethod("fitSettings", "MRCCCFit", function(x) x@settings)

.checkFlag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
    }
    invisible(x)
}

.tripletTitle <- function(lab) {
    if (all(is.na(lab[c("ligand", "receptor")]))) return(NULL)
    paste0(ifelse(is.na(lab[["ligand"]]), "ligand", lab[["ligand"]]), " -> ",
           ifelse(is.na(lab[["receptor"]]), "receptor", lab[["receptor"]]),
           if (!is.na(lab[["pathway"]])) paste0(" (", lab[["pathway"]], ")"))
}

#' @rdname MRCCCFit-class
#' @param object An [MRCCCFit-class] object.
#' @return `show()` prints a short summary of the fit and returns `object`
#'   invisibly.
#' @export
setMethod("show", "MRCCCFit", function(object) {
    s <- object@settings
    cat("MRCCCFit")
    ttl <- .tripletTitle(object@labels)
    if (!is.null(ttl)) cat(":", ttl)
    cat("\n")
    cat("  Chains: ", s$nChains, ", retained draws per chain: ", s$nKeep,
        " (nIter = ", s$nIter, ", burnIn = ", s$burnIn, ", thin = ",
        s$thin, ")\n", sep = "")
    cat("  Posterior inclusion probability: ", fmt(object@pip),
        "  (MCSE ", fmt(object@pipMCSE), ", ESS ", fmt(object@pipESS),
        ")\n", sep = "")
    fi <- object@instruments
    cat("  First-stage F: ligand = ", fmt(fi$F[1L]), " (",
        fi$nInstruments[1L], " SNPs), receptor = ", fmt(fi$F[2L]), " (",
        fi$nInstruments[2L], " SNPs)",
        if (any(fi$weak)) "  [weak: F < 10]", "\n", sep = "")
    e <- effectSizes(object)
    cat("  Standardised posterior means: beta_X = ", fmt(e$mean[1L]),
        ", beta_XZ = ", fmt(e$mean[2L]), ", beta_Z = ", fmt(e$mean[3L]),
        "\n", sep = "")
    rh <- object@diagnostics$rhat
    cat("  Max R-hat: ",
        if (all(is.na(rh))) "NA" else fmt(max(rh, na.rm = TRUE)),
        "\n", sep = "")
    invisible(object)
})

#' Summarise an MR-CCC fit
#'
#' Collects the posterior inclusion probability with its Monte Carlo error,
#' the standardised posterior means and credible intervals of the three
#' effect sizes, the instrument strength and the convergence diagnostics.
#'
#' @param object An [MRCCCFit-class] object returned by [mrccc()].
#' @param level Credible level of the reported intervals.
#' @param ... Not used.
#'
#' @return An [MRCCCSummary-class] object, which prints the summary.
#'   `as.data.frame()` extracts its table, with one row per effect
#'   (`beta_X`, `beta_XZ`, `beta_Z`) and columns `mean`, `lower` and
#'   `upper` on the standardised scale, and `rhat`.
#'
#' @seealso [mrccc()], [credibleIntervals()]
#'
#' @examples
#' sim <- simulateMRCCC(n = 200, seed = 1)
#' dat <- buildTriplet(sim$sender, sim$receiver, sim$genotypes,
#'                     ligand = "LIG1", receptor = "REC1",
#'                     pathway = sim$pathwayGenes,
#'                     covariates = c("sex", "age", "PC1", "PC2"))
#' fit <- mrccc(dat, nIter = 4000, burnIn = 1000, seed = 1)
#' s <- summary(fit)
#' s
#' as.data.frame(s)
#'
#' @export
setMethod("summary", "MRCCCFit", function(object, level = 0.95, ...) {
    ci <- credibleIntervals(object, level = level, standardized = TRUE)
    e <- effectSizes(object, standardized = TRUE)
    d <- object@diagnostics
    tab <- data.frame(
        parameter = e$parameter,
        mean = e$mean,
        lower = ci$lower,
        upper = ci$upper,
        rhat = c(d$rhat[match(c("beta_X", "beta_XZ"), d$parameter)],
                 NA_real_),
        stringsAsFactors = FALSE)
    s <- object@settings
    methods::new("MRCCCSummary",
        table = tab, diagnostics = d, instruments = object@instruments,
        info = list(title = .tripletTitle(object@labels),
                    nChains = s$nChains, nKeep = s$nKeep,
                    pip = object@pip, pipMCSE = object@pipMCSE,
                    pipESS = object@pipESS,
                    pipThreshold = s$pipThreshold, level = level))
})

#' @rdname MRCCCSummary-class
#' @param object An [MRCCCSummary-class] object.
#' @return `show()` prints the summary and returns `object` invisibly.
#' @export
setMethod("show", "MRCCCSummary", function(object) {
    info <- object@info
    tab <- object@table
    d <- object@diagnostics
    fi <- object@instruments
    cat("MR-CCC fit summary")
    if (!is.null(info$title)) cat(":", info$title)
    cat("\n")
    cat("  Chains: ", info$nChains, ", retained draws per chain: ",
        info$nKeep, "\n", sep = "")
    cat("  Posterior inclusion probability: ", fmt(info$pip),
        "  (MCSE ", fmt(info$pipMCSE), ", ESS ", fmt(info$pipESS),
        ")\n", sep = "")
    cat("  Discovery threshold: ", info$pipThreshold, "\n", sep = "")
    cat("  First-stage F: ligand = ", fmt(fi$F[1L]), ", receptor = ",
        fmt(fi$F[2L]),
        if (any(fi$weak)) "  (weak by the F > 10 heuristic)", "\n", sep = "")
    cat("\nStandardised effects (posterior mean and ",
        round(100 * info$level), "% credible interval):\n", sep = "")
    for (i in seq_len(nrow(tab))) {
        cat("  ", format(tab$parameter[i], width = 8), " ",
            fmt(tab$mean[i]), "  [", fmt(tab$lower[i]), ", ",
            fmt(tab$upper[i]), "]\n", sep = "")
    }
    cat("\nGelman-Rubin R-hat (NA unless nChains >= 2):\n")
    for (i in seq_len(nrow(d))) {
        cat("  ", format(d$parameter[i], width = 8), " ", fmt(d$rhat[i]),
            "\n", sep = "")
    }
    if (is.finite(info$pipMCSE) &&
        abs(info$pip - info$pipThreshold) < 2 * info$pipMCSE) {
        cat("\nNote: the PIP lies within two Monte Carlo standard errors of ",
            "the discovery threshold. The discovery status is not resolved ",
            "at this chain length. Increase nIter or nChains.\n", sep = "")
    }
    invisible(object)
})

#' @rdname MRCCCSummary-class
#' @param x An [MRCCCSummary-class] object.
#' @param row.names,optional Not used.
#' @return `as.data.frame()` returns the table of effects.
#' @export
as.data.frame.MRCCCSummary <- function(x, row.names = NULL,
                                       optional = FALSE, ...) {
    x@table
}

#' Plot an MR-CCC fit
#'
#' Plots the receptor-modulated ligand effect with its credible band. See
#' [plotEffectCurve()].
#'
#' @param x An [MRCCCFit-class] object.
#' @param y Not used.
#' @param ... Further arguments passed to [plotEffectCurve()].
#'
#' @return A `ggplot` object.
#'
#' @seealso [plotEffectCurve()]
#'
#' @examples
#' sim <- simulateMRCCC(n = 200, seed = 1)
#' dat <- buildTriplet(sim$sender, sim$receiver, sim$genotypes,
#'                     ligand = "LIG1", receptor = "REC1",
#'                     pathway = sim$pathwayGenes,
#'                     covariates = c("sex", "age", "PC1", "PC2"))
#' fit <- mrccc(dat, nIter = 4000, burnIn = 1000, seed = 1)
#' plot(fit)
#'
#' @export
setMethod("plot", signature(x = "MRCCCFit", y = "missing"),
          function(x, y, ...) plotEffectCurve(x, ...))
