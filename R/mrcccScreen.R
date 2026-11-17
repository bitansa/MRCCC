#' @include mrccc.R MRCCCFit-methods.R
NULL

#' Fit MR-CCC to many triplets
#'
#' Fits [mrccc()] to every triplet in a list and returns the fits, ready for
#' Bayesian false discovery rate control with [bayesFDR()], [bfdrSelect()]
#' or [screenTable()]. Triplets are fitted in parallel through
#' BiocParallel, and the chains within each triplet are run one after
#' another.
#'
#' @param triplets A list of [MRCCCData-class] objects, usually named. Each
#'   is typically made with [buildTriplet()].
#' @param ... Further arguments passed to [mrccc()], such as `nIter`,
#'   `burnIn`, `nChains` or `hyper`.
#' @param seed Optional integer. If supplied, one seed per triplet is
#'   derived from it, so the screen is reproducible whatever the back end.
#' @param BPPARAM A [BiocParallel::BiocParallelParam-class] object that
#'   controls how the triplets are distributed. The default
#'   [BiocParallel::SerialParam()] fits them one after another.
#'
#' @return A list of [MRCCCFit-class] objects with the names of `triplets`.
#'   When `triplets` is unnamed, the names are built from the ligand and
#'   receptor labels.
#'
#' @seealso [screenTable()], [bayesFDR()], [bfdrSelect()], [mrccc()]
#'
#' @examples
#' sims <- lapply(c(1, 2), function(s) simulateMRCCC(n = 200, seed = s))
#' triplets <- lapply(sims, function(sim) {
#'     buildTriplet(sim$sender, sim$receiver, sim$genotypes,
#'                  ligand = "LIG1", receptor = "REC1",
#'                  pathway = sim$pathwayGenes,
#'                  covariates = c("sex", "age", "PC1", "PC2"))
#' })
#' names(triplets) <- c("replicate1", "replicate2")
#' fits <- mrcccScreen(triplets, nIter = 4000, burnIn = 1000, seed = 1)
#' screenTable(fits)
#'
#' @export
mrcccScreen <- function(triplets, ..., seed = NULL,
                        BPPARAM = BiocParallel::SerialParam()) {
    if (methods::is(triplets, "MRCCCData")) triplets <- list(triplets)
    if (!is.list(triplets) || length(triplets) == 0L ||
        !all(vapply(triplets, methods::is, logical(1),
                    class2 = "MRCCCData"))) {
        stop("'triplets' must be a non-empty list of MRCCCData objects.",
             call. = FALSE)
    }
    if (!methods::is(BPPARAM, "BiocParallelParam")) {
        stop("'BPPARAM' must be a BiocParallelParam object, for example ",
             "BiocParallel::SerialParam().", call. = FALSE)
    }
    if (!is.null(seed)) {
        check_scalar(seed, "seed")
        if (seed != round(seed) || abs(seed) > .Machine$integer.max) {
            stop("'seed' must be a whole number no larger in absolute value ",
                 "than .Machine$integer.max.", call. = FALSE)
        }
    }
    seeds <- .chainSeeds(length(triplets), seed)
    args <- list(...)
    fitOne <- function(k) {
        do.call(mrccc, c(list(triplets[[k]]), args, list(seed = seeds[[k]])))
    }
    fits <- BiocParallel::bplapply(seq_along(triplets), fitOne,
                                   BPPARAM = BPPARAM)
    names(fits) <- if (!is.null(names(triplets))) {
        names(triplets)
    } else {
        vapply(triplets, function(t) {
            lab <- tripletLabels(t)
            paste(lab[["ligand"]], lab[["receptor"]], sep = "-")
        }, character(1))
    }
    fits
}

#' Summarise a screen of many triplets
#'
#' Collects the posterior inclusion probability, its Monte Carlo error and
#' the largest Gelman-Rubin statistic of every fit, and applies Bayesian
#' false discovery rate control across the triplets.
#'
#' @param fits A list of [MRCCCFit-class] objects, as returned by
#'   [mrcccScreen()].
#' @param alpha Target Bayesian false discovery rate for the selected set.
#'
#' @return A data frame with one row per fit, ordered by decreasing PIP, and
#'   columns `triplet` (the name of the fit), `ligand`, `receptor`,
#'   `pathway`, `pip`, `pipMCSE`, `maxRhat` (the largest Gelman-Rubin
#'   statistic, `NA` for a single chain), `bFDR` (the Bayesian false
#'   discovery rate of the PIP-ranked set ending at that triplet) and
#'   `selected` (the largest set with Bayesian FDR at most `alpha`).
#'
#' @seealso [mrcccScreen()], [bayesFDR()], [bfdrSelect()]
#'
#' @examples
#' sim <- simulateMRCCC(n = 200, seed = 1)
#' dat <- buildTriplet(sim$sender, sim$receiver, sim$genotypes,
#'                     ligand = "LIG1", receptor = "REC1",
#'                     pathway = sim$pathwayGenes,
#'                     covariates = c("sex", "age", "PC1", "PC2"))
#' fits <- mrcccScreen(list(example = dat), nIter = 4000, burnIn = 1000,
#'                     seed = 1)
#' screenTable(fits, alpha = 0.1)
#'
#' @export
screenTable <- function(fits, alpha = 0.10) {
    if (methods::is(fits, "MRCCCFit")) fits <- list(fits)
    if (!is.list(fits) || length(fits) == 0L ||
        !all(vapply(fits, methods::is, logical(1), class2 = "MRCCCFit"))) {
        stop("'fits' must be a non-empty list of MRCCCFit objects.",
             call. = FALSE)
    }
    nm <- names(fits)
    if (is.null(nm)) nm <- paste0("triplet", seq_along(fits))
    lab <- vapply(fits, function(f) tripletLabels(f), character(5))
    maxRhat <- vapply(fits, function(f) {
        rh <- convergence(f)$rhat
        if (all(is.na(rh))) NA_real_ else max(rh, na.rm = TRUE)
    }, numeric(1))
    tab <- data.frame(
        triplet = nm,
        ligand = lab["ligand", ],
        receptor = lab["receptor", ],
        pathway = lab["pathway", ],
        pip = vapply(fits, pip, numeric(1)),
        pipMCSE = vapply(fits, function(f) pipMCSE(f)[["mcse"]], numeric(1)),
        maxRhat = maxRhat,
        bFDR = bayesFDR(fits),
        selected = bfdrSelect(fits, alpha = alpha),
        row.names = NULL, stringsAsFactors = FALSE)
    tab <- tab[order(tab$pip, decreasing = TRUE), , drop = FALSE]
    rownames(tab) <- NULL
    tab
}
