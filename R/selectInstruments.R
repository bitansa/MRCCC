#' Select cis-eQTL instruments for one gene
#'
#' Selects the cis-eQTL SNPs that serve as instrumental variables for the
#' expression of one gene (a ligand in the sender cell type or a receptor in
#' the receiver cell type), following the procedure of the MR-CCC paper:
#' SNPs in a window around the gene's promoter are screened in-sample by a
#' covariate-adjusted regression of the gene's donor-level expression on
#' allele dosage, ranked by p-value, and retained greedily subject to a
#' linkage-disequilibrium (LD) bound until at most `maxSNPs` are kept.
#'
#' @param expression A `RangedSummarizedExperiment` of donor-level
#'   (pseudo-bulk) expression for one cell type, genes in rows and donors in
#'   columns. Its `rowRanges` must give the genomic location of `gene`.
#' @param gene Name of the gene (a row name of `expression`).
#' @param genotypes A `RangedSummarizedExperiment` of genotypes, SNPs in rows
#'   and donors in columns, with SNP identifiers as row names and SNP
#'   positions as `rowRanges`. The assay holds allele dosages: the number of
#'   copies of the alternative allele (0, 1 or 2), or an imputed dosage in
#'   between. Missing values are allowed.
#' @param covariates Donor covariates to adjust the screen for: either a
#'   character vector naming columns of `colData(expression)`, or a numeric
#'   matrix or data frame with donors as row names. `NULL` adjusts for an
#'   intercept only. Use the same covariates as in the model.
#' @param assay Assay of `expression` to use: a name or an index.
#' @param genotypeAssay Assay of `genotypes` to use: a name or an index.
#' @param window Half-width in base pairs of the window around the
#'   transcription start site in which candidate SNPs are sought. The default
#'   200 kb is the paper's setting.
#' @param maxSNPs Maximum number of instruments retained.
#' @param r2Max Maximum squared dosage correlation (LD \eqn{r^2}) between
#'   two retained instruments. A value of 1 disables clumping.
#' @param minCarrierFreq SNPs carried (dosage of at least 0.5) by a fraction
#'   of donors no greater than this are removed before screening.
#' @param minCandidates Minimum number of SNPs with a valid screening
#'   p-value required. With fewer, no instruments are returned.
#'
#' @details
#' Only donors present in both `expression` and `genotypes` (and with
#' covariates, when supplied as a matrix) are used. The steps are:
#' 1. *Window.* Candidate SNPs lie within `window` base pairs upstream or
#'    downstream of the transcription start site of `gene`, as given by
#'    [IRanges::promoters()] applied to its range.
#' 2. *Carrier filter.* SNPs carried by at most `minCarrierFreq` of the
#'    donors with a called genotype are removed. Remaining missing dosages
#'    are imputed by the SNP's mean dosage.
#' 3. *Screen.* For each candidate the gene's expression is regressed on the
#'    dosage and the covariates (with an intercept), and the two-sided
#'    p-value of the dosage coefficient is computed. SNPs whose dosage is
#'    constant or collinear with the covariates are skipped.
#' 4. *LD clumping.* SNPs are visited in order of increasing p-value. A SNP
#'    is kept if its squared dosage correlation with every SNP already kept
#'    is at most `r2Max`, until `maxSNPs` are kept. The strongest SNP is
#'    always kept.
#'
#' No p-value threshold is applied: the strength of the selected set is
#' assessed in the model through the first-stage F statistic reported by
#' [instrumentStrength()]. Because the screen uses the same donors as the
#' model, MR-CCC is a one-sample analysis.
#'
#' @return A `DataFrame` with one row per selected SNP, in the order of
#'   selection, and columns `snp` (identifier), `estimate` (screening
#'   coefficient per copy of the alternative allele) and `pvalue`. Its
#'   `metadata()` records the gene, the window, the number of donors, the
#'   number of SNPs in the window and the number screened. When fewer than
#'   `minCandidates` SNPs can be screened, a zero-row `DataFrame` is
#'   returned with a message.
#'
#' @seealso [buildTriplet()], which calls this function for the ligand and
#'   the receptor.
#'
#' @references
#' Sarkar B, Ni Y (2026). MR-CCC: Bayesian Mendelian Randomization for
#' Causal Cell-Cell Communication. arXiv:2604.23917.
#'
#' @examples
#' sim <- simulateMRCCC(n = 120, seed = 1)
#' sel <- selectInstruments(sim$sender, "LIG1", sim$genotypes,
#'                          covariates = c("sex", "age", "PC1", "PC2"))
#' sel
#' S4Vectors::metadata(sel)
#'
#' @export
selectInstruments <- function(expression, gene, genotypes, covariates = NULL,
                              assay = 1L, genotypeAssay = 1L,
                              window = 2e5, maxSNPs = 10L, r2Max = 0.8,
                              minCarrierFreq = 0.05, minCandidates = 5L) {
    .checkRangedSE(expression, "expression")
    .checkRangedSE(genotypes, "genotypes")
    if (!is.character(gene) || length(gene) != 1L || is.na(gene)) {
        stop("'gene' must be a single gene name.", call. = FALSE)
    }
    if (!gene %in% rownames(expression)) {
        stop("Gene '", gene, "' is not a row name of 'expression'.",
             call. = FALSE)
    }
    if (is.null(rownames(genotypes))) {
        stop("'genotypes' must have SNP identifiers as row names.",
             call. = FALSE)
    }
    check_scalar(window, "window")
    check_scalar(maxSNPs, "maxSNPs")
    check_scalar(r2Max, "r2Max")
    check_scalar(minCarrierFreq, "minCarrierFreq")
    check_scalar(minCandidates, "minCandidates")
    if (window <= 0) stop("'window' must be positive.", call. = FALSE)
    if (maxSNPs < 1 || maxSNPs != round(maxSNPs)) {
        stop("'maxSNPs' must be a positive integer.", call. = FALSE)
    }
    if (r2Max <= 0 || r2Max > 1) {
        stop("'r2Max' must lie in (0, 1].", call. = FALSE)
    }
    if (minCarrierFreq < 0 || minCarrierFreq >= 1) {
        stop("'minCarrierFreq' must lie in [0, 1).", call. = FALSE)
    }

    # ---- Donors and covariates ----
    donorIds <- intersect(colnames(expression), colnames(genotypes))
    V <- .resolveCovariates(covariates, expression)
    if (!is.null(V)) donorIds <- intersect(donorIds, rownames(V))
    if (length(donorIds) < 10L) {
        stop("Only ", length(donorIds), " donors are shared by 'expression'",
             ", 'genotypes' and the covariates. At least 10 are required. ",
             "Check that donors are identified by the same column names.",
             call. = FALSE)
    }
    y <- .assayRow(expression, gene, assay, donorIds, "expression")
    if (!is.null(V)) V <- V[donorIds, , drop = FALSE]

    # ---- Window around the transcription start site ----
    geneRange <- .featureRange(expression, gene)
    win <- GenomicRanges::promoters(geneRange, upstream = window,
                                    downstream = window)
    snpRanges <- SummarizedExperiment::rowRanges(genotypes)
    # overlapsAny() warns whenever the two objects carry different sets of
    # chromosomes, which is the normal situation here (the window lies on one
    # chromosome, the genotypes span many). Those warnings are muffled. When
    # the naming styles differ entirely, no overlap is found and the message
    # below says what to check.
    inWindow <- which(withCallingHandlers(
        IRanges::overlapsAny(snpRanges, win, ignore.strand = TRUE),
        warning = function(w) {
            if (grepl("sequence levels", conditionMessage(w))) {
                invokeRestart("muffleWarning")
            }
        }))

    meta <- list(gene = gene, window = window, nDonors = length(donorIds),
                 nInWindow = length(inWindow), nScreened = 0L)
    empty <- S4Vectors::DataFrame(snp = character(0), estimate = numeric(0),
                                  pvalue = numeric(0))

    if (length(inWindow) == 0L) {
        message("No SNPs of 'genotypes' lie within ", window, " bp of the ",
                "transcription start site of '", gene, "'. If SNPs are ",
                "expected there, check that both objects name chromosomes ",
                "the same way (for example '1' versus 'chr1').")
        S4Vectors::metadata(empty) <- meta
        return(empty)
    }

    # ---- Dosages, carrier filter and imputation ----
    D <- .dosageMatrix(genotypes, genotypeAssay, rownames(genotypes)[inWindow],
                       donorIds)
    carrier <- rowMeans(D >= 0.5, na.rm = TRUE)
    keepCarrier <- is.finite(carrier) & carrier > minCarrierFreq
    D <- D[keepCarrier, , drop = FALSE]
    D <- .imputeMeanDosage(D)

    # ---- Covariate-adjusted screen ----
    screen <- .screenSNPs(y, D, V)
    ok <- which(is.finite(screen$pvalue))
    meta$nScreened <- length(ok)
    if (length(ok) < minCandidates) {
        message("Only ", length(ok), " SNP(s) near '", gene, "' could be ",
                "screened (minCandidates = ", minCandidates, "). No ",
                "instruments returned.")
        S4Vectors::metadata(empty) <- meta
        return(empty)
    }

    # ---- LD clumping in p-value order ----
    ranked <- ok[order(screen$pvalue[ok])]
    kept <- .clumpByLD(ranked, D, keepMax = as.integer(maxSNPs),
                       r2Max = r2Max)

    out <- S4Vectors::DataFrame(snp = rownames(D)[kept],
                                estimate = screen$estimate[kept],
                                pvalue = screen$pvalue[kept])
    S4Vectors::metadata(out) <- meta
    out
}

# ---- Internal helpers ------------------------------------------------------

.checkSE <- function(x, arg) {
    if (!methods::is(x, "SummarizedExperiment")) {
        stop("'", arg, "' must be a SummarizedExperiment with features in ",
             "rows and donors in columns.", call. = FALSE)
    }
    if (is.null(colnames(x)) || anyDuplicated(colnames(x)) > 0L) {
        stop("'", arg, "' must have unique donor identifiers as column ",
             "names.", call. = FALSE)
    }
    if (is.null(rownames(x))) {
        stop("'", arg, "' must have feature names as row names.",
             call. = FALSE)
    }
    invisible(x)
}

.checkRangedSE <- function(x, arg) {
    .checkSE(x, arg)
    if (!methods::is(x, "RangedSummarizedExperiment")) {
        stop("'", arg, "' must be a RangedSummarizedExperiment whose ",
             "rowRanges give genomic positions. These are needed to find ",
             "cis-SNPs.", call. = FALSE)
    }
    invisible(x)
}

.featureRange <- function(se, feature) {
    rr <- SummarizedExperiment::rowRanges(se)
    if (!methods::is(rr, "GRanges")) {
        stop("rowRanges of the expression object must be a GRanges with one ",
             "range per gene.", call. = FALSE)
    }
    rr[match(feature, rownames(se))]
}

.assayRow <- function(se, feature, assay, donorIds, arg) {
    a <- .getAssay(se, assay, arg)
    v <- as.numeric(as.matrix(a[feature, donorIds, drop = FALSE]))
    names(v) <- donorIds
    if (any(!is.finite(v))) {
        stop("Expression of '", feature, "' in '", arg, "' contains missing ",
             "or infinite values for the selected donors.", call. = FALSE)
    }
    v
}

.getAssay <- function(se, assay, arg) {
    nms <- SummarizedExperiment::assayNames(se)
    nA <- length(SummarizedExperiment::assays(se))
    if (is.character(assay)) {
        if (length(assay) != 1L || !assay %in% nms) {
            stop("'", assay, "' is not an assay of '", arg, "'. Available: ",
                 toString(nms), ".", call. = FALSE)
        }
    } else if (!is.numeric(assay) || length(assay) != 1L ||
               assay != round(assay) || assay < 1 || assay > nA) {
        stop("The assay of '", arg, "' must be a name or an index between 1 ",
             "and ", nA, ".", call. = FALSE)
    }
    SummarizedExperiment::assay(se, assay, withDimnames = TRUE)
}

.dosageMatrix <- function(genotypes, genotypeAssay, snps, donorIds) {
    a <- .getAssay(genotypes, genotypeAssay, "genotypes")
    D <- as.matrix(a[snps, donorIds, drop = FALSE])
    if (!is.numeric(D)) {
        stop("The genotype assay must hold numeric allele dosages (0, 1 or 2 ",
             "copies of the alternative allele).", call. = FALSE)
    }
    storage.mode(D) <- "double"
    if (any(D[!is.na(D)] < 0 | D[!is.na(D)] > 2)) {
        stop("Genotype dosages must lie in [0, 2] (copies of the ",
             "alternative allele).", call. = FALSE)
    }
    D
}

.imputeMeanDosage <- function(D) {
    if (!anyNA(D)) return(D)
    mu <- rowMeans(D, na.rm = TRUE)
    idx <- which(is.na(D), arr.ind = TRUE)
    D[idx] <- mu[idx[, 1L]]
    D
}

.resolveCovariates <- function(covariates, se) {
    if (is.null(covariates)) return(NULL)
    if (is.character(covariates)) {
        cd <- SummarizedExperiment::colData(se)
        miss <- setdiff(covariates, colnames(cd))
        if (length(miss)) {
            stop("Covariate(s) ", toString(miss),
                 " not found in colData.", call. = FALSE)
        }
        df <- methods::as(cd[, covariates, drop = FALSE], "data.frame")
        rownames(df) <- colnames(se)
        return(.covariateMatrix(df))
    }
    V <- .covariateMatrix(covariates)
    if (is.null(rownames(V))) {
        stop("A covariate matrix or data frame must have donor identifiers ",
             "as row names.", call. = FALSE)
    }
    V
}

# Regression of y on each row of D (one SNP at a time) with an intercept and
# the covariates V, by the Frisch-Waugh-Lovell theorem: both sides are
# residualised on [1, V] once, after which each SNP's coefficient, standard
# error and t statistic follow in closed form. Identical to lm(y ~ g + V).
.screenSNPs <- function(y, D, V) {
    n <- length(y)
    W <- if (is.null(V)) matrix(1, n, 1L) else cbind(1, V)
    qrW <- qr(W)
    yr <- qr.resid(qrW, y)
    Gt <- t(D)
    Gr <- qr.resid(qrW, Gt)
    ss <- colSums(Gr^2)
    ssRaw <- colSums(sweep(Gt, 2L, colMeans(Gt))^2)
    valid <- ssRaw > 0 & ss > 1e-8 * pmax(ssRaw, 1e-300)
    df <- n - qrW$rank - 1L
    est <- rep(NA_real_, ncol(Gt))
    pv <- rep(NA_real_, ncol(Gt))
    if (df > 0L && any(valid)) {
        b <- colSums(Gr[, valid, drop = FALSE] * yr) / ss[valid]
        rss <- pmax(sum(yr^2) - b^2 * ss[valid], 0)
        se <- sqrt(rss / df / ss[valid])
        tt <- b / se
        est[valid] <- b
        pv[valid] <- 2 * stats::pt(-abs(tt), df)
    }
    list(estimate = est, pvalue = pv)
}

# Greedy LD clumping over rows of D, visited in the order of `cand`.
.clumpByLD <- function(cand, D, keepMax, r2Max) {
    if (r2Max >= 1) return(cand[seq_len(min(keepMax, length(cand)))])
    kept <- integer(0)
    for (s in cand) {
        if (length(kept) >= keepMax) break
        if (length(kept) == 0L) {
            kept <- s
            next
        }
        # Every candidate varies across donors (the screen skips constant
        # SNPs), so the correlations are defined.
        r <- stats::cor(D[s, ], t(D[kept, , drop = FALSE]))
        if (any(is.na(r))) next
        if (max(r^2) <= r2Max) kept <- c(kept, s)
    }
    kept
}
