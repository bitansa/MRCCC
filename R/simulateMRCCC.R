#' Simulate expression and genotype data for one triplet
#'
#' Generates donor-level pseudo-bulk expression for a sender and a receiver
#' cell type and donor genotypes, stored as `RangedSummarizedExperiment`
#' objects, under the MR-CCC data-generating process with an unmeasured
#' confounder. The output has the layout of real data and can be passed to
#' [selectInstruments()] and [buildTriplet()].
#'
#' @param n Number of donors.
#' @param nInstruments Number of causal cis-eQTL SNPs for the ligand and for
#'   the receptor.
#' @param betaX True ligand effect on the pathway at the mean receptor
#'   level.
#' @param betaXZ True ligand-by-receptor interaction effect.
#' @param betaZ True receptor main effect on the pathway.
#' @param instrumentEffect Effect on expression of one standard deviation
#'   of each causal SNP's dosage (instrument strength).
#' @param confounding Loading of an unmeasured confounder on ligand
#'   expression, receptor expression and pathway activity. Zero removes
#'   confounding.
#' @param nPathwayGenes Number of genes in the downstream pathway.
#' @param seed Optional integer. If supplied, the data are simulated under
#'   this seed and the caller's random number state is restored afterwards.
#'
#' @details
#' Donors have covariates `sex`, `age` and two genotype principal components
#' `PC1` and `PC2`, stored in the `colData` of both expression objects.
#' Genotypes are allele dosages (0, 1 or 2 copies of the alternative
#' allele). Around each of the ligand gene `LIG1` (chromosome 1, plus strand)
#' and the receptor gene `REC1` (chromosome 2, minus strand) the genotype
#' object contains `nInstruments` causal SNPs, one strongly linked proxy for
#' each causal SNP, six SNPs with no effect, one rare SNP and, outside the
#' 200 kb window, five distant SNPs. With \eqn{U} the confounder and
#' \eqn{\tilde G, \tilde H} the standardised causal dosages,
#' \deqn{X = \tilde G \pi + V\alpha + cU + \varepsilon_X, \quad
#'   Z = \tilde H \pi + V\alpha + cU + \varepsilon_Z,}
#' \deqn{Y = \beta_X X + \beta_Z Z + \beta_{XZ} XZ + V\alpha + cU +
#'   \varepsilon_Y,}
#' where every element of \eqn{\pi} is `instrumentEffect`, every element of
#' \eqn{\alpha} is 0.3, \eqn{c} is `confounding` and the errors are standard
#' normal. The sender expresses `LIG1` at \eqn{5 + X}. The receiver
#' expresses `REC1` at \eqn{4 + Z} and each pathway gene at \eqn{6 + Y} plus
#' independent noise with standard deviation 0.5, so the mean expression of
#' the pathway genes tracks \eqn{Y}. Both objects also contain unrelated
#' genes. Setting `betaX = betaXZ = 0` gives a null triplet.
#'
#' @return A list with elements
#' \describe{
#'   \item{`sender`}{`RangedSummarizedExperiment` of sender expression
#'     (assay `logcounts`).}
#'   \item{`receiver`}{`RangedSummarizedExperiment` of receiver expression
#'     (assay `logcounts`).}
#'   \item{`genotypes`}{`RangedSummarizedExperiment` of allele dosages
#'     (assay `dosage`), SNPs in rows.}
#'   \item{`pathwayGenes`}{Names of the pathway genes.}
#'   \item{`truth`}{List with the generating values `betaX`, `betaXZ`,
#'     `betaZ`, the causal SNP identifiers `ligandCausal` and
#'     `receptorCausal`, and `tauStd`, the true sign-reversal threshold
#'     \eqn{-\beta_X / \beta_{XZ}} in standard deviations of receptor
#'     expression from its mean (`NA` when `betaXZ` is zero).}
#' }
#'
#' @seealso [buildTriplet()], [selectInstruments()], [mrccc()]
#'
#' @examples
#' sim <- simulateMRCCC(n = 120, seed = 1)
#' sim$sender
#' sim$genotypes
#' sim$truth$betaXZ
#'
#' @export
simulateMRCCC <- function(n = 300L, nInstruments = 4L, betaX = 0.3,
                          betaXZ = 0.3, betaZ = 0.5, instrumentEffect = 0.8,
                          confounding = 0.7, nPathwayGenes = 8L,
                          seed = NULL) {
    for (nm in c("n", "nInstruments", "nPathwayGenes")) {
        v <- get(nm)
        check_scalar(v, nm)
        if (v < 1 || v != round(v)) {
            stop("'", nm, "' must be a positive integer.", call. = FALSE)
        }
    }
    if (n < 20) stop("'n' must be at least 20.", call. = FALSE)
    for (nm in c("betaX", "betaXZ", "betaZ", "instrumentEffect",
                 "confounding")) {
        check_scalar(get(nm), nm)
    }
    if (!is.null(seed)) {
        check_scalar(seed, "seed")
        if (seed != round(seed)) {
            stop("'seed' must be a whole number.", call. = FALSE)
        }
        withr::local_seed(as.integer(seed))
    }
    n <- as.integer(n)
    k <- as.integer(nInstruments)
    donorIds <- sprintf("donor%03d", seq_len(n))

    # ---- Donor covariates ----
    cd <- S4Vectors::DataFrame(
        sex = factor(sample(c("female", "male"), n, replace = TRUE)),
        age = round(stats::runif(n, 20, 80)),
        PC1 = stats::rnorm(n), PC2 = stats::rnorm(n),
        row.names = donorIds)
    Vs <- scale(cbind(as.numeric(cd$sex == "male"), cd$age, cd$PC1, cd$PC2))
    U <- stats::rnorm(n)

    # ---- Genotypes around the ligand and receptor genes ----
    # Every object uses the same set of chromosome levels, so that ranges
    # combine and overlap without sequence-level warnings.
    chroms <- c("chr1", "chr2", "chr5", "chr7", "chr9")
    lig <- .simulateLocus(n, k, chrom = "chr1", tss = 1e6, prefix = "rsL",
                          levels = chroms)
    rec <- .simulateLocus(n, k, chrom = "chr2", tss = 5e6, prefix = "rsR",
                          levels = chroms)

    X <- drop(scale(lig$causal) %*% rep(instrumentEffect, k) +
                  Vs %*% rep(0.3, 4L) + confounding * U + stats::rnorm(n))
    Z <- drop(scale(rec$causal) %*% rep(instrumentEffect, k) +
                  Vs %*% rep(0.3, 4L) + confounding * U + stats::rnorm(n))
    Y <- drop(betaX * X + betaZ * Z + betaXZ * X * Z +
                  Vs %*% rep(0.3, 4L) + confounding * U + stats::rnorm(n))

    # ---- Expression objects ----
    pwGenes <- sprintf("PW%02d", seq_len(nPathwayGenes))
    otherGenes <- function(prefix, m) {
        matrix(stats::rnorm(m * n, 5, 1), m, n,
               dimnames = list(paste0(prefix, seq_len(m)), NULL))
    }
    senderMat <- rbind(LIG1 = 5 + X, otherGenes("SGENE", 4L))
    pwMat <- t(vapply(pwGenes, function(g) 6 + Y + stats::rnorm(n, 0, 0.5),
                      numeric(n)))
    receiverMat <- rbind(REC1 = 4 + Z, pwMat, otherGenes("RGENE", 3L))
    colnames(senderMat) <- colnames(receiverMat) <- donorIds

    geneRanges <- function(genes, chrom, start, strand) {
        GenomicRanges::GRanges(
            seqnames = factor(chrom, levels = chroms),
            ranges = IRanges::IRanges(start = start, width = 20000L),
            strand = strand, gene_id = genes)
    }
    senderRanges <- geneRanges(
        rownames(senderMat), c("chr1", rep("chr5", 4L)),
        c(1e6, 1e6 + 1e6 * seq_len(4L)), c("+", rep("+", 4L)))
    # REC1 lies on the minus strand, so its transcription start site is the
    # end of its range, at 5 Mb.
    receiverRanges <- geneRanges(
        rownames(receiverMat),
        c("chr2", rep("chr7", nPathwayGenes), rep("chr9", 3L)),
        c(5e6 - 20000L + 1L, 1e6 * seq_len(nPathwayGenes),
          1e6 * seq_len(3L)),
        c("-", rep("+", nPathwayGenes + 3L)))
    names(senderRanges) <- rownames(senderMat)
    names(receiverRanges) <- rownames(receiverMat)

    sender <- SummarizedExperiment::SummarizedExperiment(
        assays = list(logcounts = senderMat), rowRanges = senderRanges,
        colData = cd)
    receiver <- SummarizedExperiment::SummarizedExperiment(
        assays = list(logcounts = receiverMat), rowRanges = receiverRanges,
        colData = cd)

    dosage <- rbind(lig$all, rec$all)
    colnames(dosage) <- donorIds
    snpRanges <- c(lig$ranges, rec$ranges)
    genotypes <- SummarizedExperiment::SummarizedExperiment(
        assays = list(dosage = dosage), rowRanges = snpRanges)

    list(sender = sender, receiver = receiver, genotypes = genotypes,
         pathwayGenes = pwGenes,
         truth = list(betaX = betaX, betaXZ = betaXZ, betaZ = betaZ,
                      ligandCausal = rownames(lig$causalT),
                      receptorCausal = rownames(rec$causalT),
                      tauStd = if (betaXZ != 0) {
                          (-betaX / betaXZ - mean(Z)) / stats::sd(Z)
                      } else NA_real_))
}

# Simulate the SNPs of one locus. Returns the donors-by-SNPs matrix of causal
# dosages, the SNPs-by-donors matrix of all dosages, and their GRanges.
.simulateLocus <- function(n, k, chrom, tss, prefix, levels) {
    draw <- function(m, lo = 0.15, hi = 0.5) {
        maf <- stats::runif(m, lo, hi)
        vapply(maf, function(p) stats::rbinom(n, 2L, p), numeric(n))
    }
    redraw <- function(M, frac = 0.05) {
        # Re-draw a small fraction of donors so the proxy is in strong but
        # imperfect linkage disequilibrium with its causal SNP.
        apply(M, 2L, function(g) {
            idx <- sample.int(n, max(1L, round(frac * n)))
            g[idx] <- stats::rbinom(length(idx), 2L, mean(g) / 2)
            g
        })
    }
    causal <- draw(k)
    proxy <- redraw(causal)
    null <- draw(6L)
    rare <- matrix(stats::rbinom(n, 2L, 0.01), ncol = 1L)
    distant <- draw(5L)

    near <- sort(sample(seq(-150000L, 150000L, by = 97L), k * 2L + 7L))
    pos <- c(tss + near, tss + 1e6 + 10000L * seq_len(5L))
    all <- cbind(causal, proxy, null, rare, distant)
    ids <- sprintf("%s%04d", prefix, seq_len(ncol(all)))
    # Shuffle identifiers against positions so that SNP order is not
    # informative of the causal ones.
    perm <- c(sample.int(k * 2L + 7L), k * 2L + 7L + seq_len(5L))
    colnames(all) <- ids[order(perm)]
    allT <- t(all)
    rng <- GenomicRanges::GRanges(
        seqnames = factor(rep(chrom, length(perm)), levels = levels),
        ranges = IRanges::IRanges(start = pos[perm], width = 1L))
    names(rng) <- colnames(all)
    ord <- order(pos[perm])
    list(causal = all[, seq_len(k), drop = FALSE],
         causalT = allT[seq_len(k), , drop = FALSE],
         all = allT[ord, , drop = FALSE],
         ranges = rng[ord])
}
