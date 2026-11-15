#' @include MRCCCData.R selectInstruments.R
NULL

#' Build an MR-CCC triplet from SummarizedExperiment data
#'
#' Assembles the [MRCCCData-class] object for one ligand-receptor-pathway
#' triplet from donor-level (pseudo-bulk) expression of a sender and a
#' receiver cell type and from donor genotypes. Donors are matched across
#' the three objects by column name, covariates are taken from `colData`,
#' the pathway activity score is computed from the receiver's expression of
#' the pathway genes, and, unless supplied, the cis-eQTL instruments of the
#' ligand and the receptor are chosen with [selectInstruments()].
#'
#' @param sender A `SummarizedExperiment` (a `RangedSummarizedExperiment`
#'   when instruments are to be selected) of donor-level expression in the
#'   sender cell type, genes in rows and donors in columns.
#' @param receiver A `SummarizedExperiment` of donor-level expression in the
#'   receiver cell type, laid out as `sender`.
#' @param genotypes A `SummarizedExperiment` of allele dosages, SNPs in rows
#'   (with identifiers as row names) and donors in columns. A
#'   `RangedSummarizedExperiment` with SNP positions is required when
#'   instruments are to be selected. Dosages count copies of the alternative
#'   allele (0, 1 or 2, or an imputed value in between). Missing values are
#'   imputed by the SNP's mean dosage.
#' @param ligand Name of the ligand gene (a row name of `sender`).
#' @param receptor Name of the receptor gene (a row name of `receiver`).
#' @param pathway Either a character vector of pathway genes, scored in
#'   `receiver` according to `pathwayMethod`, or a numeric vector named by
#'   donor giving a precomputed pathway activity score (for example from
#'   AUCell, UCell or GSVA).
#' @param covariates Character vector naming donor covariates in
#'   `colData(receiver)` (or, failing that, `colData(sender)`), or a numeric
#'   matrix or data frame with donors as row names. Factors and character
#'   columns are expanded to indicator columns. At least one covariate is
#'   required. The same covariates adjust the instrument screen and the
#'   model.
#' @param ligandSNPs,receptorSNPs Optional instruments: character vectors of
#'   SNP identifiers (row names of `genotypes`), or the `DataFrame`s returned
#'   by [selectInstruments()]. When `NULL`, [selectInstruments()] is called.
#' @param assay Assay of `sender` and `receiver` holding the normalised
#'   expression: a name or an index.
#' @param genotypeAssay Assay of `genotypes` holding the dosages.
#' @param pathwayMethod How a character `pathway` is scored: `"mean"`, the
#'   mean expression of the pathway genes, or `"pc1"`, their first principal
#'   component (genes centred and scaled), oriented so that higher scores
#'   mean higher mean expression.
#' @param pathwayName,senderName,receiverName Optional labels used in
#'   printed output.
#' @param ... Further arguments passed to [selectInstruments()], such as
#'   `window`, `maxSNPs` or `r2Max`.
#'
#' @details
#' Only donors present in `sender`, `receiver` and `genotypes` are used, and
#' a message reports how many were dropped. The instrument screens use the
#' same donors and covariates as the model, so the analysis is one-sample.
#' The selection tables are kept in the `metadata` slot of the result under
#' `ligandInstruments` and `receptorInstruments`.
#'
#' @return An [MRCCCData-class] object.
#'
#' @seealso [selectInstruments()], [MRCCCData()], [mrccc()],
#'   [simulateMRCCC()]
#'
#' @examples
#' sim <- simulateMRCCC(n = 120, seed = 1)
#' dat <- buildTriplet(sim$sender, sim$receiver, sim$genotypes,
#'                     ligand = "LIG1", receptor = "REC1",
#'                     pathway = sim$pathwayGenes,
#'                     covariates = c("sex", "age", "PC1", "PC2"),
#'                     pathwayName = "Simulated pathway")
#' dat
#'
#' @export
buildTriplet <- function(sender, receiver, genotypes, ligand, receptor,
                         pathway, covariates, ligandSNPs = NULL,
                         receptorSNPs = NULL, assay = 1L, genotypeAssay = 1L,
                         pathwayMethod = c("mean", "pc1"),
                         pathwayName = NA_character_,
                         senderName = NA_character_,
                         receiverName = NA_character_, ...) {
    pathwayMethod <- match.arg(pathwayMethod)
    .checkSE(sender, "sender")
    .checkSE(receiver, "receiver")
    .checkSE(genotypes, "genotypes")
    for (a in list(list(ligand, "ligand"), list(receptor, "receptor"))) {
        if (!is.character(a[[1]]) || length(a[[1]]) != 1L || is.na(a[[1]])) {
            stop("'", a[[2]], "' must be a single gene name.", call. = FALSE)
        }
    }
    if (!ligand %in% rownames(sender)) {
        stop("Ligand '", ligand, "' is not a row name of 'sender'.",
             call. = FALSE)
    }
    if (!receptor %in% rownames(receiver)) {
        stop("Receptor '", receptor, "' is not a row name of 'receiver'.",
             call. = FALSE)
    }
    if (missing(covariates) || is.null(covariates)) {
        stop("'covariates' is required: name at least one donor covariate ",
             "in colData, such as sex or age.", call. = FALSE)
    }

    # ---- Donors shared by all three objects ----
    allIds <- list(sender = colnames(sender), receiver = colnames(receiver),
                   genotypes = colnames(genotypes))
    donorIds <- Reduce(intersect, allIds)
    V <- .tripletCovariates(covariates, sender, receiver)
    donorIds <- intersect(donorIds, rownames(V))
    if (is.numeric(pathway)) {
        if (is.null(names(pathway))) {
            stop("A numeric 'pathway' score must be named by donor.",
                 call. = FALSE)
        }
        donorIds <- intersect(donorIds, names(pathway))
    }
    if (length(donorIds) < 10L) {
        stop("Only ", length(donorIds), " donors are shared by 'sender', ",
             "'receiver', 'genotypes' and the covariates. At least 10 are ",
             "required. Check that donors are identified by the same column ",
             "names in every object.", call. = FALSE)
    }
    nDropped <- vapply(allIds, function(ids) sum(!ids %in% donorIds),
                       integer(1))
    if (any(nDropped > 0L)) {
        notUsed <- toString(sprintf("%d %s", nDropped, names(nDropped)))
        message("Using the ", length(donorIds), " donors shared by all ",
                "inputs (not used: ", notUsed, ").")
    }
    sender <- sender[, donorIds]
    receiver <- receiver[, donorIds]
    genotypes <- genotypes[, donorIds]
    V <- V[donorIds, , drop = FALSE]

    # ---- Exposures and outcome ----
    x <- .assayRow(sender, ligand, assay, donorIds, "sender")
    z <- .assayRow(receiver, receptor, assay, donorIds, "receiver")
    pw <- .pathwayScore(pathway, receiver, assay, donorIds, pathwayMethod)

    # ---- Instruments ----
    selL <- .resolveInstruments(ligandSNPs, sender, ligand, genotypes, V,
                                assay, genotypeAssay, "ligand", ...)
    selR <- .resolveInstruments(receptorSNPs, receiver, receptor, genotypes,
                                V, assay, genotypeAssay, "receptor", ...)
    G <- t(.imputeMeanDosage(.dosageMatrix(genotypes, genotypeAssay,
                                           selL$snp, donorIds)))
    H <- t(.imputeMeanDosage(.dosageMatrix(genotypes, genotypeAssay,
                                           selR$snp, donorIds)))

    MRCCCData(
        ligandExpression = x, receptorExpression = z,
        pathwayActivity = pw$score, ligandGenotypes = G,
        receptorGenotypes = H, covariates = V,
        labels = c(ligand = ligand, receptor = receptor,
                   pathway = as.character(pathwayName),
                   sender = as.character(senderName),
                   receiver = as.character(receiverName)),
        metadata = list(ligandInstruments = selL$table,
                        receptorInstruments = selR$table,
                        pathwayGenes = pw$genes,
                        pathwayMethod = pw$method))
}

# ---- Internal helpers ------------------------------------------------------

.tripletCovariates <- function(covariates, sender, receiver) {
    if (!is.character(covariates)) {
        V <- .covariateMatrix(covariates)
        if (is.null(rownames(V))) {
            stop("A covariate matrix or data frame must have donor ",
                 "identifiers as row names.", call. = FALSE)
        }
        return(V)
    }
    cdR <- SummarizedExperiment::colData(receiver)
    cdS <- SummarizedExperiment::colData(sender)
    cols <- lapply(covariates, function(nm) {
        if (nm %in% colnames(cdR)) {
            stats::setNames(cdR[[nm]], colnames(receiver))
        } else if (nm %in% colnames(cdS)) {
            stats::setNames(cdS[[nm]], colnames(sender))
        } else {
            stop("Covariate '", nm, "' is not a column of colData(receiver) ",
                 "or colData(sender).", call. = FALSE)
        }
    })
    ids <- Reduce(intersect, lapply(cols, names))
    df <- data.frame(lapply(cols, function(v) v[ids]),
                     row.names = ids, stringsAsFactors = FALSE)
    names(df) <- covariates
    .covariateMatrix(df)
}

.pathwayScore <- function(pathway, receiver, assay, donorIds, method) {
    if (is.numeric(pathway)) {
        score <- as.numeric(pathway[donorIds])
        names(score) <- donorIds
        if (any(!is.finite(score))) {
            stop("The pathway score contains missing values for some donors.",
                 call. = FALSE)
        }
        return(list(score = score, genes = character(0),
                    method = "precomputed"))
    }
    if (!is.character(pathway) || length(pathway) < 1L) {
        stop("'pathway' must be a character vector of genes or a numeric ",
             "score named by donor.", call. = FALSE)
    }
    genes <- intersect(unique(pathway), rownames(receiver))
    if (length(genes) == 0L) {
        stop("None of the pathway genes is a row name of 'receiver'.",
             call. = FALSE)
    }
    if (length(genes) < length(unique(pathway))) {
        message(length(unique(pathway)) - length(genes), " of ",
                length(unique(pathway)), " pathway genes are absent from ",
                "'receiver' and were ignored.")
    }
    a <- .getAssay(receiver, assay, "receiver")
    E <- as.matrix(a[genes, donorIds, drop = FALSE])
    storage.mode(E) <- "double"
    if (any(!is.finite(E))) {
        stop("Pathway gene expression contains missing values.",
             call. = FALSE)
    }
    meanScore <- colMeans(E)
    if (method == "mean") {
        score <- meanScore
    } else {
        varying <- apply(E, 1L, stats::var) > 0
        if (sum(varying) < 2L) {
            stop("pathwayMethod = \"pc1\" needs at least two pathway genes ",
                 "that vary across donors.", call. = FALSE)
        }
        pc <- stats::prcomp(t(E[varying, , drop = FALSE]), center = TRUE,
                            scale. = TRUE)$x[, 1L]
        r <- if (stats::sd(meanScore) > 0) stats::cor(pc, meanScore) else 0
        score <- if (r < 0) -pc else pc
    }
    score <- as.numeric(score)
    names(score) <- donorIds
    list(score = score, genes = genes, method = method)
}

.resolveInstruments <- function(snps, expression, gene, genotypes, V, assay,
                                genotypeAssay, role, ...) {
    if (is.null(snps)) {
        tab <- selectInstruments(expression, gene, genotypes, covariates = V,
                                 assay = assay, genotypeAssay = genotypeAssay,
                                 ...)
        if (nrow(tab) == 0L) {
            stop("No cis-eQTL instruments could be selected for the ", role,
                 " '", gene, "'.", call. = FALSE)
        }
        return(list(snp = as.character(tab$snp), table = tab))
    }
    if (methods::is(snps, "DataFrame") || is.data.frame(snps)) {
        if (!"snp" %in% colnames(snps)) {
            stop("A table of ", role, " instruments must have a 'snp' column.",
                 call. = FALSE)
        }
        tab <- snps
        snps <- as.character(snps$snp)
    } else {
        tab <- S4Vectors::DataFrame(snp = as.character(snps))
    }
    if (!is.character(snps) || length(snps) < 1L) {
        stop("The ", role, " instruments must be SNP identifiers.",
             call. = FALSE)
    }
    miss <- setdiff(snps, rownames(genotypes))
    if (length(miss)) {
        stop(length(miss), " ", role, " instrument(s) are not row names of ",
             "'genotypes', for example '", miss[1L], "'.", call. = FALSE)
    }
    list(snp = snps, table = tab)
}
