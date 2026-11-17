## Builds inst/extdata/onek1k_nk_monocytes.rds, the real-data example of the
## MRCCC package, from the public OneK1K data deposited on Zenodo.
##
## Source data
##   Zenodo record 10.5281/zenodo.19675075 (files B_T_NK_monocytes.rda and
##   donor.rda): donor-level single-cell counts of the OneK1K cohort
##   (Yazar et al., Science 2022) and donor genotypes, as prepared for the
##   MR-CCC paper. CellPhoneDB v5 and MSigDB Reactome define the
##   ligand-receptor pairs and the pathway gene set.
##
## What the example contains
##   The NK cells (sender) to monocytes (receiver) axis of the paper, restricted
##   to the ligand-receptor pairs that CellPhoneDB classifies as "Signaling by
##   Interleukin". Donors pass the quality control of the paper. Expression is
##   pseudo-bulk counts divided by a library-size factor, exactly as in the
##   paper. Genotypes are allele dosages of the SNPs within 200 kb of the
##   transcription start site of each ligand and receptor gene, after the
##   paper's carrier-frequency filter. The receiver colData holds the donor
##   covariates and two activity scores of the Reactome interleukin gene set:
##   its first principal component and the mean expression of its genes.
##
## How to run
##   Run from the root of the analysis repository
##   (https://github.com/bitansa/MR-CCC), with the MRCCC package installed and
##   the environment variables below pointing at the downloaded data:
##     Sys.setenv(MRCCC_DATA_DIR  = "/path/to/zenodo/files")
##     Sys.setenv(CELLPHONEDB_DIR = "/path/to/cellphonedb/csvs")
##     source("/path/to/MRCCC/inst/scripts/make-onek1k-example.R")
##   The output file is written to the working directory and is then copied
##   to inst/extdata of the package.

suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(GenomicRanges)
    library(dplyr)
    library(stringr)
})

## ---- Ligand-receptor database and Reactome gene sets -----------------------
## Produces lr_matrix and key_cols, as in the paper.
source("Github Codes/build_lr_database.R")

## ---- OneK1K data ------------------------------------------------------------
dataDir <- Sys.getenv("MRCCC_DATA_DIR",
                      unset = file.path("~", "data", "mrccc"))
load(file.path(dataDir, "B_T_NK_monocytes.rda"))
load(file.path(dataDir, "donor.rda"))

## Donor-level counts, genes by donors. Donors are aligned by position across
## the expression lists, the genotype matrix and the donor table.
countsNK <- sapply(NK_Cells, colSums)
getCounts <- function(x) {
    if (!is.null(dim(x)) && length(dim(x)) == 2L) return(colSums(x))
    x
}
countsMono <- sapply(monocytes, getCounts)
validGenes <- rownames(sapply(B_Cells, colSums))
rm(B_Cells, CD4_Cells, CD8_Cells, NK_Cells, monocytes)
invisible(gc())
stopifnot(ncol(countsNK) == nrow(donor), ncol(countsMono) == nrow(donor),
          ncol(genotype) == nrow(donor))

## Donor identifiers.
donorIds <- if (!is.null(colnames(genotype))) {
    colnames(genotype)
} else {
    sprintf("donor%04d", seq_len(nrow(donor)))
}

## ---- Donor quality control and normalisation (as in the paper) -------------
libSize <- function(m) colSums(m) / stats::median(colSums(m))
ls1 <- libSize(countsNK)
ls2 <- libSize(countsMono)
keep <- (ls1 <= stats::median(ls1) + 3 * stats::mad(ls1)) & (ls1 > 0.25) &
    (ls2 <= stats::median(ls2) + 3 * stats::mad(ls2)) & (ls2 > 0.25)
donor <- donor[keep, ]
genotype <- genotype[, keep]
donorIds <- donorIds[keep]
countsNK <- countsNK[, keep]
countsMono <- countsMono[, keep]
normNK <- sweep(countsNK, 2L, libSize(countsNK), "/")
normMono <- sweep(countsMono, 2L, libSize(countsMono), "/")
colnames(normNK) <- colnames(normMono) <- donorIds
cat("Donors after quality control:", length(donorIds), "\n")

## ---- Interleukin triplets on the NK cells to monocytes axis ----------------
pairs <- lr_matrix %>%
    filter(ligand_ensembl %in% validGenes,
           receptor_ensembl %in% validGenes,
           sub(".*by\\s+", "", pathway_name) == "Interleukin") %>%
    distinct(ligand_symbol, ligand_ensembl, receptor_symbol,
             receptor_ensembl, pathway_name)
pairs <- pairs[pairs$ligand_ensembl %in% rownames(normNK) &
               pairs$receptor_ensembl %in% rownames(normMono) &
               pairs$ligand_ensembl %in% gene.ref$gene_id &
               pairs$receptor_ensembl %in% gene.ref$gene_id, ]
cat("Interleukin ligand-receptor pairs:", nrow(pairs), "\n")

## ---- Pathway activity (Reactome gene sets matching "Interleukin") ----------
pwGenes <- key_cols %>%
    filter(ensembl_gene %in% validGenes,
           str_detect(gs_name, regex("Interleukin", ignore_case = TRUE))) %>%
    pull(ensembl_gene) %>% unique()
pwGenes <- intersect(pwGenes, rownames(normMono))
pwExpr <- normMono[pwGenes, , drop = FALSE]
pwExpr <- pwExpr[apply(pwExpr, 1L, stats::var) > 0, , drop = FALSE]
pwMean <- colMeans(pwExpr)
pwPC1 <- stats::prcomp(t(pwExpr), center = TRUE, scale. = TRUE)$x[, 1L]
if (stats::cor(pwPC1, pwMean) < 0) pwPC1 <- -pwPC1
cat("Interleukin pathway genes scored:", nrow(pwExpr), "\n")

## ---- Donor covariates (as in the paper) -------------------------------------
ageBin <- cut(donor$age, breaks = c(-Inf, 30, 40, 50, 60, 70, 80, Inf),
              right = FALSE, labels = FALSE)
cd <- S4Vectors::DataFrame(
    sex = factor(donor$sex, levels = c("female", "male")),
    ageBin = as.integer(ageBin),
    PC1 = donor$PC1, PC2 = donor$PC2, PC3 = donor$PC3,
    interleukinPC1 = as.numeric(pwPC1),
    interleukinMean = as.numeric(pwMean),
    row.names = donorIds)
stopifnot(!anyNA(as.data.frame(cd)))

## ---- Expression objects -----------------------------------------------------
geneRange <- function(ens, sym) {
    idx <- match(ens, gene.ref$gene_id)
    stopifnot(!anyNA(idx))
    gr <- gene.ref[idx]
    S4Vectors::mcols(gr) <- S4Vectors::DataFrame(gene_id = ens)
    names(gr) <- sym
    gr
}
## Expression rows are named by gene symbol, matching the rowRanges names.
geneMatrix <- function(m, ens, sym) {
    out <- m[ens, , drop = FALSE]
    rownames(out) <- sym
    out
}
ligands <- unique(pairs[, c("ligand_symbol", "ligand_ensembl")])
receptors <- unique(pairs[, c("receptor_symbol", "receptor_ensembl")])
sender <- SummarizedExperiment(
    assays = list(normcounts = geneMatrix(normNK, ligands$ligand_ensembl,
                                          ligands$ligand_symbol)),
    rowRanges = geneRange(ligands$ligand_ensembl, ligands$ligand_symbol),
    colData = cd[, c("sex", "ageBin", "PC1", "PC2", "PC3")])
receiver <- SummarizedExperiment(
    assays = list(normcounts = geneMatrix(normMono,
                                          receptors$receptor_ensembl,
                                          receptors$receptor_symbol)),
    rowRanges = geneRange(receptors$receptor_ensembl,
                          receptors$receptor_symbol),
    colData = cd)

## ---- Genotypes of the cis windows ------------------------------------------
## SNPs are first restricted to the cis windows, then coded as dosages.
windows <- promoters(c(rowRanges(sender), rowRanges(receiver)),
                     upstream = 200000, downstream = 200000)
inWindow <- which(overlapsAny(vcf.ref, windows, ignore.strand = TRUE))
snpIds <- rownames(genotype)
if (is.null(snpIds)) snpIds <- names(vcf.ref)
stopifnot(!is.null(snpIds), anyDuplicated(snpIds[inWindow]) == 0L)
gt <- genotype[inWindow, , drop = FALSE]
dosage <- matrix(NA_integer_, nrow(gt), ncol(gt),
                 dimnames = list(snpIds[inWindow], donorIds))
dosage[gt == "0/0"] <- 0L
dosage[gt == "0/1"] <- 1L
dosage[gt == "1/1"] <- 2L
carrier <- rowMeans(dosage > 0, na.rm = TRUE)
keepSnp <- is.finite(carrier) & carrier > 0.05
snpRanges <- vcf.ref[inWindow[keepSnp]]
names(snpRanges) <- rownames(dosage)[keepSnp]
S4Vectors::mcols(snpRanges) <- NULL
genotypes <- SummarizedExperiment(
    assays = list(dosage = dosage[keepSnp, , drop = FALSE]),
    rowRanges = snpRanges)
cat("SNPs in the cis windows:", nrow(genotypes), "\n")

## ---- Triplets with instruments for both genes -------------------------------
covs <- c("sex", "ageBin", "PC1", "PC2", "PC3")
hasInstruments <- function(se, gene) {
    sel <- suppressMessages(MRCCC::selectInstruments(se, gene, genotypes,
                                                     covariates = covs))
    nrow(sel) > 0L
}
okLig <- vapply(ligands$ligand_symbol, hasInstruments, logical(1),
                se = sender)
okRec <- vapply(receptors$receptor_symbol, hasInstruments, logical(1),
                se = receiver)
triplets <- data.frame(
    ligand = pairs$ligand_symbol,
    receptor = pairs$receptor_symbol,
    pathway = "Signaling by Interleukin",
    stringsAsFactors = FALSE)
triplets <- triplets[okLig[triplets$ligand] & okRec[triplets$receptor], ]
rownames(triplets) <- NULL
print(triplets)

## Keep only the genes and SNPs that the remaining triplets use.
sender <- sender[unique(triplets$ligand), ]
receiver <- receiver[unique(triplets$receptor), ]
windows <- promoters(c(rowRanges(sender), rowRanges(receiver)),
                     upstream = 200000, downstream = 200000)
genotypes <- genotypes[overlapsAny(rowRanges(genotypes), windows,
                                   ignore.strand = TRUE), ]

example <- list(sender = sender, receiver = receiver,
                genotypes = genotypes, triplets = triplets)
saveRDS(example, "onek1k_nk_monocytes.rds", compress = "xz")
cat("Written onek1k_nk_monocytes.rds:",
    round(file.size("onek1k_nk_monocytes.rds") / 1e6, 2), "MB\n")
