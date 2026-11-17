#' Real-data example: OneK1K NK cells to monocytes
#'
#' Loads a small subset of the OneK1K cohort (Yazar et al., 2022) that ships
#' with the package, for demonstrating MR-CCC on real data. It covers the NK
#' cells (sender) to monocytes (receiver) axis analysed in the MR-CCC paper,
#' restricted to the ligand-receptor pairs that CellPhoneDB classifies as
#' "Signaling by Interleukin".
#'
#' @details
#' **Provenance.** The data were derived from the files deposited on Zenodo
#' for the MR-CCC paper (\doi{10.5281/zenodo.19675075}) by the script
#' `inst/scripts/make-onek1k-example.R`, which can be found with
#' `system.file("scripts", "make-onek1k-example.R", package = "MRCCC")`.
#' Donors pass the quality control of the paper. Expression is donor-level
#' (pseudo-bulk) counts divided by a library-size factor, as in the paper.
#' Genotypes are allele dosages of the SNPs within 200 kb of the
#' transcription start site of each ligand and receptor gene that are
#' carried by more than 5% of donors. The receiver `colData` holds the donor
#' covariates used in the paper (`sex`, `ageBin` for age in seven groups
#' coded 1 to 7, namely under 30, five 10-year bins and 80 or over, and the genotype principal components `PC1`, `PC2` and
#' `PC3`) and two activity scores of the Reactome interleukin gene set:
#' `interleukinPC1`, the first principal component of the genes (centred,
#' scaled, and oriented so that higher scores mean higher mean expression),
#' and `interleukinMean`, their mean expression.
#'
#' @return A list with elements
#' \describe{
#'   \item{`sender`}{`RangedSummarizedExperiment` of NK-cell expression of
#'     the ligand genes (assay `normcounts`), donors in columns.}
#'   \item{`receiver`}{`RangedSummarizedExperiment` of monocyte expression
#'     of the receptor genes (assay `normcounts`), with the covariates and
#'     pathway scores in `colData`.}
#'   \item{`genotypes`}{`RangedSummarizedExperiment` of allele dosages
#'     (assay `dosage`), SNPs in rows.}
#'   \item{`triplets`}{Data frame of the ligand-receptor-pathway triplets
#'     with instruments for both genes, with columns `ligand`, `receptor`
#'     and `pathway`.}
#' }
#'
#' @references
#' Yazar S, Alquicira-Hernandez J, Wing K, et al. (2022). Single-cell eQTL
#' mapping identifies cell type-specific genetic control of autoimmune
#' disease. Science 376(6589):eabf3041.
#'
#' Sarkar B, Ni Y (2026). MR-CCC: Bayesian Mendelian Randomization for
#' Causal Cell-Cell Communication. Bioinformatics (in press).
#' arXiv:2604.23917.
#'
#' @seealso [buildTriplet()], [mrcccScreen()]
#'
#' @examples
#' ex <- onek1kExample()
#' ex$receiver
#' ex$triplets
#'
#' @export
onek1kExample <- function() {
    path <- system.file("extdata", "onek1k_nk_monocytes.rds",
                        package = "MRCCC")
    if (!nzchar(path)) {
        stop("The example data file is not installed with MRCCC.",
             call. = FALSE)
    }
    readRDS(path)
}
