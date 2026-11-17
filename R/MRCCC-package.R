#' MRCCC: Bayesian Mendelian randomization for causal cell-cell communication
#'
#' MRCCC implements MR-CCC, a Bayesian Mendelian randomization framework for
#' inferring causal ligand-receptor signalling between cell types from
#' population-scale single-cell expression and genotype data. For each
#' ligand-receptor-pathway triplet, cis-eQTL genotypes of the same donors
#' instrument ligand expression in the sender cell type and receptor
#' expression in the receiver cell type. A working model with a
#' ligand-by-receptor interaction is then fitted to downstream pathway
#' activity by a blocked Gibbs sampler. A spike-and-slab prior on the causal block yields
#' a posterior inclusion probability that serves as the communication score.
#'
#' @section Workflow:
#' \describe{
#'   \item{[buildTriplet()]}{Assemble one triplet from donor-level
#'     `SummarizedExperiment` objects for the sender and receiver cell types
#'     and a genotype `RangedSummarizedExperiment`, selecting cis-eQTL
#'     instruments with [selectInstruments()]. [MRCCCData()] builds the same
#'     object from aligned vectors and matrices.}
#'   \item{[mrccc()]}{Fit the model. Returns an [MRCCCFit-class] object.}
#'   \item{[MRCCCFit-accessors]}{[pip()], [effectSizes()], [convergence()],
#'     [instrumentStrength()], [firstStageEffects()], [posteriorMeans()],
#'     [posteriorDraws()], [inputScale()] and [fitSettings()].}
#'   \item{Posterior summaries}{[credibleIntervals()],
#'     [signProbabilities()], [signReversal()] and [plotEffectCurve()].}
#'   \item{[mrcccScreen()], [screenTable()]}{Fit many triplets, in
#'     parallel through BiocParallel, and summarise them.}
#'   \item{[bayesFDR()], [bfdrSelect()]}{Bayesian false discovery rate
#'     control across many triplets.}
#'   \item{[onek1kExample()]}{A real-data example from the OneK1K cohort.}
#'   \item{[simulateMRCCC()]}{Simulated data in the layout of real data.}
#' }
#'
#' A worked analysis, on simulated and real data, is in
#' `vignette("MRCCC", package = "MRCCC")`.
#'
#' @references
#' Sarkar B, Ni Y (2026). MR-CCC: Bayesian Mendelian Randomization for
#' Causal Cell-Cell Communication. Bioinformatics (in press).
#' arXiv:2604.23917.
#'
#' @useDynLib MRCCC, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom methods new is slot setClass setGeneric setMethod
#'   setValidity show
#' @importFrom ggplot2 .data
#' @keywords internal
"_PACKAGE"
