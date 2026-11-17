# S4 classes of MRCCC: the input object for one ligand-receptor-pathway
# triplet (MRCCCData), the fitted model (MRCCCFit) and its summary
# (MRCCCSummary).

.labelNames <- c("ligand", "receptor", "pathway", "sender", "receiver")

#' The MRCCCData class: one ligand-receptor-pathway triplet
#'
#' An `MRCCCData` object holds everything MR-CCC needs for one
#' ligand-receptor-pathway triplet, aligned by donor: ligand expression in
#' the sender cell type, receptor expression and pathway activity in the
#' receiver cell type, the cis-eQTL genotypes that instrument the ligand and
#' the receptor, and the donor covariates. Objects are normally created from
#' `SummarizedExperiment` data with [buildTriplet()], or from aligned
#' vectors and matrices with the constructor [MRCCCData()]. The validity
#' method checks that every component refers to the same donors in the same
#' order, that genotypes are allele dosages in \[0, 2\], and that no value is
#' missing.
#'
#' @slot ligandExpression Named numeric vector, one value per donor: ligand
#'   expression in the sender cell type.
#' @slot receptorExpression Named numeric vector: receptor expression in the
#'   receiver cell type.
#' @slot pathwayActivity Named numeric vector: pathway activity score in the
#'   receiver cell type.
#' @slot ligandGenotypes Numeric matrix, donors by SNPs: allele dosages
#'   (0, 1 or 2 copies of the alternative allele, or an imputed value in
#'   between) of the cis-eQTL SNPs that instrument the ligand.
#' @slot receptorGenotypes Numeric matrix, donors by SNPs: allele dosages of
#'   the cis-eQTL SNPs that instrument the receptor.
#' @slot covariates Numeric matrix, donors by covariates: donor-level
#'   covariates (for example sex, age and genotype principal components),
#'   with factors already expanded to indicator columns.
#' @slot labels Named character vector with elements `ligand`, `receptor`,
#'   `pathway`, `sender` and `receiver` (any may be `NA`), used in printed
#'   output.
#' @slot metadata List of additional information, inherited from
#'   [S4Vectors::Annotated-class] and read with [S4Vectors::metadata()].
#'   [buildTriplet()] records there the instrument-selection tables
#'   (`ligandInstruments`, `receptorInstruments`), the pathway genes
#'   (`pathwayGenes`) and how the pathway was scored (`pathwayMethod`).
#'
#' @seealso [MRCCCData()], [buildTriplet()], [mrccc()]
#' @name MRCCCData-class
#' @rdname MRCCCData-class
#' @importClassesFrom S4Vectors Annotated
#' @exportClass MRCCCData
setClass("MRCCCData",
    contains = "Annotated",
    slots = c(
        ligandExpression   = "numeric",
        receptorExpression = "numeric",
        pathwayActivity    = "numeric",
        ligandGenotypes    = "matrix",
        receptorGenotypes  = "matrix",
        covariates         = "matrix",
        labels             = "character"
    )
)

.validMRCCCData <- function(object) {
    msg <- character(0)
    x <- object@ligandExpression
    z <- object@receptorExpression
    y <- object@pathwayActivity
    G <- object@ligandGenotypes
    H <- object@receptorGenotypes
    V <- object@covariates
    n <- length(x)

    lens <- c(ligandExpression = length(x), receptorExpression = length(z),
              pathwayActivity = length(y), ligandGenotypes = nrow(G),
              receptorGenotypes = nrow(H), covariates = nrow(V))
    if (any(lens != n)) {
        return(paste0("all components must describe the same donors. ",
                      "Lengths and row counts are: ",
                      paste(names(lens), lens, sep = " = ", collapse = ", ")))
    }
    if (n < 10L) {
        msg <- c(msg, sprintf(
            "at least 10 donors are required, but %d were found", n))
    }

    donors <- names(x)
    if (is.null(donors) || anyNA(donors) || anyDuplicated(donors) > 0L) {
        msg <- c(msg, paste0("donor names (names of 'ligandExpression') ",
                             "must be present and unique"))
    } else {
        same <- c(receptorExpression = identical(names(z), donors),
                  pathwayActivity    = identical(names(y), donors),
                  ligandGenotypes    = identical(rownames(G), donors),
                  receptorGenotypes  = identical(rownames(H), donors),
                  covariates         = identical(rownames(V), donors))
        if (!all(same)) {
            msg <- c(msg, paste0("donors are not aligned with the ligand ",
                                 "expression in: ",
                                 paste(names(same)[!same], collapse = ", ")))
        }
    }

    for (nm in c("ligandExpression", "receptorExpression", "pathwayActivity")) {
        v <- methods::slot(object, nm)
        if (any(!is.finite(v))) {
            msg <- c(msg, sprintf("'%s' contains NA, NaN or infinite values",
                                  nm))
        } else if (n >= 2L && stats::var(v) <= 0) {
            msg <- c(msg, sprintf("'%s' is constant across donors", nm))
        }
    }

    for (nm in c("ligandGenotypes", "receptorGenotypes")) {
        M <- methods::slot(object, nm)
        if (!is.numeric(M)) {
            msg <- c(msg, sprintf("'%s' must be numeric", nm))
            next
        }
        if (ncol(M) < 1L) {
            msg <- c(msg, sprintf("'%s' must contain at least one SNP", nm))
        }
        if (is.null(colnames(M)) || anyDuplicated(colnames(M)) > 0L) {
            msg <- c(msg, sprintf(paste0("'%s' must have unique SNP ",
                                         "identifiers as column names"), nm))
        }
        if (any(!is.finite(M))) {
            msg <- c(msg, sprintf("'%s' contains missing or infinite dosages",
                                  nm))
        } else if (any(M < 0 | M > 2)) {
            msg <- c(msg, sprintf(paste0("'%s' must hold allele dosages in ",
                                         "[0, 2] (copies of the alternative ",
                                         "allele)"), nm))
        }
    }

    if (!is.numeric(V)) {
        msg <- c(msg, "'covariates' must be numeric")
    } else {
        if (ncol(V) < 1L) {
            msg <- c(msg, "'covariates' must contain at least one column")
        }
        if (is.null(colnames(V))) {
            msg <- c(msg, "'covariates' must have column names")
        }
        if (any(!is.finite(V))) {
            msg <- c(msg, "'covariates' contains NA, NaN or infinite values")
        }
    }

    lab <- object@labels
    if (!identical(names(lab), .labelNames)) {
        msg <- c(msg, paste0("'labels' must be a character vector named ",
                             paste(.labelNames, collapse = ", ")))
    }

    if (length(msg)) msg else TRUE
}

setValidity("MRCCCData", .validMRCCCData)

#' The MRCCCFit class: a fitted MR-CCC model
#'
#' An `MRCCCFit` object is returned by [mrccc()]. It stores the pooled
#' posterior draws of the causal effects and the inclusion indicator, the
#' posterior inclusion probability (PIP) with its Monte Carlo error, the
#' convergence diagnostics, the first-stage instrument strength, and the
#' settings of the run. Use the accessors rather than the slots: [pip()],
#' [pipMCSE()], [effectSizes()], [credibleIntervals()],
#' [signProbabilities()], [signReversal()], [convergence()],
#' [instrumentStrength()], [firstStageEffects()], [posteriorMeans()],
#' [posteriorDraws()], [inputScale()], [fitSettings()] and
#' [plotEffectCurve()].
#'
#' @slot pip Posterior inclusion probability, pooled over chains.
#' @slot pipMCSE Monte Carlo standard error of `pip`.
#' @slot pipESS Effective sample size of the pooled indicator chain.
#' @slot draws Data frame of retained draws, one row per draw, with columns
#'   `chain`, `beta_X`, `beta_XZ`, `beta_Z`, `gamma`, `mu` (intercept) and
#'   `loglik` (joint log-likelihood of the three equations), the
#'   coefficients on the centred input scale. Read with [posteriorDraws()].
#' @slot scale List with the standard deviations (`sdX`, `sdZ`, `sdY`) and
#'   means (`meanX`, `meanZ`, `meanY`) of ligand expression, receptor
#'   expression and pathway activity, and the observed receptor range in
#'   standard deviations from its mean (`zRangeStd`). Read with
#'   [inputScale()].
#' @slot diagnostics Data frame with columns `parameter` (`beta_X`,
#'   `beta_XZ`, `gamma`, `loglik`) and `rhat`, the Gelman-Rubin statistic
#'   (`NA` for a single chain). Read with [convergence()].
#' @slot instruments Data frame with one row per exposure and columns
#'   `exposure` (`ligand`, `receptor`), `nInstruments`, `F` (first-stage
#'   partial F statistic after the covariates) and `weak` (`F < 10`). Read
#'   with [instrumentStrength()].
#' @slot firstStage Data frame with one row per instrument and columns
#'   `exposure`, `snp`, `mean`, `lower` and `upper`: the posterior mean and
#'   95% credible interval of the instrument's effect on expression per copy
#'   of the alternative allele. Read with [firstStageEffects()].
#' @slot posteriorMeans List of the posterior means of the remaining
#'   parameters: `covariateEffects` (data frame with columns `covariate`,
#'   `ligand`, `receptor` and `pathway`), `residualVariances` (named vector
#'   for the ligand, receptor and pathway equations),
#'   `inclusionProbabilityPrior` (\eqn{\rho}) and `intercept`
#'   (\eqn{\mu}). Read with [posteriorMeans()].
#' @slot settings List of the run settings: `nIter`, `burnIn`, `thin`,
#'   `nChains`, `nKeep` (retained draws per chain), `pipThreshold`, `hyper`
#'   (the hyperparameters), `g` (the g-prior scale used), `initScale` and
#'   `priorScale` (the fixed second-stage prior scale `d_X`, `d_XZ`,
#'   `d_Z`). Read with [fitSettings()].
#' @slot labels Labels copied from the [MRCCCData-class] object.
#'
#' @seealso [mrccc()]
#' @name MRCCCFit-class
#' @rdname MRCCCFit-class
#' @exportClass MRCCCFit
setClass("MRCCCFit",
    slots = c(
        pip         = "numeric",
        pipMCSE     = "numeric",
        pipESS      = "numeric",
        draws       = "data.frame",
        scale       = "list",
        diagnostics = "data.frame",
        instruments = "data.frame",
        firstStage  = "data.frame",
        posteriorMeans = "list",
        settings    = "list",
        labels      = "character"
    )
)

.validMRCCCFit <- function(object) {
    msg <- character(0)
    p <- object@pip
    if (length(p) != 1L || !is.finite(p) || p < 0 || p > 1) {
        msg <- c(msg, "'pip' must be a single number in [0, 1]")
    }
    need <- c("chain", "beta_X", "beta_XZ", "beta_Z", "gamma", "mu", "loglik")
    miss <- setdiff(need, names(object@draws))
    if (length(miss)) {
        msg <- c(msg, paste0("'draws' lacks column(s): ",
                             paste(miss, collapse = ", ")))
    }
    if (length(msg)) msg else TRUE
}

setValidity("MRCCCFit", .validMRCCCFit)

#' The MRCCCSummary class: a printed summary of an MR-CCC fit
#'
#' Returned by [summary()][summary,MRCCCFit-method]. Printing it shows the
#' posterior inclusion probability with its Monte Carlo error, the
#' standardised effects with credible intervals, the instrument strength and
#' the convergence diagnostics. `as.data.frame()` returns the table of
#' effects.
#'
#' @slot table Data frame with one row per effect and columns `parameter`
#'   (`beta_X`, `beta_XZ`, `beta_Z`), `mean`, `lower` and `upper` (posterior
#'   mean and credible interval on the standardised scale) and `rhat`. The
#'   `rhat` of `beta_Z` is `NA` because R-hat is computed for `beta_X`,
#'   `beta_XZ`, `gamma` and the log-likelihood only.
#' @slot diagnostics Data frame of Gelman-Rubin statistics, as returned by
#'   [convergence()].
#' @slot instruments Data frame of first-stage instrument strength, as
#'   returned by [instrumentStrength()].
#' @slot info List with `title` (the triplet labels), `nChains`, `nKeep`,
#'   `pip`, `pipMCSE`, `pipESS`, `pipThreshold` and `level` (the credible
#'   level).
#' @param ... Not used.
#'
#' @seealso [summary,MRCCCFit-method]
#' @name MRCCCSummary-class
#' @rdname MRCCCSummary-class
#' @exportClass MRCCCSummary
setClass("MRCCCSummary",
    slots = c(
        table       = "data.frame",
        diagnostics = "data.frame",
        instruments = "data.frame",
        info        = "list"
    )
)

