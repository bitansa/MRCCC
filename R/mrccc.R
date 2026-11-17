#' @include AllClasses.R AllGenerics.R MRCCCData.R
NULL

#' Fit the MR-CCC model to one ligand-receptor-pathway triplet
#'
#' Fits the Bayesian Mendelian randomization model for causal cell-cell
#' communication to one triplet stored in an [MRCCCData-class] object. The
#' cis-eQTL genotypes instrument ligand expression in the sender cell type
#' and receptor expression in the receiver cell type. The model asks whether
#' the genetically predicted ligand level, modulated by the genetically
#' predicted receptor level, shifts pathway activity in the receiver. Several
#' Gibbs chains are run and pooled, and the posterior inclusion probability
#' (PIP) is returned with its Monte Carlo standard error and convergence
#' diagnostics.
#'
#' @param x An [MRCCCData-class] object, usually from [buildTriplet()].
#' @param nIter Total number of Gibbs iterations per chain.
#' @param burnIn Number of initial iterations discarded per chain. Must be
#'   smaller than `nIter`.
#' @param thin Thinning interval: every `thin`-th post-burn-in iteration is
#'   retained.
#' @param nChains Number of independent chains, pooled for all posterior
#'   summaries. With two or more chains the Gelman-Rubin statistic is
#'   reported and the chains start from dispersed values. Pooling `C` chains
#'   of `nIter` iterations gives roughly the Monte Carlo precision of one
#'   chain of `C * (nIter - burnIn)` retained draws. Each chain repeats the
#'   burn-in, which is the only extra cost.
#' @param initScale Standard deviation of the normal distribution from which
#'   \eqn{\beta_X} and \eqn{\beta_{XZ}} are initialised, on the centred
#'   input scale, when `nChains > 1`. It must then be positive, so that the
#'   chains start from different values. The default 1 gives standard normal
#'   starting values. With a single chain both effects start at zero.
#' @param pipThreshold Threshold in \eqn{(0, 1)} that defines a discovery
#'   (`PIP > pipThreshold`). The default of 0.5 corresponds to the median
#'   probability model. Used for the borderline warning and by
#'   [summary()][summary,MRCCCFit-method].
#' @param hyper Prior hyperparameters, as returned by
#'   [mrcccHyperparameters()].
#' @param seed Optional integer. If supplied, the fit is reproducible: one
#'   seed per chain is derived from it, each chain runs under its own seed
#'   wherever it is evaluated, and the caller's random number state is
#'   restored afterwards. With the default `NULL`, the per-chain seeds are
#'   drawn from the current random number stream, so calling [set.seed()]
#'   beforehand also makes the fit reproducible.
#' @param verbose Logical. If `TRUE`, each chain prints a progress line every
#'   1000 iterations.
#' @param BPPARAM A [BiocParallel::BiocParallelParam-class] object that
#'   controls how the chains are run. The default
#'   [BiocParallel::SerialParam()] runs them one after another. For example,
#'   `BiocParallel::MulticoreParam(4)` runs four chains at once. Because each
#'   chain has its own seed, the result does not depend on the back end.
#' @param ... Not used.
#'
#' @details
#' **Model.** For donor \eqn{i}, with \eqn{X^* = E[X \mid G, V]} and
#' \eqn{Z^* = E[Z \mid H, V]} the genetically predicted ligand and receptor
#' levels, the outcome model is
#' \deqn{Y = \mu + \beta_X X^* + \beta_Z Z^* + \beta_{XZ} X^* Z^*
#'   + \alpha_Y^\top V + \tilde\varepsilon_Y,}
#' where \eqn{X} is ligand expression, \eqn{Z} receptor expression, \eqn{Y}
#' pathway activity, \eqn{G} and \eqn{H} the instrument dosages and \eqn{V}
#' the covariates. The causal block \eqn{(\beta_X, \beta_{XZ})} has a
#' spike-and-slab prior with inclusion indicator
#' \eqn{\gamma \in \{0, 1\}}. Under \eqn{\gamma = 0} the prior variance of
#' the block is multiplied by `nu1` (\eqn{10^{-4}} by default), so the
#' coefficients are concentrated tightly at zero. The PIP \eqn{P(\gamma = 1 \mid \mathrm{data})} is the
#' communication score. The receptor-modulated ligand effect is
#' \eqn{\beta_X + \beta_{XZ} z}, which changes sign at
#' \eqn{\tau = -\beta_X / \beta_{XZ}}. See [signReversal()].
#'
#' **Centring and scale.** Neither first stage carries an intercept, so the
#' expression, pathway, genotype and covariate data are all centred
#' internally. Effects are reported on the centred input scale and on a
#' standardised scale (per standard deviation of \eqn{X} and \eqn{Z}, in
#' standard deviations of \eqn{Y}). See [effectSizes()].
#'
#' **Priors.** The first-stage blocks and the covariate block carry Zellner
#' g-priors on their designs. The second-stage blocks are scaled by a fixed
#' diagonal computed once from the least-squares first stage, so every
#' Gibbs update is an exact conjugate step. Unless overridden in `hyper`,
#' every scale \eqn{g} is \eqn{\min(n, 100)}. See [mrcccHyperparameters()].
#'
#' **Instrument strength.** The first-stage partial F statistic of each
#' instrument set, after the covariates, is reported by
#' [instrumentStrength()]. A message is issued when either F is below 10.
#' MR-CCC does not require F > 10, but PIPs are less decisive with weak
#' instruments and longer chains may be needed.
#'
#' **Chains, Monte Carlo error and convergence.** Odd-numbered chains start
#' with \eqn{\gamma = 1} and even-numbered chains with \eqn{\gamma = 0}.
#' The Monte Carlo standard error (MCSE) of the PIP is estimated by
#' non-overlapping batch means within each chain and combined across chains.
#' It is `NA` when the indicator never changes value in any chain, so that
#' the PIP is exactly 0 or 1. A warning is issued when
#' \eqn{|\mathrm{PIP} - t| < 2\,\mathrm{MCSE}} for threshold \eqn{t},
#' because the discovery status is then not resolved at the current chain
#' length, and when any Gelman-Rubin \eqn{\hat R} exceeds
#' 1.01. \eqn{\hat R} is reported for \eqn{\beta_X}, \eqn{\beta_{XZ}},
#' \eqn{\gamma} and the joint log-likelihood. It requires at least two
#' chains. The default chain length is a starting point. Reported analyses
#' should use chains long enough that neither warning is issued.
#'
#' **Which R-hat.** \eqn{\hat R} is the classic Gelman-Rubin potential scale
#' reduction factor (Gelman and Rubin, 1992), without the correction for
#' sampling variability, which is the statistic defined and reported in the
#' paper. The threshold 1.01 is the one recommended for the rank-normalised
#' split \eqn{\hat R} of Vehtari et al. (2021). For the classic statistic,
#' the conventional threshold is 1.1, so the check here is stricter than
#' usual and errs on the side of longer chains.
#'
#' **Settings used in the paper.** The real-data analysis ran four chains of
#' 100,000 iterations per triplet with a burn-in of 2,000. The simulation
#' studies ran one chain per fit, of 20,000 iterations when the instruments
#' are strong (scenarios S1 to S3) and of 400,000 iterations under weak or
#' otherwise misspecified instruments (scenarios S4 to S7), with a burn-in of
#' 2,000. A comparison with MR-CCC should use chains of these lengths, or
#' longer, for the corresponding setting.
#'
#' @return An [MRCCCFit-class] object. See [MRCCCFit-accessors].
#'
#' @seealso [buildTriplet()], [MRCCCFit-accessors], [credibleIntervals()],
#'   [signProbabilities()], [signReversal()], [plotEffectCurve()],
#'   [mrcccScreen()], [bayesFDR()]
#'
#' @references
#' Sarkar B, Ni Y (2026). MR-CCC: Bayesian Mendelian Randomization for
#' Causal Cell-Cell Communication. Bioinformatics (in press).
#' arXiv:2604.23917.
#'
#' Gelman A, Rubin DB (1992). Inference from iterative simulation using
#' multiple sequences. Statistical Science 7(4):457-472.
#'
#' Vehtari A, Gelman A, Simpson D, Carpenter B, Burkner P-C (2021).
#' Rank-normalization, folding, and localization: an improved R-hat for
#' assessing convergence of MCMC. Bayesian Analysis 16(2):667-718.
#'
#' @examples
#' sim <- simulateMRCCC(n = 200, seed = 1)
#' dat <- buildTriplet(sim$sender, sim$receiver, sim$genotypes,
#'                     ligand = "LIG1", receptor = "REC1",
#'                     pathway = sim$pathwayGenes,
#'                     covariates = c("sex", "age", "PC1", "PC2"))
#' fit <- mrccc(dat, nIter = 4000, burnIn = 1000, seed = 1)
#' fit
#' pip(fit)
#' effectSizes(fit)
#'
#' @rdname mrccc
#' @export
setMethod("mrccc", "MRCCCData",
    function(x, nIter = 20000L, burnIn = 2000L, thin = 1L, nChains = 4L,
             initScale = 1, pipThreshold = 0.5,
             hyper = mrcccHyperparameters(), seed = NULL, verbose = FALSE,
             BPPARAM = BiocParallel::SerialParam(), ...) {
    .fitMRCCC(
        X = x@ligandExpression, Z = x@receptorExpression,
        Y = x@pathwayActivity, G = x@ligandGenotypes,
        H = x@receptorGenotypes, V = x@covariates,
        nIter = nIter, burnIn = burnIn, thin = thin, nChains = nChains,
        initScale = initScale, pipThreshold = pipThreshold, hyper = hyper,
        seed = seed, verbose = verbose, BPPARAM = BPPARAM,
        labels = x@labels)
})

# The fitting engine. Arguments are numeric vectors (X, Z, Y) and matrices
# (G, H, V) already aligned by donor, as guaranteed by the MRCCCData validity
# method. Each step is a small helper below.
.fitMRCCC <- function(X, Z, Y, G, H, V, nIter, burnIn, thin, nChains,
                      initScale, pipThreshold, hyper, seed, verbose,
                      BPPARAM = BiocParallel::SerialParam(),
                      labels = .makeLabels(NULL)) {
    colNames <- list(ligand = colnames(G), receptor = colnames(H),
                     covariates = colnames(V))
    inp <- .prepareInputs(X, Z, Y, G, H, V)
    set <- .checkFitArgs(nIter, burnIn, thin, nChains, initScale,
                         pipThreshold, hyper, seed, verbose, BPPARAM)
    cen <- .centreInputs(inp)
    instruments <- .instrumentStrength(cen)
    gVal <- if (is.null(hyper$g)) min(nrow(inp$X), 100) else hyper$g
    fits <- .runChains(cen, set, hyper, gVal, seed, BPPARAM)
    pooled <- .poolChains(fits, set$nChains)
    .convergenceMessages(pooled, set)

    methods::new("MRCCCFit",
        pip = pooled$pip, pipMCSE = pooled$mcse, pipESS = pooled$ess,
        draws = pooled$draws,
        scale = cen$scale,
        diagnostics = pooled$diagnostics,
        instruments = instruments,
        firstStage = .firstStageSummary(fits, colNames),
        posteriorMeans = .posteriorMeansSummary(fits, colNames),
        settings = list(nIter = set$nIter, burnIn = set$burnIn,
                        thin = set$thin, nChains = set$nChains,
                        nKeep = pooled$nKeep,
                        pipThreshold = set$pipThreshold, hyper = hyper,
                        g = gVal, initScale = set$initScale,
                        priorScale = fits[[1L]]$prior_scale),
        labels = labels)
}

# Coerce the inputs and check that the designs can be fitted.
.prepareInputs <- function(X, Z, Y, G, H, V) {
    inp <- list(X = as_column(X, "ligandExpression"),
                Z = as_column(Z, "receptorExpression"),
                Y = as_column(Y, "pathwayActivity"),
                G = as_design(G, "ligandGenotypes"),
                H = as_design(H, "receptorGenotypes"),
                V = as_design(V, "covariates"))
    check_design(inp$G, "ligandGenotypes")
    check_design(inp$H, "receptorGenotypes")
    check_design(inp$V, "covariates")
    n <- nrow(inp$X)
    pTotal <- ncol(inp$G) + ncol(inp$H) + ncol(inp$V)
    if (n <= pTotal + 4L) {
        stop("Too few donors for the model: ", n, " donors against ",
             pTotal, " instrument and covariate columns plus four ",
             "second-stage coefficients. Supply more donors or fewer ",
             "instruments.", call. = FALSE)
    }
    inp
}

# Validate the sampler settings and return them as integers.
.checkFitArgs <- function(nIter, burnIn, thin, nChains, initScale,
                          pipThreshold, hyper, seed, verbose, BPPARAM) {
    for (nm in c("nIter", "burnIn", "thin", "nChains", "initScale",
                 "pipThreshold")) {
        check_scalar(get(nm), nm)
    }
    if (nIter != round(nIter) || nIter < 1) {
        stop("'nIter' must be a positive integer.", call. = FALSE)
    }
    if (burnIn != round(burnIn) || burnIn < 0) {
        stop("'burnIn' must be a non-negative integer.", call. = FALSE)
    }
    if (nIter <= burnIn) {
        stop("'nIter' (", nIter, ") must be greater than 'burnIn' (",
             burnIn, ").", call. = FALSE)
    }
    if (thin != round(thin) || thin < 1) {
        stop("'thin' must be an integer greater than or equal to 1.",
             call. = FALSE)
    }
    if (nChains != round(nChains) || nChains < 1) {
        stop("'nChains' must be an integer greater than or equal to 1.",
             call. = FALSE)
    }
    if (initScale < 0) {
        stop("'initScale' must be non-negative.", call. = FALSE)
    }
    if (nChains > 1 && initScale == 0) {
        stop("'initScale' must be positive when nChains > 1, so that the ",
             "chains start from different values.", call. = FALSE)
    }
    if (pipThreshold <= 0 || pipThreshold >= 1) {
        stop("'pipThreshold' must lie strictly between 0 and 1.",
             call. = FALSE)
    }
    if (!methods::is(hyper, "mrcccHyperparameters")) {
        stop("'hyper' must be created by mrcccHyperparameters().",
             call. = FALSE)
    }
    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
    }
    if (!is.null(seed)) {
        check_scalar(seed, "seed")
        if (seed != round(seed) || abs(seed) > .Machine$integer.max) {
            stop("'seed' must be a whole number no larger in absolute value ",
                 "than .Machine$integer.max.", call. = FALSE)
        }
    }
    if (!methods::is(BPPARAM, "BiocParallelParam")) {
        stop("'BPPARAM' must be a BiocParallelParam object, for example ",
             "BiocParallel::SerialParam().", call. = FALSE)
    }
    nKeepExpected <- (nIter - burnIn) %/% thin
    if (nKeepExpected < 2L) {
        stop("The chosen 'nIter' (", nIter, "), 'burnIn' (", burnIn,
             ") and 'thin' (", thin, ") retain only ", nKeepExpected,
             " draw(s) per chain. At least two are required.",
             call. = FALSE)
    }
    list(nIter = as.integer(nIter), burnIn = as.integer(burnIn),
         thin = as.integer(thin), nChains = as.integer(nChains),
         initScale = initScale, pipThreshold = pipThreshold,
         verbose = verbose)
}

# Centre every input (neither first stage has an intercept) and record the
# scale used for standardised effects.
.centreInputs <- function(inp) {
    sds <- vapply(inp[c("X", "Z", "Y")], function(m) stats::sd(m[, 1L]),
                  numeric(1))
    if (any(!is.finite(sds)) || any(sds <= 0)) {
        stop("Ligand expression, receptor expression and pathway activity ",
             "must each vary across donors.", call. = FALSE)
    }
    means <- vapply(inp[c("X", "Z", "Y")], function(m) mean(m[, 1L]),
                    numeric(1))
    out <- list(X = inp$X - means[["X"]], Z = inp$Z - means[["Z"]],
                Y = inp$Y - means[["Y"]], G = center_columns(inp$G),
                H = center_columns(inp$H), V = center_columns(inp$V),
                Graw = inp$G, Hraw = inp$H, Vraw = inp$V)
    check_rank(cbind(out$G, out$V), "ligand instruments and covariates")
    check_rank(cbind(out$H, out$V), "receptor instruments and covariates")
    out$scale <- list(sdX = sds[["X"]], sdZ = sds[["Z"]], sdY = sds[["Y"]],
                      meanX = means[["X"]], meanZ = means[["Z"]],
                      meanY = means[["Y"]],
                      zRangeStd = range(out$Z[, 1L] / sds[["Z"]]))
    out
}

# First-stage partial F statistics, with a message for weak instruments.
.instrumentStrength <- function(cen) {
    instruments <- data.frame(
        exposure = c("ligand", "receptor"),
        nInstruments = c(ncol(cen$G), ncol(cen$H)),
        F = c(first_stage_F(cen$X[, 1L], cen$Graw, cen$Vraw),
              first_stage_F(cen$Z[, 1L], cen$Hraw, cen$Vraw)),
        stringsAsFactors = FALSE)
    instruments$weak <- is.finite(instruments$F) & instruments$F < 10
    if (any(instruments$weak)) {
        wk <- instruments[instruments$weak, , drop = FALSE]
        weakText <- sprintf("%s (F = %s)", wk$exposure,
                            vapply(wk$F, fmt, character(1)))
        message("First-stage F below 10 for the ",
                gsub(", ", " and ", toString(weakText), fixed = TRUE),
                ": the instruments are weak by the two-stage least-squares ",
                "heuristic. MR-CCC does not require F > 10, but posterior ",
                "inclusion probabilities are less decisive in this regime ",
                "and longer chains may be needed. See instrumentStrength().")
    }
    instruments
}

# One seed per chain, derived from 'seed' or drawn from the current stream.
.chainSeeds <- function(nChains, seed) {
    draw <- function() sample.int(.Machine$integer.max, nChains)
    if (is.null(seed)) {
        draw()
    } else {
        withr::with_seed(as.integer(seed), draw(),
                         .rng_kind = "Mersenne-Twister",
                         .rng_normal_kind = "Inversion",
                         .rng_sample_kind = "Rejection")
    }
}

# Run the chains, each under its own seed, through BiocParallel.
.runChains <- function(cen, set, hyper, gVal, seed, BPPARAM) {
    seeds <- .chainSeeds(set$nChains, seed)
    initRaw <- if (set$nChains > 1L) set$initScale else 0
    runOne <- function(cc) {
        withr::with_seed(seeds[[cc]], .rng_kind = "Mersenne-Twister",
                         .rng_normal_kind = "Inversion",
                         .rng_sample_kind = "Rejection", mr_ccc_gibbs(
            cen$X, cen$Z, cen$Y, cen$G, cen$H, cen$V,
            n_iter = set$nIter, burn_in = set$burnIn, thin = set$thin,
            a_sigma = hyper$aSigma, b_sigma = hyper$bSigma,
            a_rho = hyper$aRho, b_rho = hyper$bRho, nu1 = hyper$nu1,
            gG = gVal, gV = gVal, gH = gVal, gZ = gVal, gBeta = gVal,
            ridge = hyper$ridge,
            init_gamma = if (cc %% 2L == 0L) 0L else 1L,
            init_scale = initRaw, verbose = set$verbose,
            save_first_stage = TRUE))
    }
    # The caller's random number state is left as it was, whatever the back
    # end does with the stream.
    withr::with_preserve_seed(
        BiocParallel::bplapply(seq_len(set$nChains), runOne,
                               BPPARAM = BPPARAM))
}

# Pool the retained draws and compute the PIP, its Monte Carlo error and the
# Gelman-Rubin statistics.
.poolChains <- function(fits, nChains) {
    nKeep <- as.integer(fits[[1L]]$n_keep)
    chainOf <- function(name) lapply(fits, function(f) as.numeric(f[[name]]))
    ch <- list(beta_X = chainOf("Beta_X_draws"),
               beta_XZ = chainOf("Beta_XZ_draws"),
               gamma = chainOf("gamma_draws"),
               loglik = chainOf("loglik_draws"))
    draws <- data.frame(
        chain = rep(seq_len(nChains), each = nKeep),
        beta_X = unlist(ch$beta_X, use.names = FALSE),
        beta_XZ = unlist(ch$beta_XZ, use.names = FALSE),
        beta_Z = unlist(chainOf("Beta_Z_draws"), use.names = FALSE),
        gamma = unlist(ch$gamma, use.names = FALSE),
        mu = unlist(chainOf("mu_draws"), use.names = FALSE),
        loglik = unlist(ch$loglik, use.names = FALSE))
    mc <- mcse_pooled(lapply(ch$gamma, mcse_batch))
    rhat <- if (nChains >= 2L) {
        vapply(ch, gelman_rubin, numeric(1))
    } else {
        rep(NA_real_, length(ch))
    }
    list(draws = draws, nKeep = nKeep, nChains = nChains,
         pip = mean(draws$gamma), mcse = unname(mc[["mcse"]]),
         ess = unname(mc[["ess"]]),
         constantGamma = length(unique(draws$gamma)) == 1L,
         diagnostics = data.frame(parameter = names(ch), rhat = unname(rhat),
                                  stringsAsFactors = FALSE))
}

# Warnings and messages about convergence and Monte Carlo error.
.convergenceMessages <- function(pooled, set) {
    rh <- pooled$diagnostics$rhat
    bad <- !is.na(rh) & (is.infinite(rh) | rh > 1.01)
    if (any(bad)) {
        warning("Gelman-Rubin R-hat exceeds 1.01 for ",
                toString(pooled$diagnostics$parameter[bad]),
                ": the chains may not have mixed. Increase 'nIter' and ",
                "'burnIn'.",
                call. = FALSE)
    }
    if (pooled$nChains == 1L) {
        message("Only one chain was run, so R-hat is unavailable. Use ",
                "nChains >= 2 to assess convergence.")
    }
    # When every chain holds the indicator at the same value throughout, the
    # PIP is 0 or 1 with no observed Monte Carlo variation and the MCSE is
    # reported as NA without comment.
    if (!is.finite(pooled$mcse) && !pooled$constantGamma) {
        message("The Monte Carlo uncertainty of the PIP could not be ",
                "estimated (too few retained draws per chain, or a constant ",
                "indicator chain), so the borderline check was skipped.")
    }
    if (is.finite(pooled$mcse) &&
        abs(pooled$pip - set$pipThreshold) < 2 * pooled$mcse) {
        # The paper pools about 400,000 draws per triplet: four chains of
        # 100,000 iterations, or one chain of 400,000.
        target <- as.integer(ceiling(400000 / pooled$nChains / 1000) * 1000)
        suggestion <- if (set$nIter < target) {
            paste0("Increase 'nIter' to ", target, ", so that the ",
                   pooled$nChains, " chain(s) pool about 400000 draws as ",
                   "in the paper")
        } else {
            paste("Increase 'nChains', which pools additional draws at",
                  "proportional cost")
        }
        warning("The posterior inclusion probability (", fmt(pooled$pip),
                ") lies within 2 * MCSE (MCSE = ", fmt(pooled$mcse),
                ") of the discovery threshold (", set$pipThreshold,
                "): the discovery status of this triplet is not resolved at ",
                "this chain length. ", suggestion, ".", call. = FALSE)
    }
    invisible(NULL)
}

# Posterior mean and 95% interval of each instrument's first-stage effect,
# pooled over chains.
.firstStageSummary <- function(fits, colNames) {
    one <- function(field, snps, exposure) {
        D <- do.call(rbind, lapply(fits, function(f) {
            as.matrix(f[[field]])
        }))
        data.frame(
            exposure = exposure,
            snp = if (is.null(snps)) {
                paste0(exposure, "SNP", seq_len(ncol(D)))
            } else {
                snps
            },
            mean = colMeans(D),
            lower = apply(D, 2L, stats::quantile, probs = 0.025,
                          names = FALSE),
            upper = apply(D, 2L, stats::quantile, probs = 0.975,
                          names = FALSE),
            row.names = NULL, stringsAsFactors = FALSE)
    }
    rbind(one("Pi_X_draws", colNames$ligand, "ligand"),
          one("Pi_Z_draws", colNames$receptor, "receptor"))
}

# Posterior means of the nuisance parameters, averaged over chains.
.posteriorMeansSummary <- function(fits, colNames) {
    avg <- function(field) {
        Reduce(`+`, lapply(fits, function(f) as.numeric(f[[field]]))) /
            length(fits)
    }
    cov <- colNames$covariates
    if (is.null(cov)) cov <- paste0("covariate", seq_along(avg("Alpha_Y_mean")))
    list(
        covariateEffects = data.frame(
            covariate = cov, ligand = avg("Alpha_X_mean"),
            receptor = avg("Alpha_Z_mean"), pathway = avg("Alpha_Y_mean"),
            row.names = NULL, stringsAsFactors = FALSE),
        residualVariances = c(ligand = avg("sigma_X_sq_mean"),
                              receptor = avg("sigma_Z_sq_mean"),
                              pathway = avg("sigma_Y_sq_mean")),
        inclusionProbabilityPrior = avg("rho_mean"),
        intercept = avg("mu_mean"))
}
