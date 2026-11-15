# Internal helpers for MRCCC. None of these functions is exported.

# ---- Input coercion --------------------------------------------------------

#' Coerce a numeric vector or one-column matrix to an n x 1 matrix
#'
#' @param x Object supplied by the user.
#' @param name Argument name used in error messages.
#' @return Numeric matrix with one column.
#' @noRd
#' @keywords internal
as_column <- function(x, name) {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (!is.numeric(x)) {
        stop("'", name, "' must be numeric (a vector or an n x 1 matrix).",
             call. = FALSE)
    }
    if (is.matrix(x)) {
        if (ncol(x) != 1L) {
            stop("'", name, "' must have exactly one column. It has ",
                 ncol(x), ".",
                 call. = FALSE)
        }
    } else {
        x <- matrix(as.numeric(x), ncol = 1L)
    }
    storage.mode(x) <- "double"
    dimnames(x) <- NULL
    x
}

#' Coerce a design matrix argument to a numeric matrix with at least one column
#'
#' @param x Object supplied by the user.
#' @param name Argument name used in error messages.
#' @return Numeric matrix.
#' @noRd
#' @keywords internal
as_design <- function(x, name) {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (!is.numeric(x)) {
        stop("'", name, "' must be a numeric matrix.", call. = FALSE)
    }
    if (!is.matrix(x)) x <- matrix(as.numeric(x), ncol = 1L)
    if (ncol(x) < 1L) {
        stop("'", name, "' must have at least one column.", call. = FALSE)
    }
    storage.mode(x) <- "double"
    dimnames(x) <- NULL
    x
}

#' Check that an argument is a single finite number
#'
#' @param x Value to check.
#' @param name Argument name used in error messages.
#' @return Invisibly, `x`.
#' @noRd
#' @keywords internal
check_scalar <- function(x, name) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
        stop("'", name, "' must be a single finite number.", call. = FALSE)
    }
    invisible(x)
}

#' Check that a design matrix is usable by the sampler
#'
#' Three conditions are checked, each of which produces a fitted model that is
#' numerically defined but scientifically meaningless, and each of which is
#' easier to diagnose here than from an implausible posterior:
#'
#' * a CONSTANT column. In `V` this duplicates the model's own intercept
#'   `mu`. In `G` or `H` it is a monomorphic variant carrying no information
#'   about the exposure. Either way the column contributes nothing and makes
#'   the cross-product matrix singular up to the ridge.
#' * DUPLICATED columns, which are perfectly collinear. Two identical
#'   genotype columns arise easily when a variant is listed twice in a
#'   summary file.
#' * MORE COLUMNS THAN ROWS, which saturates the corresponding regression.
#'
#' @param x Numeric matrix, already coerced by `as_design()`.
#' @param name Argument name used in error messages.
#' @return Invisibly, `x`.
#' @noRd
#' @keywords internal
check_design <- function(x, name) {
    const <- which(apply(x, 2L, function(cl) {
        v <- stats::var(cl)
        !is.finite(v) || v <= 0
    }))
    if (length(const) > 0L) {
        extra <- if (identical(name, "covariates")) {
            paste0(" The model already includes an intercept, so a constant ",
                   "covariate is redundant. Remove it.")
        } else {
            paste0(" A constant column in '", name, "' is a monomorphic ",
                   "variant ",
                   "and carries no information about the exposure. Remove it.")
        }
        stop("'", name, "' has constant column(s) at position(s) ",
             toString(const), ".", extra, call. = FALSE)
    }

    if (ncol(x) > 1L) {
        dup <- which(duplicated(t(x)))
        if (length(dup) > 0L) {
            stop("'", name, "' has duplicated column(s) at position(s) ",
                 toString(dup),
                 ", which are perfectly collinear with an earlier column. ",
                 "Remove the repeats.", call. = FALSE)
        }
    }

    if (ncol(x) >= nrow(x)) {
        stop("'", name, "' has ", ncol(x), " columns and only ", nrow(x),
             " rows. The corresponding regression is saturated. Supply fewer ",
             "columns than donors.", call. = FALSE)
    }

    invisible(x)
}

#' Check that an object is an MRCCCFit
#'
#' @param fit Object to check.
#' @return Invisibly, `fit`.
#' @noRd
#' @keywords internal
check_fit <- function(fit) {
    if (!methods::is(fit, "MRCCCFit")) {
        stop("'x' must be an MRCCCFit object returned by mrccc().",
             call. = FALSE)
    }
    invisible(fit)
}

#' Standardisation factors for the three effect sizes
#'
#' The raw coefficients are expressed per unit of the centred inputs. The
#' standardised coefficients are expressed per standard deviation of X (and
#' Z) in standard-deviation units of Y.
#'
#' @param fit An `MRCCCFit` object.
#' @return Named numeric vector with elements `beta_X`, `beta_XZ`, `beta_Z`.
#' @noRd
#' @keywords internal
std_factors <- function(fit) {
    s <- fit@scale
    c(beta_X  = s$sdX / s$sdY,
      beta_XZ = s$sdX * s$sdZ / s$sdY,
      beta_Z  = s$sdZ / s$sdY)
}

# ---- Centring and rank -----------------------------------------------------

#' Subtract each column's mean
#'
#' The design matrices G, H and V enter the first stages without an
#' intercept, so they must be centred like X, Z and Y.
#'
#' @param M Numeric matrix.
#' @return `M` with every column centred.
#' @noRd
#' @keywords internal
center_columns <- function(M) {
    M - matrix(colMeans(M), nrow(M), ncol(M), byrow = TRUE)
}

#' Stop when a centred design matrix is rank deficient
#'
#' @param W Numeric matrix.
#' @param name Name used in the error message.
#' @return `invisible(TRUE)`. The function is called for its error.
#' @noRd
#' @keywords internal
check_rank <- function(W, name) {
    r <- qr(W)$rank
    if (r < ncol(W)) {
        stop("'", name, "' is rank deficient after centring (",
             ncol(W) - r, " redundant column(s)). For example, one variant ",
             "coded on both alleles, or a full set of dummy variables. ",
             "Remove the redundant columns.", call. = FALSE)
    }
    invisible(TRUE)
}

# ---- Instrument strength ---------------------------------------------------

#' First-stage partial F statistic
#'
#' F test of the instrument block after the covariates, from ordinary least
#' squares: the residual sum of squares of `y ~ 1 + V` is compared with that
#' of `y ~ 1 + V + W`. This is the conventional measure of instrument
#' strength. The threshold F > 10 is a two-stage least-squares heuristic and
#' is reported for orientation rather than enforced.
#'
#' @param y Numeric vector: the exposure.
#' @param W Numeric matrix: instruments (n x p_W).
#' @param V Numeric matrix: covariates (n x p_V).
#' @return A single numeric value, or `NA` if the residual degrees of freedom
#'   are not positive or either fit is degenerate.
#' @noRd
#' @keywords internal
first_stage_F <- function(y, W, V) {
    n   <- length(y)
    p_V <- ncol(V)
    p_W <- ncol(W)
    df2 <- n - 1L - p_V - p_W
    if (df2 <= 0L) return(NA_real_)

    one  <- rep(1, n)
    rss0 <- sum(stats::lm.fit(cbind(one, V),    y)$residuals^2)
    rss1 <- sum(stats::lm.fit(cbind(one, V, W), y)$residuals^2)
    if (!is.finite(rss0) || !is.finite(rss1) || rss1 <= 0) return(NA_real_)

    ((rss0 - rss1) / p_W) / (rss1 / df2)
}

# ---- Monte Carlo error and convergence diagnostics -------------------------

#' Batch-means Monte Carlo standard error and effective sample size
#'
#' Non-overlapping batch means for a single chain. The chain is divided into
#' `n_batch` consecutive batches. The variance of the batch means, scaled by
#' the batch length, estimates the asymptotic variance of the chain mean.
#'
#' @param x Numeric vector: one chain.
#' @param n_batch Number of batches.
#' @return Named numeric vector with elements `mcse` and `ess`. Both are
#'   `NA` when the chain is too short or constant.
#' @noRd
#' @keywords internal
mcse_batch <- function(x, n_batch = 30L) {
    x <- as.numeric(x)
    n <- length(x)
    b <- floor(n / n_batch)
    if (n < 4L || b < 2L) return(c(mcse = NA_real_, ess = NA_real_))
    k  <- floor(n / b)
    bm <- vapply(seq_len(k),
                 function(i) mean(x[((i - 1L) * b + 1L):(i * b)]),
                 numeric(1))
    sigma2 <- b * stats::var(bm)
    vx     <- stats::var(x)
    if (!is.finite(sigma2) || sigma2 <= 0) {
        return(c(mcse = NA_real_, ess = NA_real_))
    }
    c(mcse = sqrt(sigma2 / n), ess = n * vx / sigma2)
}

#' Combine per-chain batch-means estimates across independent chains
#'
#' The pooled mean is the average of the chain means, so its variance is the
#' sum of the per-chain variances divided by the squared number of chains:
#' `mcse_pooled = sqrt(sum(mcse_c^2)) / C`. Effective sample sizes add.
#' Chains whose per-chain estimate is `NA` (constant chains, whose true Monte
#' Carlo error is zero) contribute zero to the numerator but are still
#' counted in `C`.
#'
#' @param mc_list List of outputs of `mcse_batch()`.
#' @return Named numeric vector with elements `mcse` and `ess`.
#' @noRd
#' @keywords internal
mcse_pooled <- function(mc_list) {
    ok <- vapply(mc_list, function(z) is.finite(z[["mcse"]]), logical(1))
    if (!any(ok)) return(c(mcse = NA_real_, ess = NA_real_))
    se <- vapply(mc_list[ok], function(z) z[["mcse"]], numeric(1))
    es <- vapply(mc_list[ok], function(z) z[["ess"]],  numeric(1))
    c(mcse = sqrt(sum(se^2)) / length(mc_list), ess = sum(es))
}

#' Gelman-Rubin potential scale reduction factor
#'
#' For `m` chains of equal length `n`, `W` is the mean within-chain variance
#' and `B` is `n` times the variance of the chain means. The statistic is
#' `sqrt((((n - 1) / n) * W + B / n) / W)`. When every chain is constant
#' (`W == 0`) the result is `Inf` if the chains sit at different values
#' (`B > 0`) and `NA` if they are identical.
#'
#' @param chains List of numeric vectors of equal length.
#' @return A single numeric value, or `NA` when fewer than two chains are
#'   supplied.
#' @noRd
#' @keywords internal
gelman_rubin <- function(chains) {
    m <- length(chains)
    if (m < 2L) return(NA_real_)
    n <- length(chains[[1L]])
    if (n < 2L) return(NA_real_)
    means <- vapply(chains, mean, numeric(1))
    vars  <- vapply(chains, stats::var,  numeric(1))
    W <- mean(vars)
    B <- n * stats::var(means)
    if (!is.finite(W)) return(NA_real_)
    if (W <= 0) return(if (is.finite(B) && B > 0) Inf else NA_real_)
    var_hat <- ((n - 1) / n) * W + B / n
    sqrt(var_hat / W)
}

# ---- Formatting ------------------------------------------------------------

#' Format a number for console output
#'
#' @param x Numeric value.
#' @param digits Number of significant digits.
#' @return Character string.
#' @noRd
#' @keywords internal
fmt <- function(x, digits = 3L) {
    if (length(x) != 1L || is.na(x)) return("NA")
    if (is.infinite(x)) return(if (x > 0) "Inf" else "-Inf")
    as.character(signif(x, digits))
}
