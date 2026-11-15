#' Bayesian false discovery rate of nested PIP-ranked sets
#'
#' For a vector of posterior inclusion probabilities from many
#' ligand-receptor-pathway triplets, computes the Bayesian false discovery
#' rate (bFDR) of every nested set obtained by declaring the top-\eqn{k}
#' triplets by PIP, for \eqn{k = 1, \ldots, m}.
#'
#' @param pip Posterior inclusion probabilities, one per triplet: a
#'   numeric vector in \eqn{[0, 1]}, or a list of [MRCCCFit-class] objects.
#'   Missing values are not allowed.
#'
#' @details
#' For a declared set \eqn{S}, the posterior expected number of false
#' discoveries is \eqn{\sum_{j \in S} (1 - \mathrm{PIP}_j)}, so the Bayesian
#' false discovery rate is
#' \deqn{\mathrm{bFDR}(S) = |S|^{-1} \sum_{j \in S} (1 - \mathrm{PIP}_j).}
#' Ranking the triplets by decreasing PIP and accumulating gives, for each
#' possible cut-off, the bFDR incurred by declaring everything at or above
#' it. Because the PIPs are posterior probabilities, this controls the
#' expected proportion of false discoveries without a separate multiplicity
#' correction. The value returned for triplet \eqn{j} is the bFDR of the
#' smallest nested set that contains \eqn{j}. It is non-decreasing in the
#' rank of \eqn{j}. Ties in PIP are broken by input order.
#'
#' @return A numeric vector of the same length as `pip`, in the order of the
#'   input, giving the bFDR of the nested set that ends at each triplet.
#'
#' @seealso [bfdrSelect()], [mrccc()], [pip()]
#'
#' @references
#' Newton MA, Noueiry A, Sarkar D, Ahlquist P (2004). Detecting differential
#' gene expression with a semiparametric hierarchical mixture method.
#' Biostatistics 5(2):155-176.
#'
#' @examples
#' pip <- c(0.9, 0.8, 0.5)
#' bayesFDR(pip)
#' # Unsorted input is handled. Results are returned in input order.
#' bayesFDR(c(0.5, 0.9, 0.8))
#'
#' @export
bayesFDR <- function(pip) {
    pip <- check_pip(pip)
    if (length(pip) == 0L) return(numeric(0))
    o   <- order(pip, decreasing = TRUE)
    f   <- cumsum(1 - pip[o]) / seq_along(pip)
    out <- numeric(length(pip))
    out[o] <- f
    out
}

#' Largest PIP-ranked set with Bayesian FDR at most alpha
#'
#' Selects the largest nested set of triplets, ranked by decreasing posterior
#' inclusion probability, whose Bayesian false discovery rate does not exceed
#' `alpha`. The rule is monotone: if the top-\eqn{k} set qualifies, so does
#' every smaller top-\eqn{k'} set. If no set qualifies, nothing is selected.
#'
#' @param pip Posterior inclusion probabilities, one per triplet: a
#'   numeric vector in \eqn{[0, 1]}, or a list of [MRCCCFit-class] objects.
#' @param alpha Target Bayesian false discovery rate in \eqn{(0, 1)}.
#'
#' @return A logical vector of the same length as `pip`, in input order,
#'   that is `TRUE` for the selected triplets.
#'
#' @seealso [bayesFDR()]
#'
#' @examples
#' pip <- c(0.9, 0.8, 0.5, 0.2)
#' bayesFDR(pip)
#' bfdrSelect(pip, alpha = 0.2)
#' bfdrSelect(pip, alpha = 0.05)
#'
#' @export
bfdrSelect <- function(pip, alpha) {
    pip <- check_pip(pip)
    check_scalar(alpha, "alpha")
    if (alpha <= 0 || alpha >= 1) {
        stop("'alpha' must lie strictly between 0 and 1.", call. = FALSE)
    }
    keep <- rep(FALSE, length(pip))
    if (length(pip) == 0L) return(keep)
    o <- order(pip, decreasing = TRUE)
    f <- cumsum(1 - pip[o]) / seq_along(pip)
    ok <- which(f <= alpha)
    if (length(ok) > 0L) {
        k <- max(ok)
        keep[o[seq_len(k)]] <- TRUE
    }
    keep
}

#' Validate a vector of posterior inclusion probabilities
#'
#' @param pip Object supplied by the user.
#' @return Numeric vector.
#' @noRd
#' @keywords internal
check_pip <- function(pip) {
    if (is.list(pip) && length(pip) > 0L &&
        all(vapply(pip, methods::is, logical(1), class2 = "MRCCCFit"))) {
        pip <- vapply(pip, function(f) f@pip, numeric(1))
    }
    if (!is.numeric(pip)) {
        stop("'pip' must be a numeric vector of posterior inclusion ",
             "probabilities or a list of MRCCCFit objects.", call. = FALSE)
    }
    pip <- as.numeric(pip)
    if (any(!is.finite(pip))) {
        stop("'pip' contains NA, NaN or infinite values.", call. = FALSE)
    }
    if (any(pip < 0 | pip > 1)) {
        stop("'pip' must lie in [0, 1].", call. = FALSE)
    }
    pip
}
