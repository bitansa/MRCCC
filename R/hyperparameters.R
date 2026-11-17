#' Hyperparameters of the MR-CCC prior
#'
#' Constructs and validates the list of prior hyperparameters used by
#' [mrccc()]. The defaults are those used throughout the MR-CCC analyses and
#' are weakly informative for inputs on a scale of roughly unit variance:
#' with several hundred donors the likelihood then dominates the
#' inverse-gamma and beta priors. The inverse-gamma priors are not invariant
#' to the units of `X`, `Z` and `Y`, so inputs on very different scales
#' should be standardised before fitting.
#'
#' @param aSigma,bSigma Shape and scale of the inverse-gamma prior on the
#'   three residual variances \eqn{\sigma_X^2, \sigma_Z^2, \sigma_Y^2}. Both
#'   must be positive.
#' @param aRho,bRho Shape parameters of the beta prior on the prior
#'   inclusion probability \eqn{\rho}. Both must be positive.
#' @param nu1 Spike variance multiplier in \eqn{(0, 1)}. The spike component
#'   of the prior on \eqn{(\beta_X, \beta_{XZ})} has the slab covariance
#'   scaled by `nu1`, so small values enforce a near-zero spike.
#' @param g Common g-prior scale applied to every coefficient block
#'   (first-stage sender effects, first-stage receiver effects, covariate
#'   effects, receptor main effect, and the causal block). If `NULL` (the
#'   default), [mrccc()] sets \eqn{g = \min(n, 100)}, where \eqn{n} is the
#'   number of donors. A positive number overrides this rule.
#' @param ridge Positive diagonal ridge added before every matrix inversion
#'   for numerical stability.
#'
#' @return A named list of class `mrcccHyperparameters` with elements
#'   `aSigma`, `bSigma`, `aRho`, `bRho`, `nu1`, `g` and `ridge`.
#'
#' @details
#' Each first-stage block \eqn{\theta} with observed design \eqn{D} (the
#' genotype matrices \eqn{G}, \eqn{H} or the covariates \eqn{V}) has the
#' Zellner g-prior \eqn{\theta \sim N(0, g\,\sigma^2 (D^\top D)^{-1})}. The
#' second-stage blocks \eqn{(\beta_X, \beta_{XZ})} and \eqn{\beta_Z} use the
#' same scale \eqn{g} but against a fixed diagonal computed once from the
#' least-squares first stage, because their regressors \eqn{X^*} and
#' \eqn{Z^*} are themselves parameters of the model. See the *Priors* section
#' of [mrccc()]. The default \eqn{g = \min(n, 100)} follows the
#' unit-information convention for small samples and caps the prior
#' dispersion for large ones. The spike-and-slab prior on the causal block
#' uses scale \eqn{g} in the slab and \eqn{\nu_1 g} in the spike.
#'
#' @seealso [mrccc()]
#'
#' @examples
#' hyper <- mrcccHyperparameters()
#' str(hyper)
#'
#' # A wider spike (less strict shrinkage under gamma = 0) and a fixed
#' # g-prior scale of 50, which narrows the slab when there are more than
#' # 50 donors
#' mrcccHyperparameters(nu1 = 1e-3, g = 50)
#'
#' @export
mrcccHyperparameters <- function(aSigma = 3, bSigma = 2, aRho = 3, bRho = 1,
                                 nu1 = 1e-4, g = NULL, ridge = 1e-8) {
    check_scalar(aSigma, "aSigma")
    check_scalar(bSigma, "bSigma")
    check_scalar(aRho, "aRho")
    check_scalar(bRho, "bRho")
    check_scalar(nu1, "nu1")
    check_scalar(ridge, "ridge")
    if (aSigma <= 0 || bSigma <= 0) {
        stop("'aSigma' and 'bSigma' must be positive.", call. = FALSE)
    }
    if (aRho <= 0 || bRho <= 0) {
        stop("'aRho' and 'bRho' must be positive.", call. = FALSE)
    }
    if (nu1 <= 0 || nu1 >= 1) {
        stop("'nu1' must lie strictly between 0 and 1.", call. = FALSE)
    }
    if (ridge <= 0) {
        stop("'ridge' must be positive.", call. = FALSE)
    }
    if (!is.null(g)) {
        check_scalar(g, "g")
        if (g <= 0) stop("'g' must be positive or NULL.", call. = FALSE)
    }
    structure(
        list(aSigma = as.numeric(aSigma), bSigma = as.numeric(bSigma),
             aRho = as.numeric(aRho), bRho = as.numeric(bRho),
             nu1 = as.numeric(nu1),
             g = if (is.null(g)) NULL else as.numeric(g),
             ridge = as.numeric(ridge)),
        class = "mrcccHyperparameters")
}
