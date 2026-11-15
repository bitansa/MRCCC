#' @include AllClasses.R AllGenerics.R
NULL

#' Create an MRCCCData object from aligned donor-level data
#'
#' Constructs the input object for one ligand-receptor-pathway triplet from
#' donor-level vectors and matrices. Donors are matched by name rather than
#' by position, so components supplied in a different donor order are
#' aligned automatically and components that refer to different donors are
#' rejected. Most users will create the object
#' from `SummarizedExperiment` data with [buildTriplet()] instead.
#'
#' @param ligandExpression Numeric vector (or one-column matrix), one value
#'   per donor: ligand expression in the sender cell type, typically
#'   normalised pseudo-bulk expression.
#' @param receptorExpression Numeric vector: receptor expression in the
#'   receiver cell type.
#' @param pathwayActivity Numeric vector: activity score of the downstream
#'   pathway in the receiver cell type.
#' @param ligandGenotypes Numeric matrix or data frame, donors in rows and
#'   SNPs in columns: allele dosages (0, 1 or 2 copies of the alternative
#'   allele, or an imputed dosage between 0 and 2) of the cis-eQTL SNPs that
#'   instrument the ligand. See [selectInstruments()].
#' @param receptorGenotypes Numeric matrix or data frame, donors by SNPs:
#'   allele dosages of the cis-eQTL SNPs that instrument the receptor.
#' @param covariates Numeric matrix or data frame, donors in rows: donor
#'   covariates such as sex, age and genotype principal components. Factor,
#'   character and logical columns of a data frame are expanded to indicator
#'   columns (treatment contrasts). At least one covariate is required.
#' @param labels Optional named character vector or list with any of the
#'   elements `ligand`, `receptor`, `pathway`, `sender` and `receiver`, used
#'   in printed output.
#' @param metadata Optional list of additional information.
#'
#' @details
#' **Donor matching.** If every component carries donor names (names of the
#' vectors, row names of the matrices), the donors of `ligandExpression`
#' define the order and every other component is reordered to match. A
#' component that lacks any of these donors is an error. If no component
#' carries donor names, all must have the same length and are taken to be
#' in the same order, and donors are named `donor1`, `donor2`, and so on.
#' Supplying names for some components but not others is an error.
#'
#' **Units.** Genotypes are individual-level allele dosages of the same
#' donors whose expression is analysed. MR-CCC is a one-sample method and
#' does not use external eQTL summary statistics. Expression and pathway
#' activity may be on any scale (for example log-normalised pseudo-bulk
#' expression). [mrccc()] centres all inputs and reports effects on both
#' the centred input scale and a standardised scale.
#'
#' @return An [MRCCCData-class] object.
#'
#' @seealso [buildTriplet()], [mrccc()], [MRCCCData-accessors]
#'
#' @examples
#' set.seed(1)
#' n <- 50
#' ids <- sprintf("D%02d", seq_len(n))
#' G <- matrix(rbinom(n * 3, 2, 0.3), n, 3,
#'             dimnames = list(ids, c("rs1", "rs2", "rs3")))
#' H <- matrix(rbinom(n * 2, 2, 0.4), n, 2,
#'             dimnames = list(ids, c("rs4", "rs5")))
#' cov <- data.frame(sex = sample(c("female", "male"), n, TRUE),
#'                   age = round(runif(n, 20, 80)), row.names = ids)
#' x <- setNames(drop(G %*% c(0.5, 0.3, 0.2)) + rnorm(n), ids)
#' z <- setNames(drop(H %*% c(0.4, 0.4)) + rnorm(n), ids)
#' y <- setNames(0.3 * x + 0.2 * x * z + rnorm(n), ids)
#' dat <- MRCCCData(ligandExpression = x, receptorExpression = z,
#'                  pathwayActivity = y, ligandGenotypes = G,
#'                  receptorGenotypes = H, covariates = cov,
#'                  labels = c(ligand = "LIG1", receptor = "REC1"))
#' dat
#'
#' @export
MRCCCData <- function(ligandExpression, receptorExpression, pathwayActivity,
                      ligandGenotypes, receptorGenotypes, covariates,
                      labels = NULL, metadata = list()) {
    req <- c(ligandExpression = missing(ligandExpression),
             receptorExpression = missing(receptorExpression),
             pathwayActivity = missing(pathwayActivity),
             ligandGenotypes = missing(ligandGenotypes),
             receptorGenotypes = missing(receptorGenotypes),
             covariates = missing(covariates))
    if (any(req)) {
        stop("Argument(s) ", toString(names(req)[req]),
             " must be supplied.", call. = FALSE)
    }
    args <- list(ligandExpression = ligandExpression,
                 receptorExpression = receptorExpression,
                 pathwayActivity = pathwayActivity,
                 ligandGenotypes = ligandGenotypes,
                 receptorGenotypes = receptorGenotypes,
                 covariates = covariates)
    missingArgs <- vapply(names(args), function(nm) {
        is.null(args[[nm]])
    }, logical(1))
    if (any(missingArgs)) {
        stop("Argument(s) ", toString(names(args)[missingArgs]),
             " must be supplied.", call. = FALSE)
    }

    vecs <- lapply(names(args)[seq_len(3L)],
                   function(nm) .asDonorVector(args[[nm]], nm))
    names(vecs) <- names(args)[seq_len(3L)]
    G <- .asDonorMatrix(ligandGenotypes, "ligandGenotypes", "ligandSNP")
    H <- .asDonorMatrix(receptorGenotypes, "receptorGenotypes", "receptorSNP")
    V <- .covariateMatrix(covariates)

    # Donor alignment by name, or by position when no component is named.
    ids <- list(names(vecs[[1]]), names(vecs[[2]]), names(vecs[[3]]),
                rownames(G), rownames(H), rownames(V))
    names(ids) <- names(args)
    hasIds <- !vapply(ids, is.null, logical(1))
    if (all(hasIds)) {
        donorIds <- ids[[1]]
        if (anyDuplicated(donorIds) > 0L) {
            stop("Donor names of 'ligandExpression' are not unique.",
                 call. = FALSE)
        }
        for (nm in names(ids)[-1]) {
            miss <- setdiff(donorIds, ids[[nm]])
            if (length(miss)) {
                stop("'", nm, "' lacks ", length(miss), " of the ",
                     length(donorIds), " donors of 'ligandExpression' ",
                     "(for example '", miss[1], "'). Supply all components ",
                     "for the same donors, or use buildTriplet(), which ",
                     "matches donors across objects.", call. = FALSE)
            }
        }
        vecs <- lapply(vecs, function(v) v[donorIds])
        G <- G[donorIds, , drop = FALSE]
        H <- H[donorIds, , drop = FALSE]
        V <- V[donorIds, , drop = FALSE]
    } else if (!any(hasIds)) {
        n <- length(vecs[[1]])
        lens <- c(vapply(vecs, length, integer(1)), nrow(G), nrow(H), nrow(V))
        if (any(lens != n)) {
            found <- toString(sprintf("%s = %d", names(args), lens))
            stop("Without donor names all components must have the same ",
                 "number of donors. Found ", found, ".", call. = FALSE)
        }
        donorIds <- paste0("donor", seq_len(n))
        vecs <- lapply(vecs, function(v) stats::setNames(v, donorIds))
        rownames(G) <- rownames(H) <- rownames(V) <- donorIds
    } else {
        stop("Donor names are present for ",
             toString(names(ids)[hasIds]), " but missing for ",
             toString(names(ids)[!hasIds]),
             ". Name the donors of every component, or of none.",
             call. = FALSE)
    }

    methods::new("MRCCCData",
        ligandExpression   = vecs$ligandExpression,
        receptorExpression = vecs$receptorExpression,
        pathwayActivity    = vecs$pathwayActivity,
        ligandGenotypes    = G,
        receptorGenotypes  = H,
        covariates         = V,
        labels             = .makeLabels(labels),
        metadata           = as.list(metadata))
}

# ---- Internal coercion helpers ---------------------------------------------

.asDonorVector <- function(x, arg) {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.matrix(x)) {
        if (ncol(x) != 1L) {
            stop("'", arg, "' must be a vector or a one-column matrix. It has ",
                 ncol(x), " columns.", call. = FALSE)
        }
        x <- stats::setNames(x[, 1L], rownames(x))
    }
    if (!is.numeric(x)) {
        stop("'", arg, "' must be numeric.", call. = FALSE)
    }
    out <- as.numeric(x)
    names(out) <- names(x)
    out
}

.asDonorMatrix <- function(x, arg, prefix) {
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.null(dim(x))) {
        nm <- names(x)
        x <- matrix(x, ncol = 1L, dimnames = list(nm, NULL))
    }
    if (!is.numeric(x)) {
        stop("'", arg, "' must be a numeric matrix of allele dosages ",
             "(donors in rows, SNPs in columns).", call. = FALSE)
    }
    storage.mode(x) <- "double"
    if (is.null(colnames(x))) {
        colnames(x) <- paste0(prefix, seq_len(ncol(x)))
    }
    x
}

.covariateMatrix <- function(x) {
    if (is.matrix(x)) {
        if (!is.numeric(x)) {
            stop("'covariates' must be numeric, or a data frame whose ",
                 "non-numeric columns are expanded to indicators.",
                 call. = FALSE)
        }
        storage.mode(x) <- "double"
        if (is.null(colnames(x))) {
            colnames(x) <- paste0("covariate", seq_len(ncol(x)))
        }
        return(x)
    }
    if (is.numeric(x) && is.null(dim(x))) {
        return(matrix(as.numeric(x), ncol = 1L,
                      dimnames = list(names(x), "covariate1")))
    }
    if (methods::is(x, "DataFrame")) x <- methods::as(x, "data.frame")
    if (!is.data.frame(x)) {
        stop("'covariates' must be a numeric matrix or a data frame.",
             call. = FALSE)
    }
    if (ncol(x) < 1L) {
        stop("'covariates' must contain at least one column.", call. = FALSE)
    }
    if (anyNA(x)) {
        bad <- names(x)[vapply(x, anyNA, logical(1))]
        stop("Covariate(s) ", toString(bad),
             " contain missing values. Impute them or remove the affected ",
             "donors first.", call. = FALSE)
    }
    for (nm in names(x)) {
        if (is.character(x[[nm]]) || is.logical(x[[nm]])) {
            x[[nm]] <- factor(x[[nm]])
        }
        if (is.factor(x[[nm]])) {
            x[[nm]] <- droplevels(factor(x[[nm]], ordered = FALSE))
            if (nlevels(x[[nm]]) < 2L) {
                stop("Covariate '", nm, "' takes a single value across ",
                     "donors and carries no information. Remove it.",
                     call. = FALSE)
            }
        } else if (!is.numeric(x[[nm]])) {
            stop("Covariate '", nm, "' must be numeric, logical, character ",
                 "or a factor.", call. = FALSE)
        }
    }
    rn <- if (.hasRealRownames(x)) rownames(x) else NULL
    M <- stats::model.matrix(~ ., data = x)
    M <- M[, colnames(M) != "(Intercept)", drop = FALSE]
    attr(M, "assign") <- NULL
    attr(M, "contrasts") <- NULL
    rownames(M) <- rn
    storage.mode(M) <- "double"
    M
}

.hasRealRownames <- function(df) {
    # data.frame() stores automatic row names as integers. Only explicit
    # names identify donors.
    rn <- attr(df, "row.names")
    is.character(rn)
}

.makeLabels <- function(labels) {
    out <- stats::setNames(rep(NA_character_, length(.labelNames)),
                           .labelNames)
    if (is.null(labels)) return(out)
    labels <- unlist(labels)
    if (is.null(names(labels)) || any(!names(labels) %in% .labelNames)) {
        stop("'labels' must be named with any of: ",
             toString(.labelNames), ".", call. = FALSE)
    }
    out[names(labels)] <- as.character(labels)
    out
}

# ---- Accessors -------------------------------------------------------------

#' Accessors for MRCCCData objects
#'
#' Extract the components of an [MRCCCData-class] object. All components
#' are aligned by donor.
#'
#' @param x An [MRCCCData-class] object (or, for `tripletLabels()`, an
#'   [MRCCCFit-class] object).
#'
#' @return `donors()` returns the donor identifiers.
#'   `ligandExpression()`, `receptorExpression()` and `pathwayActivity()`
#'   return named numeric vectors. `ligandGenotypes()` and
#'   `receptorGenotypes()` return donors-by-SNPs dosage matrices.
#'   `covariates()` returns the donors-by-covariates matrix.
#'   `tripletLabels()` returns the named character vector of labels.
#'
#' @seealso [MRCCCData()], [buildTriplet()]
#'
#' @examples
#' sim <- simulateMRCCC(n = 120, seed = 1)
#' dat <- buildTriplet(sim$sender, sim$receiver, sim$genotypes,
#'                     ligand = "LIG1", receptor = "REC1",
#'                     pathway = sim$pathwayGenes,
#'                     covariates = c("sex", "age", "PC1", "PC2"))
#' head(donors(dat))
#' colnames(ligandGenotypes(dat))
#' head(covariates(dat))
#' tripletLabels(dat)
#'
#' @name MRCCCData-accessors
#' @rdname MRCCCData-accessors
NULL

#' @rdname MRCCCData-accessors
#' @export
setMethod("donors", "MRCCCData", function(x) names(x@ligandExpression))

#' @rdname MRCCCData-accessors
#' @export
setMethod("ligandExpression", "MRCCCData", function(x) x@ligandExpression)

#' @rdname MRCCCData-accessors
#' @export
setMethod("receptorExpression", "MRCCCData",
          function(x) x@receptorExpression)

#' @rdname MRCCCData-accessors
#' @export
setMethod("pathwayActivity", "MRCCCData", function(x) x@pathwayActivity)

#' @rdname MRCCCData-accessors
#' @export
setMethod("ligandGenotypes", "MRCCCData", function(x) x@ligandGenotypes)

#' @rdname MRCCCData-accessors
#' @export
setMethod("receptorGenotypes", "MRCCCData", function(x) x@receptorGenotypes)

#' @rdname MRCCCData-accessors
#' @export
setMethod("covariates", "MRCCCData", function(x) x@covariates)

#' @rdname MRCCCData-accessors
#' @export
setMethod("tripletLabels", "MRCCCData", function(x) x@labels)

#' @rdname MRCCCData-class
#' @param object An [MRCCCData-class] object.
#' @return `show()` prints a description of the triplet and returns `object`
#'   invisibly.
#' @export
setMethod("show", "MRCCCData", function(object) {
    lab <- object@labels
    showLab <- function(v) if (is.na(v)) "(unnamed)" else v
    cat("MRCCCData: one ligand-receptor-pathway triplet\n")
    cat("  ligand:   ", showLab(lab[["ligand"]]),
        if (!is.na(lab[["sender"]])) {
            paste0("  (sender: ", lab[["sender"]], ")")
        },
        "\n", sep = "")
    cat("  receptor: ", showLab(lab[["receptor"]]),
        if (!is.na(lab[["receiver"]])) {
            paste0("  (receiver: ", lab[["receiver"]], ")")
        },
        "\n", sep = "")
    cat("  pathway:  ", showLab(lab[["pathway"]]), "\n", sep = "")
    cat("  donors:   ", length(object@ligandExpression), "\n", sep = "")
    cat("  instruments: ", ncol(object@ligandGenotypes),
        " SNP(s) for the ligand, ", ncol(object@receptorGenotypes),
        " for the receptor\n", sep = "")
    cv <- colnames(object@covariates)
    cat("  covariates (", length(cv), "): ",
        paste(utils::head(cv, 6L), collapse = ", "),
        if (length(cv) > 6L) ", ...", "\n", sep = "")
    invisible(object)
})
