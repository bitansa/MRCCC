# S4 generics of MRCCC. Methods are defined next to the classes they act on.

#' @rdname mrccc
#' @export
setGeneric("mrccc", function(x, ...) standardGeneric("mrccc"))

# ---- MRCCCData accessors ---------------------------------------------------

#' @rdname MRCCCData-accessors
#' @export
setGeneric("donors", function(x) standardGeneric("donors"))

#' @rdname MRCCCData-accessors
#' @export
setGeneric("ligandExpression", function(x) standardGeneric("ligandExpression"))

#' @rdname MRCCCData-accessors
#' @export
setGeneric("receptorExpression",
           function(x) standardGeneric("receptorExpression"))

#' @rdname MRCCCData-accessors
#' @export
setGeneric("pathwayActivity", function(x) standardGeneric("pathwayActivity"))

#' @rdname MRCCCData-accessors
#' @export
setGeneric("ligandGenotypes", function(x) standardGeneric("ligandGenotypes"))

#' @rdname MRCCCData-accessors
#' @export
setGeneric("receptorGenotypes",
           function(x) standardGeneric("receptorGenotypes"))

#' @rdname MRCCCData-accessors
#' @export
setGeneric("covariates", function(x) standardGeneric("covariates"))

#' @rdname MRCCCData-accessors
#' @export
setGeneric("tripletLabels", function(x) standardGeneric("tripletLabels"))

# ---- MRCCCFit accessors ----------------------------------------------------

#' @rdname MRCCCFit-accessors
#' @export
setGeneric("pip", function(x) standardGeneric("pip"))

#' @rdname MRCCCFit-accessors
#' @export
setGeneric("pipMCSE", function(x) standardGeneric("pipMCSE"))

#' @rdname MRCCCFit-accessors
#' @export
setGeneric("effectSizes", function(x, ...) standardGeneric("effectSizes"))

#' @rdname MRCCCFit-accessors
#' @export
setGeneric("convergence", function(x) standardGeneric("convergence"))

#' @rdname MRCCCFit-accessors
#' @export
setGeneric("instrumentStrength",
           function(x) standardGeneric("instrumentStrength"))

#' @rdname MRCCCFit-accessors
#' @export
setGeneric("posteriorDraws", function(x, ...) standardGeneric("posteriorDraws"))

#' @rdname MRCCCFit-accessors
#' @export
setGeneric("firstStageEffects",
           function(x) standardGeneric("firstStageEffects"))

#' @rdname MRCCCFit-accessors
#' @export
setGeneric("posteriorMeans", function(x) standardGeneric("posteriorMeans"))

#' @rdname MRCCCFit-accessors
#' @export
setGeneric("inputScale", function(x) standardGeneric("inputScale"))

#' @rdname MRCCCFit-accessors
#' @export
setGeneric("fitSettings", function(x) standardGeneric("fitSettings"))

# ---- Posterior summaries ---------------------------------------------------

#' @rdname credibleIntervals
#' @export
setGeneric("credibleIntervals",
           function(x, ...) standardGeneric("credibleIntervals"))

#' @rdname signProbabilities
#' @export
setGeneric("signProbabilities",
           function(x, ...) standardGeneric("signProbabilities"))

#' @rdname signReversal
#' @export
setGeneric("signReversal", function(x, ...) standardGeneric("signReversal"))

#' @rdname plotEffectCurve
#' @export
setGeneric("plotEffectCurve",
           function(x, ...) standardGeneric("plotEffectCurve"))
