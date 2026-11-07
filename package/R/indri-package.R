#' indri: scalable additive genetic covariance estimation
#'
#' Fits a multivariate animal model on many small trait subsets and reassembles
#' the subset estimates into one positive semi-definite additive genetic
#' covariance matrix. See `NOTES.md` in the repository root for the estimator,
#' the engine contract and the milestones.
#'
#' @references
#' Meyer, H. V., Casale, F. P., Stegle, O., and Birney, E. (2018). LiMMBo:
#' a simple, scalable approach for linear mixed models in high-dimensional
#' genetic association studies. bioRxiv, preprint.
#' \doi{10.1101/255497}
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats optim sd
#' @importFrom utils modifyList
NULL
