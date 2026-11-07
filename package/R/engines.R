#' Block engines
#'
#' The per-subset fit is a plugin. An engine is a list with a `name`, a
#' `capabilities` list and a `fit` function with the signature
#' `fit(Y, X, a_forms, control)` returning at least `G`, `E`, and one of `draws`
#' (an `s` by `s` by `D` array) or `se` (an `s` by `s` matrix).
#'
#' LiMMBo delegates its subset fits to limix; indri delegates to stanimal,
#' sommer or ASReml-R, because the right engine depends on the number of
#' individuals, on whether the data have gaps, on whether a licence is
#' available, and on whether the caller wants a posterior or a likelihood.
#'
#' @param name Engine name.
#'
#' @return An engine list.
#' @export
indri_engine <- function(name = "stanimal") {
  if (!is.character(name) || length(name) != 1L || is.na(name))
    stop("`engine` must be a single engine name.", call. = FALSE)
  switch(name,
    stanimal = .indri_engine_stanimal(),
    sommer = .indri_engine_sommer(),
    asreml = .indri_engine_asreml(),
    stop("Unknown engine `", name, "`. Implemented: ",
         paste(indri_engines(), collapse = ", "), ".", call. = FALSE)
  )
}

#' @rdname indri_engine
#' @return `indri_engines()` returns the names of the implemented engines.
#' @export
indri_engines <- function() {
  c("stanimal", "sommer", "asreml")
}
