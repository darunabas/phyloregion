#' Progress Bar
#'
#' This function shows the progress bar of analysis.
#'
#' @param x a vector (list) or an expression object.
#' @param FUN the function to be applied to each element of x.
#' @param ... Further arguments passed to or from other methods.
#' @rdname progress
#' @return For lapply, sapply(simplify = FALSE) and replicate(simplify = FALSE), a list.
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @examples
#' \donttest{
#' m <- sapply(1:50000, function(x) list(rnorm(1000)))
#' n <- progress(m, sd)
#' }
#' @export
progress <- function(x, FUN, ...) {
   env <- environment()
   pb_Total <- length(x)
   counter <- 0
   pb <- txtProgressBar(min = 0, max = pb_Total, style = 3,
                        width = getOption("width")/2L)

   # wrapper around FUN
   wrapper <- function(...){
      curVal <- get("counter", envir = env)
      assign("counter", curVal +1 ,envir=env)
      setTxtProgressBar(get("pb", envir=env), curVal +1)
      FUN(...)
   }
   r <- lapply(x, wrapper, ...)
   close(pb)
   r
}
