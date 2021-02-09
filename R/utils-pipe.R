#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

set_cores <- function() {
  if (xfun::is_windows()) {
    1
  } else {
    as.integer(parallel::detectCores())
  }
}