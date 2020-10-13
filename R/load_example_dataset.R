#' Load Example Dataset for Analysis
#'
#' @param id A dataset identifier.
#'
#' @return a `list`.
#' @export
#'
#' @examples
#' x <- load_example_dataset()
#' head(x$expr)
#' head(x$pdata)
#'
#' @testexamples
#' expect_is(x, "list")
load_example_dataset <- function(id = "GSE85871") {
  message("Loading dataset ", id, "...")
  if (id == "GSE85871") {

    expr <- data.table::fread(system.file(
      "extdata", "GSE85871_expr.tsv.gz", package = "TCMR", mustWork = TRUE
    ), data.table = FALSE) %>%
      tibble::column_to_rownames("symbol")
    pdata <- data.table::fread(system.file(
      "extdata", "GSE85871_pdata.tsv.gz", package = "TCMR", mustWork = TRUE
    ), data.table = FALSE)
    message("Done.")
    return(list(
      expr = expr,
      pdata = pdata
    ))
  }
}
