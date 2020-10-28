#' Load all deg
#'
#' @param id A dataset identifier
#'
#' @return a `list`.
#' @export
#'
#' @examples
#' deg.all <- load_example_deg()
#'
load_example_deg <- function(id = "GSE85871") {
  message("Loading deg ", id, "...")
  if (id == "GSE85871") {

    deg.all <- readRDS('./data-raw/deg.all.RDS')

    message("Done.")

    return(deg.all)
  }
}
