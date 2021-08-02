#' Load Analyzed DEG Results
#'
#' Currently, the DEG results are for GSE85871.
#'
#' @param id A vector for subset of `id` column in [AnalyzedDEG] or just a subset of [AnalyzedDEG].
#'
#' @return A `data.frame` or a list of `data.frame`.
#' @export
#'
#' @examples
#'
#' # Firstly, check the dataset
#' data("AnalyzedDEG")
#' head(AnalyzedDEG)
#'
#' # Pick up what you want to load
#' # A subset of AnalyzedDEG is recommended,
#' # you can use %>%
#' # Otherwise specify a subset of id column
#'
#' # e.g. head 5 sets
#' head5_reports <- AnalyzedDEG %>%
#'   head(n = 5) %>%
#'   load_analyzedDEG()
#'
#' # only the second
#' one_report <- load_analyzedDEG(2)
#' @testexamples
#' expect_is(head5_reports, "list")
#' expect_is(one_report, "data.frame")
load_analyzedDEG <- function(id) {
  if (is.data.frame(id)) {
    filename <- id$filename
    if (is.null(filename)) {
      stop("When input a data.frame, should be a subset of AnalyzedDEG!")
    }
  } else if (is.numeric(id)) {
    #AnalyzedDEG <- get("AnalyzedDEG", envir = as.environment("package:TCMR"))
    AnalyzedDEG <- get("AnalyzedDEG")
    filename <- AnalyzedDEG[AnalyzedDEG$id %in% as.integer(id), ]$filename
  } else {
    stop("Invalid input for 'id', check the documentation!")
  }

  res <- list()
  for (i in seq_along(filename)) {
    res[[i]] <- readRDS(system.file("extdata", filename[i], package = "TCMR", mustWork = TRUE))
  }

  if (length(res) == 1L) {
    res[[1]]
  } else {
    res
  }
}
