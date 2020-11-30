#' Load Analyzed Pathways Results
#'
#' Currently, the top ten pathways results are for GSE85871.
#'
#' @param id A vector for subset of `id` column in [Analyzedsigpath] or just a subset of [Analyzedsigpath].
#'
#' @return A `data.frame` or a list of `data.frame`.
#' @export
#'
#' @examples
#'
#' # Firstly, check the dataset
#' data("Analyzedsigpath")
#' head(Analyzedsigpath)
#'
#' # Pick up what you want to load
#' # A subset of Analyzedsigpath is recommended,
#' # you can use %>%
#' # Otherwise specify a subset of id column
#'
#' # e.g. head 5 sets
#' head5_reports <- Analyzedsigpath %>%
#'   head(n = 5) %>%
#'   load_analyzedsigpath()
#'
#' # only the second
#' one_report <- load_analyzedsigpath(2)
#'


load_analyzedsigpath <- function(id) {
  if (is.data.frame(id)) {
    filename <- id$filename
    if (is.null(filename)) {
      stop("When input a data.frame, should be a subset of Analyzedsigpath!")
    }
  } else if (is.numeric(id)) {
    Analyzedsigpath <- get("Analyzedsigpath", envir = as.environment("package:TCMR"))
    filename <- Analyzedsigpath[Analyzedsigpath$id %in% as.integer(id), ]$filename
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






