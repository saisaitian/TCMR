#' Load Example Dataset for Analysis
#'
#' @param id A dataset identifier.
#'
#' @return a `list`.
#' @export
#'
#' @examples
#' x <- tcm.LoadExampleDataset()
#' head(x$expr)
#' head(x$pdata)
#' @testexamples
#' expect_is(x, "list")
tcm.LoadExampleDataset <- function(id = "GSE85871") {
  message("Loading dataset ", id, "...")
  if (id == "GSE85871") {
    expr <- data.table::fread(system.file(
      "extdata", "GSE85871_expr.tsv.gz",
      package = "TCMR", mustWork = TRUE
    ), data.table = FALSE) %>%
      tibble::column_to_rownames("symbol")
    pdata <- data.table::fread(system.file(
      "extdata", "GSE85871_pdata.tsv.gz",
      package = "TCMR", mustWork = TRUE
    ), data.table = FALSE)
    message("Done.")
    return(list(
      expr = expr,
      pdata = pdata
    ))
  }
}


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
#'   tcm.LoadAnalyzedDEG()
#'
#' # only the second
#' one_report <- tcm.LoadAnalyzedDEG(2)
#' @testexamples
#' expect_is(head5_reports, "list")
#' expect_is(one_report, "data.frame")
tcm.LoadAnalyzedDEG <- function(id) {
  if (is.data.frame(id)) {
    filename <- id$filename
    if (is.null(filename)) {
      stop("When input a data.frame, should be a subset of AnalyzedDEG!")
    }
  } else if (is.numeric(id)) {
    # AnalyzedDEG <- get("AnalyzedDEG", envir = as.environment("package:TCMR"))
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


#' Load Analyzed Pathways Results
#'
#' Currently, the top ten pathways results are for GSE85871.
#'
#' @param id A vector for subset of `id` column in [AnalyzedSigPathway] or just a subset of [AnalyzedSigPathway].
#'
#' @return A `data.frame` or a list of `data.frame`.
#' @export
#'
#' @examples
#'
#' # Firstly, check the dataset
#' data("AnalyzedSigPathway")
#' head(AnalyzedSigPathway)
#'
#' # Pick up what you want to load
#' # A subset of AnalyzedSigPathway is recommended,
#' # you can use %>%
#' # Otherwise specify a subset of id column
#'
#' # e.g. head 5 sets
#' head5_reports <- AnalyzedSigPathway %>%
#'   head(n = 5) %>%
#'   tcm.LoadAnalyzedSigPathway()
#'
#' # only the second
#' one_report <- tcm.LoadAnalyzedSigPathway(2)
tcm.LoadAnalyzedSigPathway <- function(id) {
  if (is.data.frame(id)) {
    filename <- id$filename
    if (is.null(filename)) {
      stop("When input a data.frame, should be a subset of AnalyzedSigPathway!")
    }
  } else if (is.numeric(id)) {
    AnalyzedSigPathway <- get("AnalyzedSigPathway", envir = as.environment("package:TCMR"))
    filename <- AnalyzedSigPathway[AnalyzedSigPathway$id %in% as.integer(id), ]$filename
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

# Substitute func names
# sed -i "" "s/tcm.MatrixDotplot/tcm.MatrixDotplot/g" `grep "tcm.MatrixDotplot" -rl ./*`
