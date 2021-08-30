.zenodo_dir <- "https://zenodo.org/record/5336709/files"
.zenodo_local <- getOption("zenodo", default = path.expand("~/.tcmR"))

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
    expr <- data.table::fread(tcm.LoadData("GSE85871_expr.tsv.gz", just_query = TRUE), data.table = FALSE) %>%
      tibble::column_to_rownames("symbol")
    pdata <- data.table::fread(tcm.LoadData("GSE85871_pdata.tsv.gz", just_query = TRUE), data.table = FALSE)
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
    res[[i]] <- readRDS(tcm.LoadData(filename[i], just_query = TRUE))
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


# Utils to get data from remote -------------------------------------------

#' Load builtin/remote dataset
#'
#' @param name a dataset name, e.g., "brca_disease", "TCM_expr", "tfdata".
#' Check <https://zenodo.org/record/5336709> for remote data.
#' @param just_query if `TRUE`, (download and) load data. Only set `TRUE` for
#' data in `rda` format.
#'
#' @return a data or a string represting data path.
#' @export
#'
#' @examples
#' tcm.LoadData("SMILES")
tcm.LoadData <- function(name, just_query = FALSE) {
  stopifnot(length(name) == 1)
  name2 <- if (just_query) name else paste0(name, ".rda")
  data_path <- file.path(.zenodo_local, name2)
  if (!dir.exists(dirname(data_path))) dir.create(dirname(data_path), showWarnings = FALSE, recursive = TRUE)

  # Keep in remote "brca_disease", "TCM_expr", "tfdata"

  # builtin datasets
  available_datasets <- c(
    "TCM_pdata",
    "SMILES",
    "signature_tme",
    "data_logFC",
    "apfp",
    "AnalyzedSigPathway",
    "AnalyzedDEG"
  )
  if (name %in% available_datasets) {
    # The data is builtin
    data(list = name, package = "UCSCXenaShiny", envir = environment())
  } else {
    if (!file.exists(data_path)) {
      data_url <- file.path(.zenodo_dir, name2)
      message("Loading data from remote: ", data_url, ", please wait...")
      name <- FALSE
      tryCatch(
        {
          download.file(data_url, data_path)
          message("Data has been saved to ", data_path)
        },
        error = function(e) {
          message("Data load failed, please check your input and the internet.\n NULL will be returned.")
          if (file.exists(data_path)) unlink(data_path, recursive = TRUE, force = TRUE)
          name <<- TRUE
        }
      )
      if (isTRUE(name)) {
        return(invisible(NULL))
      }
    }
    if (!just_query) {
      tryCatch(
        load(data_path, envir = environment()),
        error = function(e) {
          message("Data load failed, probably due to broken download file, please try again.\n This time NULL will be returned.")
          if (file.exists(data_path)) unlink(data_path, recursive = TRUE, force = TRUE)
          name <<- TRUE
        }
      )
    }
    if (isTRUE(name)) {
      return(invisible(NULL))
    }
  }

  if (!just_query) {
    return(get(setdiff(ls(), c("name2", "name", "data_path", "data_url", "available_datasets", "just_query"))))
  } else {
    return(data_path)
  }
}

# Substitute func names
# sed -i "" "s/tcm.MatrixDotplot/tcm.MatrixDotplot/g" `grep "tcm.MatrixDotplot" -rl ./*`
