#' Differential Expression Analysis
#'
#' - `deg_caller` is used to run two group differential expression analysis.
#' - `deg_batch_caller` is used to run batch DEG analysis.
#'
#' @param data An expression dataset in `data.frame` format whose rows indicate genes and columns indicate samples.
#' @param group A character vector specifying 1 of 2 groups which samples belong to.
#' @param level Levels for group, should be length-2 vector, the first indicates the reference group.
#' At default, it is `NULL`, the group for the first sample will be used as reference group.
#'
#' @return A `data.frame` or a `list` of `data.frame`.
#' @export
#'
#' @examples
#' data <- load_example_dataset()
#'
#' # Get a subset of expression dataset
#' ix <- c(1, 2, 61, 62)
#' expr <- data$expr[, ix]
#' group <- data$pdata$perturbagen[ix]
#'
#' # Run DEG analysis
#' report <- deg_caller(expr, group = group, level = group[c(3, 1)])
#'
#' # Run batch DEG analysis
#' reports <- deg_batch_caller(expr, groups = group, ref_group = group[3])
#'
#' identical(report, reports[[1]])
#' @testexamples
#' expect_is(report, "data.table")
#' expect_is(reports, "list")
deg_caller <- function(data, group, level = NULL) {
  stopifnot(is.data.frame(data), !is.null(rownames(data)))

  if (!is.character(group) | length(unique(group)) != 2) {
    stop("'group' can only be a character vector with two groups!")
  }

  if (is.null(level)) {
    level <- unique(group)
  }

  message("Constructing design matrix...")
  design <- stats::model.matrix(~factor(group, levels = level))

  message("Running DEG analysis with limma...")
  fit <- limma::lmFit(data, design)
  fit <- limma::eBayes(fit)

  message("Reporting results...")
  deg <- limma::topTable(fit, coef = 2, number = Inf)

  message("Done.")
  # The rownames could be any identifiers instead of gene symbols
  return(data.table::as.data.table(deg, keep.rownames = "identifier"))
}

#' @rdname deg_caller
#' @param groups A character vector specifying the group which samples belong to. Should at least 2 groups.
#' @param ref_group A string specifying the reference group.
#' @export
deg_batch_caller <- function(data, groups, ref_group) {
  stopifnot(is.data.frame(data), !is.null(rownames(data)), length(ref_group) == 1L)

  contrast_groups <- setdiff(unique(groups), ref_group)
  reports <- purrr::map(contrast_groups,
                        ~deg_caller(
                          data[, groups %in% c(., ref_group)],
                          group = groups[groups %in% c(., ref_group)],
                          level = c(ref_group, .)))
  names(reports) <- paste(contrast_groups, ref_group, sep = ":")
  reports
}


