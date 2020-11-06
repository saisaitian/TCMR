#' Differential Expression Analysis
#'
#' - `deg_caller` is used to run two group differential expression analysis.
#' - `deg_batch_caller` is used to run batch DEG analysis.
#'
#' @param data An expression dataset in `data.frame` format whose rows indicate genes and columns indicate samples.
#' @param group A character vector specifying 1 of 2 groups which samples belong to.
#' @param ref_group A string specifying the reference group.
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
#' report <- deg_caller(expr, group = group, ref_group = group[3])
#' head(report)
#'
#' # Run batch DEG analysis
#' reports <- deg_batch_caller(expr, groups = group, ref_group = group[3])
#' head(reports)
#'
#' identical(report, reports[[1]])
#' @testexamples
#' expect_is(report, "data.table")
#' expect_is(reports, "list")
#' expect_identical(report, reports[[1]])
deg_caller <- function(data, group, ref_group) {
  stopifnot(is.data.frame(data), !is.null(rownames(data)))

  uniq_group <- unique(group)
  if (!is.character(uniq_group) | length(uniq_group) != 2) {
    stop("'group' can only be a character vector with two groups!")
  }

  level <- c(ref_group, setdiff(uniq_group, ref_group))

  stable <- table(group)
  message("Info: ", level[2], " vs ", level[1], " (reference group)")
  message(paste0("N: ", paste0(paste(names(stable), stable, sep = ":#"), collapse = "  ")))

  message("Constructing design matrix...")
  design <- stats::model.matrix(~ factor(group, levels = level))

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
#' @export
deg_batch_caller <- function(data, groups, ref_group) {
  stopifnot(is.data.frame(data), !is.null(rownames(data)), length(ref_group) == 1L)

  contrast_groups <- setdiff(unique(groups), ref_group)
  reports <- purrr::map(
    contrast_groups,
    ~ deg_caller(
      data[, groups %in% c(., ref_group)],
      group = groups[groups %in% c(., ref_group)],
      ref_group = ref_group
    )
  )
  names(reports) <- paste(contrast_groups, ref_group, sep = ":")
  reports
}
