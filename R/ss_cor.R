#' Cor method for for Signature Search
#'
#' @param input A `data.frame` containing `logFC` value, which `rownames` should be gene symbols.
#' @param data A `data.frame` contains `logFC` value data from 103 compounds.
#' @param method A string specifying correlation calculation method, should be one of "pearson", "spearman", "kendall".
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' data("data_logFC")
#' query2 <- data_logFC[1:60, 1, drop = FALSE]
#' aa <- ss_cor(query2, data_logFC, method = "pearson")
#' aaa <- ss_cor(query2, data_logFC, method = "spearman")
#' aaaa <- ss_cor(query2, data_logFC, method = "kendall")
ss_cor <- function(input, data, method = c("pearson", "spearman", "kendall")) {
  stopifnot(is.data.frame(input))
  method <- match.arg(method)
  num <- sum(rownames(input) %in% rownames(data))
  if (num <= 10) {
    stop("the commom gene less than 10")
  }
  message(paste(
    num, "/", length(rownames(data)),
    "genes in input share identifiers with reference database"
  ))
  if (is.null(input)) {
    stop(" Input is NULL !")
  }

  res_list <- NULL
  common <- intersect(rownames(input), rownames(data))
  input2 <- input[common, , drop = F]
  data2 <- data[common, ]

  for (i in 1:ncol(data2)) {
    tt <- stats::cor.test(input2[, 1], data2[, i], method = method)
    p <- tt$p.value
    cor <- tt$estimate
    direction <- cor
    direction[cor >= 0] <- "up"
    direction[cor < 0] <- "down"
    res <- data.frame(
      tcm = colnames(data2)[i],
      direction = direction,
      raw_score = cor,
      p = p,
      Nset = length(common),
      stringsAsFactors = FALSE
    )
    res_list <- c(res_list, list(res))
  }
  result <- do.call(rbind, res_list)
  result$fdr <- stats::p.adjust(result$p, "fdr")
  result$scaled_score <- .S(result$raw_score)

  result <- result[order(abs(result$scaled_score), decreasing = TRUE), ]
  result <- result[, c(
    "tcm", "direction", "raw_score",
    "scaled_score", "p", "fdr", "Nset"
  )]
  names(result) <- c("tcm", "direction", "raw_score", "scaled_score", "pvalue", "fdr", "Nset")
  rownames(result) <- NULL
  return(result)
}
