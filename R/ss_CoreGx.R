#' CoreGx method for Signature Search
#'
#' @inherit ss_cor
#' @export
#'
#' @examples
#' data("data_logFC")
#' query2 <- data_logFC[1:60, 1, drop = FALSE]
#' txp <- ss_CoreGx(query2, data_logFC[1:3])
#' @testexamples
#' expect_is(txp, "data.frame")

ss_CoreGx <- function(input, data) {
  stopifnot(is.data.frame(input))
  num <- sum(rownames(input) %in% rownames(data))
  if (num <= 10) {
    stop("the commom gene less than 10")
  }
  message(paste(
    num, "/", length(rownames(data)),
    "genes in input share identifiers with reference database"
  ))

  res_list <- NULL
  common <- intersect(rownames(input), rownames(data))
  input2 <- input[common, , drop = F]

  res_list <- parallel::mclapply(
    seq_len(ncol(data)),
    function(i) {
      tt <- suppressWarnings(CoreGx::connectivityScore(data[, i, drop = F], input2, method = "fgsea", nperm = 1e4, nthread = 1))
      direction <- tt[1]
      direction[direction >= 0] <- "up"
      direction[direction < 0] <- "down"
      res <- data.frame(
        tcm = colnames(data)[i],
        direction = direction,
        raw_score = tt[1],
        p = tt[2],
        Nset = length(common),
        stringsAsFactors = FALSE
      )
    },
    mc.cores = parallel::detectCores()
  )

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
