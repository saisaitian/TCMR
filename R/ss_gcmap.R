#' gcmap Method for Signature Search
#'
#' @inherit ss_cor
#' @inheritParams ss_cmap
#' @param higher A cutoff value. If `logFC` larger than or equal to `higher` will be included
#' in the gene set with  `+1`.
#' @param lower A cutoff value. If `logFC` smaller than or equal to `higher` will be included
#' in the gene set to `-1`, others should be set to `0`.
#'
#' @export
#' @importFrom ExperimentHub ExperimentHub
#' @examples
#' data("data_logFC")
#' query2 <- data_logFC[1:60, 1, drop = FALSE]
#' gcmap_kk <- ss_gcmap(input = query2, data = data_logFC)
ss_gcmap <- function(input, data, higher = 1, lower = -1, cores = 1L) {
  stopifnot(is.data.frame(input))
  num <- sum(rownames(input) %in% rownames(data))
  if (num <= 10) {
    stop("the commom gene less than 10")
  }
  message(paste(
    num, "/", length(rownames(data)),
    "genes in input share identifiers with reference database"
  ))
  data2 <- ifelse(data > higher, 1, ifelse(data > -lower, 0, -1))
  ## subset objects to shared genes
  matched.features <- match(rownames(input), rownames(data2))
  matched.sets <- data2[stats::na.omit(matched.features), ]

  ## extract scores for each gene set
  sets.up <- parallel::mclapply(
    seq(ncol(matched.sets)),
    function(x) which(matched.sets[, x] == 1),
    mc.cores = cores
  )

  sets.down <- parallel::mclapply(
    seq(ncol(matched.sets)),
    function(x) which(matched.sets[, x] == -1),
    mc.cores = cores
  )

  ## transform experiment to (reverse) ranks
  rank.matrix <- apply(input, 2, function(x) {
    length(x) - rank(x) + 1
  })

  ## calculate connectivity score
  raw.score <- apply(rank.matrix, 2, function(r) {
    vapply(seq_along(sets.up), function(n) {
      .s(r[sets.up[[n]]], r[sets.down[[n]]], length(r))
    }, FUN.VALUE = numeric(1))
  })

  raw.score <- matrix(raw.score, ncol = ncol(input))

  eh <- suppressMessages(ExperimentHub::ExperimentHub())
  CSnull <- suppressMessages(eh[["EH3234"]])
  CSnull[CSnull[, "Freq"] == 0, "Freq"] <- 1
  myrounding <- max(nchar(as.character(CSnull[, "WTCS"]))) - 3
  gcmap_round <- round(as.numeric(raw.score), myrounding)

  CS_pval <- vapply(gcmap_round, function(x) {
    sum(CSnull[abs(CSnull[, "WTCS"]) > abs(x), "Freq"]) / sum(CSnull[, "Freq"])
  }, FUN.VALUE = numeric(1))
  CS_fdr <- stats::p.adjust(CS_pval, "fdr")

  score <- matrix(raw.score, ncol = ncol(input))
  ## store results
  results <- data.frame(
    tcm = colnames(data2),
    direction = ifelse(score[, 1] >= 0, "up", "down"),
    raw_score = score[, 1],
    Pval = CS_pval,
    FDR = CS_fdr,
    Nset = colSums(as.matrix(abs(matched.sets)))
  )


  ## Apply scaling of scores to full data set
  results[, "scaled_score"] <- .S(results$raw_score)

  results <- results[order(abs(results$scaled_score), decreasing = TRUE), ]
  results <- results[, c(
    "tcm", "direction", "raw_score",
    "scaled_score", "Pval", "FDR", "Nset"
  )]

  results <- results[order(abs(results$scaled_score), decreasing = TRUE), ]
  names(results) <- c("tcm", "direction", "raw_score", "scaled_score", "pvalue", "fdr", "Nset")

  return(results)
}
