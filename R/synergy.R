#' Drug Pair Seeker
#'
#' @param data the logFC matrix of TCM
#' @param disease the disease of you intersted
#' @param cutoff DEG cutoff
#'
#' @return the data.frame
#' @export
#'
#' @examples
#' data("data_logFC")
#' data("tcga_disease")
#' tmp <- tcm.SeekDrugPair(data = data_logFC, disease = res_RNA2)
tcm.SeekDrugPair <- function(data = NULL, disease = NULL, cutoff = 1) {
  message(paste0("\n", ">>> Calculating ", "drug pairs to treat the disease"))

  # obtain common genes -----------------------------------------------------

  idx <- intersect(rownames(data), rownames(disease))

  data <- data[idx, ]

  DEG_disease <- disease[idx, , drop = F] %>%
    subset(abs(log2FoldChange) > cutoff)

  # disease_down ------------------------------------------------------------

  down_disease <- DEG_disease %>%
    subset(log2FoldChange < 0) %>%
    rownames()

  # disease_up ------------------------------------------------------------

  up_disease <- DEG_disease %>%
    subset(log2FoldChange > 0) %>%
    rownames()

  # calculate drug pair -----------------------------------------------------

  aa <- data.frame(t(utils::combn(names(data), 2)), check.rows = F, check.names = F, stringsAsFactors = F)

  colnames(aa) <- c("DrugA", "DrugB")

  aa$synergy_score <- 0

  for (i in 1:nrow(aa)) {
    drug_a <- aa[i, 1]

    drug_b <- aa[i, 2]

    drug_a <- data[, drug_a]

    names(drug_a) <- rownames(data)

    down_drug_a <- drug_a[order(-drug_a)] %>%
      tail(500) %>%
      names()

    up_drug_a <- drug_a[order(-drug_a)] %>%
      head(500) %>%
      names()

    drug_b <- data[, drug_b]

    names(drug_b) <- rownames(data)

    down_drug_b <- drug_b[order(-drug_b)] %>%
      tail(500) %>%
      names()

    up_drug_b <- drug_b[order(-drug_b)] %>%
      head(500) %>%
      names()

    synergyscore <- length(intersect(up_disease, down_drug_a)) +
      length(intersect(up_disease, down_drug_b)) +
      length(intersect(down_disease, up_drug_a)) +
      length(intersect(down_disease, up_drug_b)) -
      length(intersect(up_disease, up_drug_a)) -
      length(intersect(up_disease, up_drug_b)) -
      length(intersect(down_disease, down_drug_a)) -
      length(intersect(down_disease, down_drug_b))

    aa[i, 3] <- synergyscore
  }

  message("Done")

  aa <- aa[order(-synergyscore), ]

  return(aa)
}
