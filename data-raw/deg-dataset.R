devtools::load_all()
library(purrr)

data <- tcm.LoadExampleDataset()

groups <- data$pdata$perturbagen

# DMSO are the controls

batch_list <- list(
  batch1 = 1:62,
  batch2 = 63:130,
  batch3 = 131:212
)

reports <- purrr::map(batch_list, function(b) {
  batch_groups <- groups[b]
  df <- data$expr[, b]
  tcm.RunBatchDEG(df, batch_groups, ref_group = "DMSO")
})

head(str(reports))
reports2 <- purrr::flatten(reports)

AnalyzedDEG <- data.frame(
  id = 1:103,
  batch = rep(c("batch1", "batch2", "batch3"), sapply(reports, length)),
  vs = names(reports2),
  filename = paste0("GSE85871-DEG-", 1:103, ".rds")
)


purrr::map2(reports2, AnalyzedDEG$filename, function(DEG, filename) {
  saveRDS(DEG, file = paste0("inst/extdata/", filename))
})

data <- reports2[[1]][, c(1, 2)]

for (i in 2:103) {
  tmp <- reports2[[i]][, 2]
  data <- cbind.data.frame(data, tmp)
}

names(data)[2:104] <- gsub(":DMSO", "", as.character(AnalyzedDEG$vs))

str(data)
data <- data.frame(data, stringsAsFactors = F)

rownames(data) <- data$identifier

data <- data[, -1]
data_logFC <- data

names(data_logFC)[c(76, 91, 96)] <- c("beta_ecdysterone", "6_gingerol", "(+)2-(1-hydroxyl-4-oxocyclohexyl) ethyl caffeate")

AnalyzedDEG$vs[c(76, 91, 96)] <- c("beta_ecdysterone", "6_gingerol", "(+)2-(1-hydroxyl-4-oxocyclohexyl) ethyl caffeate")

AnalyzedDEG <- AnalyzedDEG[1:102, ]

data_logFC <- data_logFC[, 1:102]

usethis::use_data(data_logFC, overwrite = TRUE)

usethis::use_data(AnalyzedDEG, overwrite = TRUE)
