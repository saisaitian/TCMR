devtools::load_all()
library(purrr)

data <- load_example_dataset()

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
  deg_batch_caller(df, batch_groups, ref_group = "DMSO")
})

head(str(reports))
reports2 <- purrr::flatten(reports)

str(reports2)

AnalyzedDEG <- data.frame(
  id = 1:103,
  batch = rep(c("batch1", "batch2", "batch3"), sapply(reports, length)),
  vs = names(reports2),
  filename = paste0("GSE85871-DEG-", 1:103, ".rds")
)

purrr::map2(reports2, AnalyzedDEG$filename, function(DEG, filename) {
  saveRDS(DEG, file = paste0("inst/extdata/", filename))
})

usethis::use_data(AnalyzedDEG, overwrite = TRUE)
