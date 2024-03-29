library(purrr)
data("AnalyzedDEG")
one_report <- tcm.LoadAnalyzedDEG(2)
# DMSO are the controls

degreports <- purrr::map(1:103, function(i) {
  one_report <- tcm.LoadAnalyzedDEG(i)
  deg <- subset(one_report, abs(logFC) > 1 & P.Value < 0.05)
  aa <- getsigpath(deg, p.cutoff = 0.05, p.adj.cutoff = 1, n.path = 10)
  return(aa)
})


Analyzedsigpath <- data.frame(
  id = 1:103,
  batch = rep(c("batch1", "batch2", "batch3"), sapply(reports, length)),
  vs = names(reports2),
  filename = paste0("GSE85871-sigpath-", 1:103, ".rds")
)

purrr::map2(degreports, Analyzedsigpath$filename, function(DEG, filename) {
  saveRDS(DEG, file = paste0("inst/extdata/", filename))
})

AnalyzedSigPathway <- Analyzedsigpath

AnalyzedSigPathway$vs <- gsub("尾-", "", AnalyzedSigPathway$vs)

AnalyzedSigPathway <- AnalyzedSigPathway[1:102, ]

usethis::use_data(AnalyzedSigPathway, overwrite = TRUE)
