library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(xlsx)
library(tidyverse)
options(stringsAsFactors = F)
Sys.setenv(LANGUAGE = "en")
# 感兴趣的signature
gmtfile <- "G:/小丫/FigureYa136fgsea/h.all.v7.0.symbols.gmt"
msigdb.path <- system.file("extdata", "c2.cp.kegg.v7.2.symbols.gmt",
  package = "TCMR", mustWork = TRUE
)

hallmark <- read.gmt(gmtfile)

hallmark <- read.gmt(msigdb.path)
# hallmark.list <- hallmark %>%
#   split(.$term) %>%
#   lapply("[[", 2)
data("data_logFC")

geneList <- data_logFC[, 1]
names(geneList) <- rownames(data_logFC)
geneList <- sort(geneList, decreasing = TRUE) # ranked gene set

clustergsea <- suppressWarnings(clusterProfiler::GSEA(
  geneList = geneList, TERM2GENE = hallmark, verbose = F,
  minGSSize = 15, maxGSSize = 500, nPerm = 10000, seed = TRUE, pvalueCutoff = 1
))


gsea.list <- suppressWarnings(clusterProfiler::GSEA(
  geneList = geneList,
  TERM2GENE = msigdb,
  nPerm = nPerm,
  minGSSize = minGSSize,
  maxGSSize = maxGSSize,
  seed = TRUE,
  verbose = FALSE,
  pvalueCutoff = 1
))

















kkk <- clustergsea@result

kkk2 <- kkk[, c(1, 6, 7, 4, 5)]
names(kkk2) <- c("pathway", "pval", "padj", "ES", "NES")

kkk2 <- data.table::data.table(kkk2)

topPathwaysUp <- kkk2[ES > 0][head(order(pval), n = 10), pathway]

topPathwaysDown <- kkk2[ES < 0][head(order(pval), n = 10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
length(topPathways)

plotGseaTable(hallmark.list[topPathways], logfc, kkk2,
  gseaParam = 0.5
)
