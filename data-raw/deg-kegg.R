
devtools::load_all()
library(purrr)
library(clusterProfiler)
library(ggplot2)
data("AnalyzedDEG")
one_report <- tcm.LoadAnalyzedDEG(1)
names(one_report)

deg <- subset(one_report, abs(logFC) > 1 & P.Value < 0.05)

getsigpath <- function(data,
                       p.cutoff = 0.05,
                       p.adj.cutoff = NULL,
                       sortby = "pvalue",
                       n.path = 10) {
  stopifnot(is.data.frame(data), !is.null(rownames(data)))
  identifiers <- suppressWarnings(bitr(data$identifier, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db"))
  message("The number of gene is: ", paste0(nrow(identifiers), collapse = "  "), " for next analysis")
  kk <- suppressWarnings(enrichKEGG(
    gene = identifiers$ENTREZID,
    organism = "hsa",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  ))

  kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

  df <- as.data.frame(kk)

  if (!is.null(p.cutoff)) {
    df <- df[which(df$pvalue < p.cutoff), ]
  }

  if (!is.null(p.adj.cutoff)) {
    df <- df[which(df$p.adjust < p.adj.cutoff), ]
  }

  if (!is.null(p.cutoff) & !is.null(p.adj.cutoff)) {
    df <- df[which(df$pvalue < p.cutoff & df$p.adjust < p.adj.cutoff), ]
  }

  if (nrow(df) > n.path) {
    pathway <- df[1:n.path, ]
  } else {
    pathway <- df
  }

  sortdf <- pathway[order(pathway[, sortby], pathway$GeneRatio), ]

  sortdf$Description <- factor(sortdf$Description, levels = sortdf$Description)

  sortdf$GeneRatio <- DOSE::parse_ratio(sortdf$GeneRatio)

  return(sortdf)
}

aa <- getsigpath(deg, p.cutoff = 0.05, p.adj.cutoff = 0.5, n.path = 10)
aaa <- getsigpath(deg, p.cutoff = 0.05, p.adj.cutoff = 0.5, sortby = "p.adjust", n.path = 10)
