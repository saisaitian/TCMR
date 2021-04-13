#' Get Enriched Signal Pathways
#'
#' @param data A `data.frame` from a subset of [deg_caller] result.
#' @param p.cutoff the p value cutoff.
#' @param p.adj.cutoff the cutoff value of p.adjust.
#' @param n.path A integer value to indicate how many top pathways sorted by p value should be identified; 10 by default.
#' @param sortby A string value to indicate the pathway was sorted by p value or p.adjust value. Allowed value contains c('pvalue', 'p.adjust').
#' @importFrom clusterProfiler bitr
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom clusterProfiler enrichGO
#' @importFrom clusterProfiler setReadable
#' @importFrom DOSE parse_ratio
#' @return  A `data.frame`.
#'
#' @export
#'
#' @examples
#' data("AnalyzedDEG")
#' one_report <- load_analyzedDEG(1)
#' deg <- subset(one_report, abs(logFC) > 1 & P.Value < 0.05)
#' path1 <- get_sig_pathway(deg, p.cutoff = 0.05, p.adj.cutoff = 0.5, n.path = 10)
#' path2 <- get_sig_pathway(deg, p.cutoff = 0.05, p.adj.cutoff = 0.5, sortby = "p.adjust", n.path = 10)
#' @testexamples
#' expect_is(deg, "data.frame")
#' expect_is(path1, "data.frame")
#' expect_is(path2, "data.frame")

get_sig_pathway <- function(data,
                            p.cutoff = 0.05,
                            p.adj.cutoff = 0.2,
                            sortby = "pvalue",
                            n.path = 10) {
  stopifnot(is.data.frame(data), !is.null(rownames(data)))
  identifiers <- suppressWarnings(clusterProfiler::bitr(data$identifier, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db"))
  message("The number of gene is: ", paste0(nrow(identifiers), collapse = "  "), " for next analysis")
  kk <- suppressWarnings(clusterProfiler::enrichKEGG(
    gene = identifiers$ENTREZID,
    organism = "hsa",
    pvalueCutoff = p.cutoff,
    qvalueCutoff = p.adj.cutoff
  ))

  df <- as.data.frame(
    setReadable(kk, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "ENTREZID")
  )

  df <- df[order(df[, sortby], df$GeneRatio), ]
  df$Description <- factor(df$Description, levels = df$Description)
  df$GeneRatio <- DOSE::parse_ratio(df$GeneRatio)

  if (nrow(df) > n.path) {
    df <- df[1:n.path, ]
  }

  return(df)
}
