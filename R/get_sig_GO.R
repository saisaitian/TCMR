
#' Get Enriched GO Sets
#'
#' @inheritParams get_sig_pathway
#' @inheritParams clusterProfiler::enrichGO
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' data("AnalyzedDEG")
#' one_report <- load_analyzedDEG(1)
#' deg <- subset(one_report, abs(logFC) > 1 & P.Value < 0.05)
#'
#' go <- get_sig_GO(deg, p.cutoff = 0.05, p.adj.cutoff = 0.5, n.path = 10)
#' go
#' dot <- dotplot(go)
#' dot
#' bar <- barplot(go)
#' bar
#' @testexamples
#' expect_is(go, "data.frame")
#' expect_is(dot, "ggplot")
#' expect_is(bar, "ggplot")
get_sig_GO <- function(data,
                       p.cutoff = 0.05,
                       p.adj.cutoff = NULL,
                       sortby = "pvalue",
                       ont = "BP",
                       n.path = NULL) {
  stopifnot(is.data.frame(data), !is.null(rownames(data)))
  identifiers <- suppressWarnings(clusterProfiler::bitr(data$identifier, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db"))
  message("The number of gene is: ", paste0(nrow(identifiers), collapse = "  "), " for next analysis")
  kk <- suppressWarnings(clusterProfiler::enrichGO(
    gene = identifiers$ENTREZID,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = ont,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  ))
  df <- as.data.frame(kk)
 if(is.numeric(p.cutoff)&&is.numeric(p.adj.cutoff)){
  df <- df[df$pvalue<p.cutoff,]
  df <- df[df$p.adjust<p.adj.cutoff,]
 }

  if(is.numeric(p.cutoff)&&is.null(p.adj.cutoff)){
    df <- df[df$pvalue<p.cutoff,]
  }
  if(is.null(p.cutoff)&&is.numeric(p.adj.cutoff)){
    df <- df[df$p.adjust<p.adj.cutoff,]
  }
  df <- df[order(df[, sortby], df$GeneRatio), ]
  df$Description <- factor(df$Description, levels = df$Description)
  df$GeneRatio <- DOSE::parse_ratio(df$GeneRatio)

if(is.numeric(n.path)){

  if (nrow(df) > n.path) {
    df <- df[1:n.path, ]
  }

}
  return(df)
}
