
#' Title
#'
#' @param data
#' @param p.cutoff
#' @param p.adj.cutoff
#' @param sortby
#' @param n.path
#'
#' @return
#' @export
#'
#' @examples
#' data("AnalyzedDEG")
#' one_report <- load_analyzedDEG(1)
#' deg <- subset(one_report,abs(logFC)>1&P.Value<0.05)

#' aa=getsigGO(deg,p.cutoff = 0.05,p.adj.cutoff = 0.5,n.path = 10)

#' getdotplot(aa)

#' getbarplot(aa)

getsigGO <- function(data,
                       p.cutoff = 0.05,
                       p.adj.cutoff = NULL,
                       sortby = 'pvalue',
                       ont = 'BP',
                       n.path = 10){

  stopifnot(is.data.frame(data), !is.null(rownames(data)))
  identifiers = suppressWarnings(clusterProfiler::bitr(data$identifier, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))
  message("The number of gene is: ",paste0(nrow(identifiers), collapse = "  "),' for next analysis')
  kk <- suppressWarnings(enrichGO(gene= identifiers$ENTREZID,
                                    OrgDb = org.Hs.eg.db,
                                    ont = ont,
                                    pvalueCutoff = 1,
                                    qvalueCutoff = 1,
                                    readable = TRUE))

  df <- as.data.frame(kk)

  if(!is.null(p.cutoff) ){
    df <- df[which(df$pvalue < p.cutoff),]
  }

  if(!is.null(p.adj.cutoff) ){
    df <- df[which(df$p.adjust < p.adj.cutoff),]
  }

  if(!is.null(p.cutoff)& !is.null(p.adj.cutoff)){
    df <- df[which(df$pvalue < p.cutoff & df$p.adjust < p.adj.cutoff),]
  }

  if(nrow(df) > n.path) {
    pathway<- df[1:n.path,]
  } else {
    pathway <- df
  }

  sortdf <- pathway[order(pathway[,sortby],pathway$GeneRatio),]

  sortdf$Description <- factor(sortdf$Description,levels = sortdf$Description)

  sortdf$GeneRatio <- DOSE::parse_ratio(sortdf$GeneRatio)

  return(sortdf)

}





