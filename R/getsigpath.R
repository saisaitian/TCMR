#' Get Differential Pathway
#' - `getsigpath` is used to obtain differential pathway with p < p.cutoff .
#' @param data  A `data.frame` format data from `deg_caller` function.
#' @param p.cutoff the cutoff value of p.
#' @param p.adj.cutoff the cutoff value of p.adjust.
#' @param n.path A integer value to indicate how many top pathways sorted by pvalue should be identified; 10 by default.
#' @param sortby A string value to indicate the pathway was sorted by pvalue or p.adjust value. Allowed value contains c('pvalue', 'p.adjust').
#' @importFrom clusterProfiler bitr
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom clusterProfiler setReadable
#' @importFrom DOSE parse_ratio
#' @return  A `data.frame`.
#'
#' @export
#'
#' @examples
#' data("AnalyzedDEG")
#' one_report <- load_analyzedDEG(1)
#' deg <- subset(one_report,abs(logFC)>1&P.Value<0.05)
#' aa=getsigpath(deg,p.cutoff = 0.05,p.adj.cutoff = 0.5,n.path = 10)
#' aaa=getsigpath(deg,p.cutoff = 0.05,p.adj.cutoff = 0.5,sortby ='p.adjust' ,n.path = 10)
#'
getsigpath <- function(data,
                       p.cutoff = 0.05,
                       p.adj.cutoff = NULL,
                       sortby = 'pvalue',
                       n.path = 10){

  stopifnot(is.data.frame(data), !is.null(rownames(data)))
  identifiers = suppressWarnings(clusterProfiler::bitr(data$identifier, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))
  message("The number of gene is: ",paste0(nrow(identifiers), collapse = "  "),' for next analysis')
  kk <- suppressWarnings(enrichKEGG(gene= identifiers$ENTREZID,
                                    organism     = 'hsa',
                                    pvalueCutoff = 1,
                                    qvalueCutoff = 1))

  kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

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
