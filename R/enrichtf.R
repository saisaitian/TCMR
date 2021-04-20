#' TF Enrichment for Genes
#'
#' @param gene The special gene for enrichment
#' @param tfdata Transcription factor and gene correspondence table
#' @param colby Mapping colors by `pvalue` or `p.adjust`
#' @param colours Color to show
#' @param n An integer value to indicate how many top tfs sorted by GeneRatio should be identified; 20 by default.
#' @param res.path A string value to indicate the path for saving the results for functional pathways.
#'
#' @return A `ggplot` object
#' @export
#'
#' @examples
#' data("tfdata")
#' gene <- tfdata$gene[100:500]
#' tf_enrich(gene, tfdata = tfdata, colours = c("red", "blue"), colby = "p.adjust")
#' tf_enrich(gene, tfdata = tfdata, colours = c("red", "blue"), colby = "pvalue")
tf_enrich <- function(gene,
                      tfdata = tfdata,
                      colby = "pvalue",
                      res.path = getwd(),
                      colours = c("red", "yellow"),
                      n = 20) {
  if (!is.element(colby, c("pvalue", "p.adjust"))) {
    stop("the argument of distance should be one of pvalue, p.adjust.")
  }
  egmt <- clusterProfiler::enricher(gene, TERM2GENE = tfdata)

  egmt2 <- data.frame(egmt)
  utils::write.table(egmt2, file = paste0(res.path, "/tf_enrich_results.txt"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
  egmt2$GeneRatio <- DOSE::parse_ratio(egmt2$GeneRatio)

  sortdf <- egmt2[order(egmt2$GeneRatio), ]
  sortdf$Description <- factor(sortdf$Description, levels = sortdf$Description)
  sortdf <- sortdf[1:n, ]

  if (colby == "pvalue") {
    pvalue <- sortdf[, colby]
    p <- ggplot(sortdf, aes(GeneRatio, Description,
      colour = pvalue
    )) + # 横坐标、纵坐标、颜色代表p-value
      geom_point(aes(size = Count)) + # 圆点的大小代表组内基因数
      scale_color_gradientn(colours = colours) + # 可以自己改颜色
      # 圆点的大小代表Number of significant genes
      scale_size_continuous(range = c(2, 10)) +
      theme_bw(15) +
      ylab("Transcription factor") +
      theme(
        axis.title.x = element_text(size = 12, face = "plain", color = "black", hjust = 0.5),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 12, face = "plain", color = "black")
      )
  } else {
    p.adjust <- sortdf[, colby]

    p <- ggplot(sortdf, aes(GeneRatio, Description,
      colour = p.adjust
    )) + # 横坐标、纵坐标、颜色代表p-value
      geom_point(aes(size = Count)) + # 圆点的大小代表组内基因数
      scale_color_gradientn(colours = colours) + # 可以自己改颜色
      # 圆点的大小代表Number of significant genes
      scale_size_continuous(range = c(2, 10)) +
      theme_bw(15) +
      ylab("Transcription factor") +
      theme(
        axis.title.x = element_text(size = 12, face = "plain", color = "black", hjust = 0.5),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 12, face = "plain", color = "black")
      )
  }

  p
}
