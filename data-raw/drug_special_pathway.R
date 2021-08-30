
#' #' Title drug special pathway analysis
#'
#' @param index the index of colunm in data_logFC matrix
#' @param group the compare variable including group and class
#' @param num the number of pathways to show
#' @param colorby Mapping colors by pvalue or p.adjust
#' @param res.path A string value to indicate the path for saving the results for functional pathways.
#' @param plot  A logic value to indicate if to plot pathways; TRUE by default.
#'
#' @return data frame
#' @export
#'
#' @examples
#' data(data_logFC)
#' drug_special_pathway(
#'   group = "drug",
#'   data = data_logFC
#' )
drug_special_pathway <- function(index = c(1, 3, 5),
                                 data = data_logFC,
                                 group = "drug",
                                 num = 10,
                                 colorby = "pvalue",
                                 res.path = getwd(),
                                 plot = T) {
  message(paste0("\n", ">>> Calculating drug special KEGG pathway"))

  tmp <- data[, index]

  compname <- names(tmp)

  tmp <- tibble::rownames_to_column(tmp)

  newdata <- reshape2::melt(tmp, id.vars = "rowname")

  newdata <- newdata[abs(newdata$value) > 1, ]

  names(newdata)[1] <- "SYMBOL"

  eg <- suppressWarnings(suppressMessages(clusterProfiler::bitr(newdata$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = "org.Hs.eg.db"
  )))

  newdata <- merge(eg, newdata, by = "SYMBOL")

  newdata$class <- "up"

  newdata$class[newdata$value < 0] <- "down"

  names(newdata)[3] <- "drug"

  if (length(group) == 2) {
    formula_res <- suppressWarnings(compareCluster(ENTREZID ~ eval(parse(text = group[1])) + drug + class,
      data = newdata,
      fun = "enrichKEGG",
      pvalueCutoff = 0.05
    ))
  }

  if (group == "drug") {
    formula_res <- suppressWarnings(clusterProfiler::compareCluster(ENTREZID ~ drug,
      data = newdata,
      fun = "enrichKEGG",
      pvalueCutoff = 0.05
    ))
  }


  if (group == "class") {
    formula_res <- suppressWarnings(clusterProfiler::compareCluster(ENTREZID ~ class,
      data = newdata,
      fun = "enrichKEGG",
      pvalueCutoff = 0.05
    ))
  }

  formula_res <- DOSE::setReadable(formula_res, org.Hs.eg.db, keyType = "ENTREZID")

  outfile <- file.path(res.path, paste("compare_KEGG", "_result.", ".txt", sep = ""))

  write.table(as.data.frame(formula_res), file = outfile, row.names = FALSE, sep = "\t", quote = FALSE)

  if (plot == T) {
    enrichplot::tcm.EnrichDotplot(formula_res, color = colorby, showCategory = num)
  }
}
