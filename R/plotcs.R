
#' Title Visualize Connectivity score
#'
#' @param data A data.frame from results of five methods
#' @param x x axis
#' @param y y axis
#' @param colby color sort by a variable
#' @param color color to show
#' @importFrom forcats fct_reorder
#' @importFrom stats median order.dendrogram
#' @importFrom graphics par
#' @importFrom ggplot2 scale_colour_gradientn theme_minimal unit
#' @return a ggplot
#' @export
#'
#' @examples
#' data("data_logFC")
#' query2 <- data_logFC[1:60, 1, drop = FALSE]
#' gcmap_kk <- ss_gcmap(input = query2, data = data_logFC)
#' test <- head(gcmap_kk, 20)
#' plotcs(test, x = "raw_score", y = "tcm", colby = "pvalue", color = c("blue", "red"))
plotcs <- function(data,
                   x = "scaled_score",
                   y = "tcm",
                   colby = "pvalue",
                   color = c("blue", "red")) {
  score <- data[, x]
  y <- data[, y]
  Pval <- data[, colby]
  p <- ggplot2::ggplot(data, aes(score, forcats::fct_reorder(y, score))) +
    geom_segment(aes(xend = 0, yend = y), linetype = 2) +
    geom_point(aes(col = Pval, size = abs(score))) +
    scale_colour_gradientn(colours = color) +
    scale_size_continuous(range = c(2, 6)) +
    ylab(NULL) +
    theme_minimal() +
    theme(panel.background = element_rect(
      colour = "black",
      size = 0.5
    )) +
    labs(
      x = "Score", y = "Compounds",
      size = "Score",
      col = colby
    )
  p
}
