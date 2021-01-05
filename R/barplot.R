#' Barplot of enrichResult
#' @param data  enrichResult object which is `data.frame` format.
#' @param color color to show
#' @param linetype linetype
#' @param linecol  line color
#' @param pointcol point color
#' @param filcol bar color
#' @param colby  mapping colors by pvalue or p.adjust
#' @importFrom graphics barplot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_fill_continuous
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_rect
#' @importFrom ggplot2 element_text
#'
#' @return A `ggplot` object
#' @export
#'
#' @examples
#'
#' data("AnalyzedSigPathway")
#' one_report <- load_AnalyzedSigPathway(1)
#' barplot(one_report)
barplot <- function(data,
                    color = c("blue", "red"),
                    linetype = "dashed",
                    linecol = "chocolate",
                    pointcol = "black",
                    filcol = "red",
                    colby = "pvalue") {
  if (!is.element(colby, c("pvalue", "p.adjust"))) {
    stop("the argument of distance should be one of pvalue, p.adjust.")
  }
  if (colby == "pvalue") {
    ylab <- expression("-log"[10] * "(P Value)")
  }
  if (colby == "p.adjust") {
    ylab <- expression("-log"[10] * "(P.adjust)")
  }

  mx <- max(-log10(data[, colby]))
  colbys <- data[, colby]
  ggplot(data = data, mapping = aes(
    x = Description,
    y = -log10(colbys),
    fill = -log10(colbys),
    group = 1
  )) +
    geom_bar(stat = "identity") +
    scale_fill_continuous(low = color[1], high = color[2]) +
    geom_line(aes(x = Description, y = -log10(colbys)),
      stat = "identity",
      linetype = linetype, color = linecol, size = 1.2
    ) +
    geom_point(shape = 21, color = pointcol, fill = filcol, size = 5) +
    scale_x_discrete(limits = rev(data$Description)) +
    ylim(0, mx) +
    labs(x = "", y = ylab) +
    coord_flip() +
    geom_text(aes(label = data$Count), hjust = -1, size = 5) +
    labs(fill = ylab) +
    theme_bw(15) +
    theme(
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "white"),
      panel.background = element_blank()
    ) +
    theme(
      title = element_text(
        size = 12, color = "black",
        face = "plain", hjust = 0.2, lineheight = 0.2
      ),
      axis.title.x = element_text(size = 12, face = "plain", color = "black", hjust = 0.5),
      axis.title.y = element_text(size = 14, color = "black", hjust = 0.5, angle = 45),
      axis.text.x = element_text(size = 8, color = "black"),
      axis.text.y = element_text(size = 12, face = "plain", color = "black")
    ) +
    theme(legend.position = c(1, 0.2), legend.justification = c(1, 0.2))
}
