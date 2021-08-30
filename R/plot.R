#' Barplot for enrichResult
#'
#' @param data  enrichResult result in `data.frame` format
#' @param color Color to show
#' @param linetype Linetype
#' @param linecol  Line color
#' @param pointcol Point color
#' @param filcol Bar color
#' @param colby  Mapping colors by pvalue or p.adjust
#' @importFrom graphics barplot
#' @importFrom ggplot2 ggplot geom_bar aes scale_fill_continuous theme geom_point geom_line ylim labs geom_text theme_bw
#' scale_x_discrete coord_flip element_line element_blank element_rect element_text
#'
#' @return A `ggplot` object
#' @export
#'
#' @examples
#' data("AnalyzedSigPathway")
#' one_report <- tcm.LoadAnalyzedSigPathway(1)
#' tcm.EnrichBarplot(one_report)
tcm.EnrichBarplot <- function(data,
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



#' Visualize Connectivity score
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
#' gcmap_kk <- tcm.SSwithGCMAP(input = query2, data = data_logFC)
#' test <- head(gcmap_kk, 20)
#' cs_plot(test, x = "raw_score", y = "tcm", colby = "pvalue", color = c("blue", "red"))
cs_plot <- function(data,
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
    theme(
      panel.background = element_rect(
        colour = "black",
        size = 0.5
      ),
      title = element_text(
        size = 12, color = "black",
        face = "plain", hjust = 0.2, lineheight = 0.2
      ),
      axis.title.x = element_text(size = 12, face = "plain", color = "black", hjust = 0.5),
      axis.title.y = element_text(size = 14, color = "black"),
      axis.text.x = element_text(size = 12, face = "plain", color = "black"),
      axis.text.y = element_text(size = 12, face = "plain", color = "black")
    ) +
    labs(
      x = "Score",
      size = "Score",
      col = colby
    )
  p
}


#' Dotplot of Special Genes in Special Compounds of logFC Matrix
#'
#' @param data Special genes expression data of special compounds
#' @param color Color to distinguish between different trends
#'
#' @return A `ggplot` object
#' @export
#' @importFrom reshape2 melt
#' @examples
#' data <- as.matrix(data_logFC[1:6, ])
#' data <- data.frame(t(data), check.rows = FALSE, check.names = FALSE)
#' data$compound <- rownames(data)
#' data <- data[3:13, ]
#' gene_dot(data)
gene_dot <- function(data,
                     color = c("red", "blue")) {
  tmp <- reshape2::melt(data = data, id.vars = "compound")

  tmp$type <- ifelse(tmp$value > 0, "UP", "Down")
  tmp$expression <- abs(tmp$value)

  ggplot2::ggplot(tmp, aes_string(x = "variable", y = "compound", colour = "type")) +
    geom_point(aes(size = expression)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) +
    scale_color_manual(values = c(UP = "blue", Down = "red")) +
    labs(x = NULL, y = NULL) +
    theme(
      title = element_text(
        size = 12, color = "black",
        face = "plain", hjust = 0.2, lineheight = 0.2
      ),
      axis.title.x = element_text(size = 12, face = "plain", color = "black", hjust = 0.5),
      axis.title.y = element_text(size = 14, color = "black", hjust = 0.5, angle = 45),
      axis.text.x = element_text(size = 8, color = "black"),
      axis.text.y = element_text(size = 12, face = "plain", color = "black")
    )
}

#' Title Gene Expression Plot
#'
#' @param data Expression data of compounds
#' @param num Numbers of cluster
#' @param colors Color to show expression level
#' @param shownames Logic value to determine whether to display the gene name
#'
#' @return A plot
#' @export
#' @examples
#' tcm.CircleHeatmap(data = as.matrix(data_logFC[1:6, ]))
#' tcm.CircleHeatmap(data = as.matrix(data_logFC[1:6, ]), colors = c("blue", "white", "red"))
#' tcm.CircleHeatmap(data = as.matrix(data_logFC[1:10, ]), colors = c("skyblue", "white", "red"))
tcm.CircleHeatmap <- function(data, num = 3, colors = c("blue", "white", "red"), shownames = TRUE) {
  lable1 <- row.names(data)
  lable2 <- substring(colnames(data), 1, 10)
  dend <- stats::as.dendrogram(stats::hclust(stats::dist(t(as.matrix(data)))))
  dend <- dend %>% dendextend::set("branches_k_color", k = num)
  col_fun <- circlize::colorRamp2(
    breaks = c(min(data), median(data), max(data)),
    colors = colors
  )
  col_mat <- col_fun(data)
  mat <- data[, order.dendrogram(dend)]
  nr <- nrow(mat)
  nc <- ncol(mat)
  par(mar = c(0, 0, 0, 0))
  circlize::circos.clear()

  circlize::circos.par(
    canvas.xlim = c(-1.3, 1.3),
    canvas.ylim = c(-1.3, 1.3),
    cell.padding = c(0, 0, 0, 0),
    gap.degree = 90
  )
  factors <- "a"
  circlize::circos.initialize(factors, xlim = c(0, ncol(mat)))

  circlize::circos.track(
    ylim = c(0, nr), bg.border = NA, track.height = ifelse(nr < 8, 0.1 * nr, 5 / nr),
    panel.fun = function(x, y) {
      for (i in 1:nr) {
        circlize::circos.rect(
          xleft = 1:nc - 1, ybottom = rep(nr - i, nc),
          xright = 1:nc, ytop = rep(nr - i + 1, nc),
          border = "white",
          col = col_mat[i, ]
        )
        circlize::circos.text(
          x = nc,
          y = (nr + 0.4) - i,
          labels = lable1[i],
          facing = "downward", niceFacing = TRUE,
          cex = 0.6,
          adj = c(-0.2, 0)
        )
      }
    }
  )
  if (shownames) {
    suppressMessages(for (i in 1:nc) {
      circlize::circos.text(
        x = i - 0.4,
        y = nr + 1,
        labels = lable2[i],
        facing = "clockwise", niceFacing = TRUE,
        cex = 0.6, adj = c(0, 0)
      )
    })
  } else {
    suppressMessages(for (i in 1:nc) {
      circlize::circos.text(
        x = i - 0.4,
        y = nr + 1,
        labels = NULL,
        facing = "clockwise", niceFacing = TRUE,
        cex = 0.6, adj = c(0, 0)
      )
    })
  }



  max_height <- max(attr(dend, "height"))
  circlize::circos.track(
    ylim = c(0, max_height), bg.border = NA, track.height = 0.2,
    panel.fun = function(x, y) {
      circlize::circos.dendrogram(
        dend = dend,
        max_height = max_height
      )
    }
  )
  circlize::circos.clear()

  lgd <- ComplexHeatmap::Legend(
    at = c(-5, -2, 0, 2, 5), col_fun = col_fun,
    title_position = "topcenter", title = "Z-score"
  )
  ComplexHeatmap::draw(lgd, x = unit(0.7, "npc"), y = unit(0.7, "npc"))
}


#' Dotplot of enrichResult
#'
#' @param data enrichResult object which is `data.frame` format.
#' @param fill dot color
#' @param color color to show
#'
#' @return A `ggplot` object
#' @export
#' @importFrom ggplot2 aes_
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @examples
#' data("AnalyzedSigPathway")
#' one_report <- tcm.LoadAnalyzedSigPathway(1)
#' tcm.EnrichDotplot(one_report)
tcm.EnrichDotplot <- function(data,
                              fill = "pvalue",
                              color = c("blue", "red")) {
  if (!is.element(fill, c("pvalue", "p.adjust"))) {
    stop("the argument of distance should be one of pvalue, p.adjust.")
  }

  fills <- data[, fill]
  min_x <- min(data$GeneRatio)
  max_x <- max(data$GeneRatio)

  if (fill == "pvalue") {
    collab <- expression("-log"[10] * "(P Value)")
  }
  if (fill == "p.adjust") {
    collab <- expression("-log"[10] * "(P.adjust)")
  }

  ggplot(data, aes(GeneRatio, Description, colour = -log10(fills))) +
    geom_point(aes(size = Count)) +
    scale_color_gradientn(colours = color) +
    scale_size_continuous(range = c(2, 10)) +
    scale_x_continuous(limits = c(min_x, max_x)) +
    theme_bw(15) +
    ylab("") +
    labs(size = "Count", colour = collab) +
    theme(legend.position = c(1, 0.5), legend.justification = c(1, 0.5)) +
    theme(legend.background = element_blank()) +
    theme(legend.key = element_blank()) +
    theme(
      title = element_text(
        size = 12, color = "black",
        face = "plain", hjust = 0.2, lineheight = 0.2
      ),
      axis.title.x = element_text(size = 12, face = "plain", color = "black", hjust = 0.5),
      axis.title.y = element_text(size = 14, color = "black", hjust = 0.5, angle = 45),
      axis.text.x = element_text(size = 8, color = "black"),
      axis.text.y = element_text(size = 12, face = "plain", color = "black")
    )
}
