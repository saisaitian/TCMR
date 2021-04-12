#' Run GSEA
#' @param data a deg data.frame
#' @param res.path A string value to indicate the path for saving the results for functional pathways.
#' @param dirct A string value to indicate the direction of identifying significant pathway. Allowed values contain c('up', 'down'); `up` means up-regulated pathway, and `down` means down-regulated pathway; "up" by default.
#' @param n.path A integer value to indicate how many top pathways sorted by NES should be identified; 6 by default.
#' @param msigdb.path A string value to indicate ABSOULUTE PATH/NAME of MSigDB file (GMT file with gene symbols) downloaded from \url{https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#H}.
#' @param nPerm A integer value to indicate the number of permutations; 1000 by default and 10000 will be better for reproducibility.
#' @param minGSSize A integer value to indicate minimal size of each geneSet for analyzing; 10 by default
#' @param maxGSSize A integer value to indicate maximal size of each geneSet for analyzing; 500 by default.
#' @param p.cutoff A numeric value to indicate the nominal p value for identifying significant pathways; pvalue < 0.05 by default.
#' @param p.adj.cutoff A numeric value to indicate the adjusted p value for identifying significant pathways; padj < 0.05 by default.
#'
#' @return A plot.
#' @export
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 margin
#' @importFrom ggplot2 geom_linerange
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom ggplot2 labs
#' @importFrom cowplot plot_grid
#' @examples
#' one_report <- load_analyzedDEG(2)
#' msigdb.path <- system.file("extdata", "c2.cp.kegg.v6.2.symbols.gmt",
#'   package = "TCMR", mustWork = TRUE
#' )
#' runGSEA(
#'   data = one_report,
#'   dirct = "up",
#'   msigdb.path = msigdb.path,
#'   n.path = 3,
#'   p.cutoff = 0.05
#' )
runGSEA <- function(data = NULL,
                    res.path = getwd(),
                    dirct = "up",
                    n.path = 6,
                    msigdb.path = NULL,
                    nPerm = 1000,
                    minGSSize = 10,
                    maxGSSize = 500,
                    p.cutoff = 0.05,
                    p.adj.cutoff = NULL) {
    geneList <- data$logFC
    names(geneList) <- data$identifier
  geneList <- sort(geneList, decreasing = TRUE) # ranked gene set

  # run gsea
  msigdb <- try(clusterProfiler::read.gmt(msigdb.path), silent = TRUE)
  if (class(msigdb) == "try-error") {
    stop("please provide correct ABSOLUTE PATH for MSigDB file.")
  }

  gsea.list <- suppressWarnings(clusterProfiler::GSEA(
    geneList = geneList,
    TERM2GENE = msigdb,
    nPerm = nPerm,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    seed = TRUE,
    verbose = FALSE,
    pvalueCutoff = 1
  ))
  gsea.dat <- as.data.frame(gsea.list)
  gseaidList <- c()
  if (dirct == "up") {
    if (!is.null(p.cutoff)) {
      gsea.dat <- gsea.dat[which(gsea.dat$pvalue < p.cutoff), ]
    }

    if (!is.null(p.adj.cutoff)) {
      gsea.dat <- gsea.dat[which(gsea.dat$p.adjust < p.adj.cutoff), ]
    }

    if (!is.null(p.cutoff) & !is.null(p.adj.cutoff)) {
      gsea.dat <- gsea.dat[which(gsea.dat$pvalue < p.cutoff & gsea.dat$p.adjust < p.adj.cutoff), ]
    }
    utils::write.table(gsea.dat, file = paste0(res.path, "/gsea_up_results.txt"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
    gsea.dat <- gsea.dat[order(gsea.dat$NES, decreasing = TRUE), ]

    if (nrow(gsea.dat) > n.path) {
      gsea.dat <- gsea.dat[1:n.path, ]
    } else {
      gsea.dat <- gsea.dat[1:n.path, ]
    }

    geneSetID <- rownames(gsea.dat)
  }

  if (dirct == "down") {
    if (!is.null(p.cutoff)) {
      gsea.dat <- gsea.dat[which(gsea.dat$pvalue < p.cutoff), ]
    }

    if (!is.null(p.adj.cutoff)) {
      gsea.dat <- gsea.dat[which(gsea.dat$p.adjust < p.adj.cutoff), ]
    }

    if (!is.null(p.cutoff) & !is.null(p.adj.cutoff)) {
      gsea.dat <- gsea.dat[which(gsea.dat$pvalue < p.cutoff & gsea.dat$p.adjust < p.adj.cutoff), ]
    }

    utils::write.table(gsea.dat, file = paste0(res.path, "/gsea_down_results.txt"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
    gsea.dat <- gsea.dat[order(gsea.dat$NES, decreasing = TRUE), ]
    if (nrow(gsea.dat) > n.path) {
      gsea.dat <- gsea.dat[1:n.path, ]
    } else {
      gsea.dat <- gsea.dat[1:n.path, ]
    }
    geneSetID <- rownames(gsea.dat)
  }
  message("GSEA done...")
  gsdata <- do.call(rbind, lapply(geneSetID, .gsInfo, object = gsea.list))

  gsym.fc.id.sorted <- names(geneList)

  gsdata$gsym <- rep(gsym.fc.id.sorted, nrow(gsea.dat))

  mycol <- c("darkgreen", "chocolate4", "blueviolet", "#223D6C", "#D20A13", "#088247", "#58CDD9", "#7A142C", "#5D90BA", "#431A3D", "#91612D", "#6E568C", "#E0367A", "#D8D155", "#64495D", "#7CC767")
  mycol <- mycol[4:(3 + length(geneSetID))]
  p.res <- ggplot(gsdata, aes_(x = ~x)) +
    xlab(NULL) +
    geom_line(aes_(y = ~runningScore, color = ~Description), size = 1) +
    scale_color_manual(values = mycol) +
    geom_hline(yintercept = 0, lty = "longdash", lwd = 0.2) +
    ylab("Enrichment\n Score") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(
      legend.position = "top", legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent")
    ) +
    theme(
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      plot.margin = margin(t = .2, r = .2, b = 0, l = .2, unit = "cm")
    )

  rel_heights <- c(1.5, .5, 1.5)

  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin = ~ymin, ymax = ~ymax, color = ~Description)) +
    xlab(NULL) +
    ylab(NULL) +
    scale_color_manual(values = mycol) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(
      legend.position = "none",
      plot.margin = margin(t = -.1, b = 0, unit = "cm"),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line.x = element_blank()
    ) +
    # scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand = c(0, 0))

  df2 <- p.res$data
  df2$y <- p.res$data$geneList[df2$x]
  df2$gsym <- p.res$data$gsym[df2$x]

  p.pos <- ggplot(
    df2,
    aes_string(x = "x", y = "y", fill = "Description", color = "Description", label = "gsym")
  ) +
    geom_segment(
      data = df2, aes_(x = ~x, xend = ~x, y = ~y, yend = 0),
      color = "grey"
    ) +
    geom_bar(position = "dodge", stat = "identity") +
    scale_fill_manual(values = mycol, guide = FALSE) +
    scale_color_manual(values = mycol, guide = FALSE) +

    # scale_x_continuous(expand=c(0,0)) +
    geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
    ylab("Ranked list\n metric") +
    xlab("Rank in ordered dataset") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 12, face = "bold"),
      panel.grid = element_blank()
    ) +
    theme(plot.margin = margin(t = -.1, r = .2, b = .2, l = .2, unit = "cm"))

  plotlist <- list(p.res, p2, p.pos)
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(
      axis.line.x = element_line(),
      axis.ticks.x = element_line(),
      axis.text.x = element_text(size = 12, face = "bold")
    )
  cowplot::plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
}



.gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  if (is.numeric(geneSetID)) {
    geneSetID <- object@result[geneSetID, "ID"]
  }
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- .gseaScores(geneList, geneSet, exponent, fortify = TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore)) / 20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}




.gseaScores <- function(geneList, geneSet, exponent = 1, fortify = FALSE) {
  geneSet <- intersect(geneSet, names(geneList))

  N <- length(geneList)
  Nh <- length(geneSet)

  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet ## logical

  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit / NR)

  Pmiss[!hits] <- 1 / (N - Nh)
  Pmiss <- cumsum(Pmiss)

  runningES <- Phit - Pmiss

  ## ES is the maximum deviation from zero of Phit-Pmiss
  max.ES <- max(runningES)
  min.ES <- min(runningES)
  if (abs(max.ES) > abs(min.ES)) {
    ES <- max.ES
  } else {
    ES <- min.ES
  }

  df <- data.frame(
    x = seq_along(runningES),
    runningScore = runningES,
    position = as.integer(hits)
  )

  if (fortify == TRUE) {
    return(df)
  }

  df$gene <- names(geneList)
  res <- list(ES = ES, runningES = df)
  return(res)
}



