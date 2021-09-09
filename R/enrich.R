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
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom cowplot plot_grid
#' @examples
#' one_report <- tcm.LoadAnalyzedDEG(2)
#' msigdb.path <- system.file("extdata", "c2.cp.kegg.v6.2.symbols.gmt",
#'   package = "TCMR", mustWork = TRUE
#' )
#' p <- tcm.EnrichGSEA(
#'   data = one_report,
#'   dirct = "up",
#'   msigdb.path = msigdb.path,
#'   n.path = 6,
#'   p.cutoff = 0.05
#' )
tcm.EnrichGSEA <- function(data = NULL,
                           res.path = tempdir(),
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
  p1 <- ggplot(gsdata, aes_(x = ~x)) +
    xlab(NULL) +
    scale_x_continuous(expand = c(0, 0)) +
    geom_line(aes_(y = ~runningScore, color = ~Description), size = 1) +
    scale_color_manual(values = mycol) +
    geom_hline(yintercept = 0, lty = "longdash", lwd = 0.3) +
    ylab("Enrichment Score") +
    theme_bw() +
    theme(
      legend.position = c(0.8, 0.7), legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent")
    ) +
    theme(
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      plot.margin = margin(t = .2, r = .2, b = 0, l = .2, unit = "cm")
    )

  rel_heights <- c(1.5, .4, 0.2)

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
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  #
  #   v <- seq(1, sum(gsdata$position), length.out=10)
  #   inv <- findInterval(rev(cumsum(gsdata$position)), v)
  #   if (min(inv) == 0) inv <- inv + 1
  #
  #   col <- c(rev(RColorBrewer::brewer.pal(5, "Blues")), RColorBrewer::brewer.pal(5, "Reds"))
  #
  #   ymin <- min(p2$data$ymin)
  #   yy <- max(p2$data$ymax - p2$data$ymin) * .3
  #   xmin <- which(!duplicated(inv))
  #   xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
  #   d <- data.frame(ymin = ymin, ymax = yy,
  #                   xmin = xmin,
  #                   xmax = xmax,
  #                   col = col[unique(inv)])
  # p3 <- ggplot()+ geom_rect(
  #   aes_(xmin=~xmin,
  #        xmax=~xmax,
  #        ymin=~ymin,
  #        ymax=~ymax,
  #        fill=~I(col)),
  #   data=d,
  #   inherit.aes=FALSE)+
  #   xlab(NULL) +
  #   ylab(NULL) +
  #   theme_bw() +
  #   theme(panel.grid = element_blank()) +
  #   theme(
  #     legend.position = "none",
  #     plot.margin = margin(t = -.1, b = 0, unit = "cm"),
  #     axis.ticks = element_blank(),
  #     axis.text = element_blank(),
  #     axis.line.x = element_blank()
  #   ) +
  #   scale_x_continuous(expand = c(0, 0))+
  #   geom_text()
  #
  mydf <- data.frame(sales = 1:100)

  p3 <- ggplot(mydf) +
    geom_tile(aes(x = 1, y = sales, fill = sales, alpha = .9)) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 50) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    theme(
      legend.position = "none",
      plot.margin = margin(t = -.1, b = 0, unit = "cm"),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line.x = element_blank()
    ) +
    xlab(NULL) +
    ylab(NULL) +
    geom_text(aes(x = 1, y = 92, label = "Negative"), size = 4.5) +
    geom_text(aes(x = 1, y = 8, label = "Positive"), size = 4.5) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
  #
  # df2 <- p1$data
  # df2$y <- p1$data$geneList[df2$x]
  # df2$gsym <- p1$data$gsym[df2$x]
  #
  # p4 <- ggplot(df2,
  #                 aes_string(x = "x",
  #                            y = "y",
  #                            fill = "Description",
  #                            color = "Description",
  #                            label = "gsym")) +
  #
  #   geom_bar(position = "dodge", stat = "identity") +
  #   scale_fill_manual(values = mycol, guide = FALSE) +
  #   scale_color_manual(values = mycol, guide = FALSE) +
  #
  #   # scale_x_continuous(expand=c(0,0)) +
  #   geom_hline(yintercept = 0, lty = 2, lwd = 0.2) +
  #   ylab("Ranked list metric") +
  #   xlab("Rank in ordered dataset") +
  #   theme_bw() +
  #   theme(
  #     axis.title.y = element_text(size = 15, face = "bold"),
  #     axis.title.x = element_text(size = 15, face = "bold"),
  #     axis.text.y = element_text(size = 10, face = "bold"),
  #     panel.grid = element_blank()
  #   ) +
  #   theme(plot.margin = margin(t = -.1, r = .2, b = .2, l = .2, unit = "cm"))

  plotlist <- list(p1, p2, p3)
  n <- length(plotlist)

  hm <- cowplot::plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)

  return(list(gsea.list = gsea.list, fig = hm))
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


#' Get Enriched GO Sets
#'
#' @inheritParams tcm.EnrichPathway
#' @inheritParams clusterProfiler::enrichGO
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' data("AnalyzedDEG")
#' one_report <- tcm.LoadAnalyzedDEG(1)
#' deg <- subset(one_report, abs(logFC) > 1 & P.Value < 0.05)
#' go <- tcm.EnrichGO(deg, p.cutoff = 0.05, p.adj.cutoff = 0.5, n.path = 10)
#' @testexamples
#' expect_is(go, "data.frame")

tcm.EnrichGO <- function(data,
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
  if (is.numeric(p.cutoff) && is.numeric(p.adj.cutoff)) {
    df <- df[df$pvalue < p.cutoff, ]
    df <- df[df$p.adjust < p.adj.cutoff, ]
  }

  if (is.numeric(p.cutoff) && is.null(p.adj.cutoff)) {
    df <- df[df$pvalue < p.cutoff, ]
  }
  if (is.null(p.cutoff) && is.numeric(p.adj.cutoff)) {
    df <- df[df$p.adjust < p.adj.cutoff, ]
  }
  df <- df[order(df[, sortby], df$GeneRatio), ]
  df$Description <- factor(df$Description, levels = df$Description)
  df$GeneRatio <- DOSE::parse_ratio(df$GeneRatio)

  if (is.numeric(n.path)) {
    if (nrow(df) > n.path) {
      df <- df[1:n.path, ]
    }
  }
  return(df)
}


#' Get Enriched Signal Pathways
#'
#' @param data A `data.frame` from a subset of [tcm.RunDEG] result.
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
#' one_report <- tcm.LoadAnalyzedDEG(1)
#' deg <- subset(one_report, abs(logFC) > 1 & P.Value < 0.05)
#' path <- tcm.EnrichPathway(deg, p.cutoff = 1, p.adj.cutoff = 1, n.path = 10)
tcm.EnrichPathway <- function(data,
                              p.cutoff = 0.05,
                              p.adj.cutoff = 0.2,
                              sortby = "pvalue",
                              n.path = 10) {
  options(warn = -1)
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
#' \dontrun{
#' tfdata <- tcm.LoadData("tfdata")
#' gene <- tfdata$gene[100:500]
#' tcm.EnrichTF(gene, tfdata = tfdata, colours = c("red", "blue"), colby = "p.adjust")
#' tcm.EnrichTF(gene, tfdata = tfdata, colours = c("red", "blue"), colby = "pvalue")
#' }
tcm.EnrichTF <- function(gene,
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
  utils::write.table(egmt2, file = paste0(res.path, "/tcm.EnrichTF_results.txt"), row.names = TRUE, col.names = NA, sep = "\t", quote = FALSE)
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
