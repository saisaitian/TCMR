library(plyr)
library(stringr)
library(ape)
library(GOSemSim)
library(ggtree)
library(scales)
library(cowplot)

library(simplifyEnrichment)
# go <- tcm.EnrichGO(deg, p.cutoff =NULL, p.adj.cutoff = 0.05, n.path = 20)

GOcluster <- function(data,
                      num = 4,
                      fontsize = 4,
                      width = 0.3,
                      low = "red",
                      high = "white") {
  GOid <- data$ID
  mat <- simplifyEnrichment::GO_similarity(GOid, ont = "BP")
  rownames(mat) <- go$Description
  colnames(mat) <- go$Description
  Pmat <- data.frame(pvalue = data$pvalue, p.adjust = data$p.adjust, stringsAsFactors = F)
  rownames(Pmat) <- go$Description
  tree <- ape::nj(as.dist(1 - mat))
  node <- seq(min(tree$edge[, 1]), length.out = num)

  gtree <- tidytree::groupClade(tree, .node = node)

  pbase <- ggtree::ggtree(
    gtree,
    aes(color = group)
  )
  pnode <- pbase +

    ggtree::geom_tiplab(size = fontsize, align = TRUE) +
    ggplot2::coord_cartesian(xlim = c(-.1, 1.5))

  p <- ggtree::gheatmap(pnode, smallmat,
    offset = .75,
    width = width,
    colnames_angle = 90, hjust = 0,
    low = low, high = high
  )

  p
}


# GOcluster(go)
