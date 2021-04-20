












library(ComplexHeatmap, quietly = T)
library(pheatmap, quietly = T)
library(circlize, quietly = T)

mat <- data_logFC
mat <- as.matrix(mat)
# mat<-log10(mat+1)
cluster_rows <- T
cluster_cols <- T
scale <- "none"
a <- pheatmap(mat, cluster_rows = T, cluster_cols = T, scale = "none", silent = T)

sample_num <- ncol(mat)
if (cluster_col == "TRUE" && cluster_row == "TRUE") {
  samplenames <- rev(colnames(mat[a$tree_row$order, a$tree_col$order]))
  genenames <- rownames(mat[a$tree_row$order, a$tree_col$order])
}
if (cluster_col == "TRUE" && cluster_row == "FALSE") {
  samplenames <- rev(colnames(mat[, a$tree_col$order]))
  genenames <- rownames(mat[, a$tree_col$order])
}
if (cluster_col == "FALSE" && cluster_row == "TRUE") {
  samplenames <- rev(colnames(mat[a$tree_row$order, ]))
  genenames <- rownames(mat[a$tree_row$order, ])
}
if (cluster_col == "FALSE" && cluster_row == "FALSE") {
  samplenames <- rev(colnames(mat))
  genenames <- rownames(mat)
}
if (gcex == 0) {
  genenames <- rep("", ncol(mat))
}
if (scex == 0) {
  samplenames <- rep("", nrow(mat))
}

if (scale == "none" && cluster_col == "TRUE" && cluster_row == "TRUE") {
  mat_cluster <- mat[a$tree_row$order, a$tree_col$order]
}
if (scale == "row" && cluster_col == "TRUE" && cluster_row == "TRUE") {
  mat_scale <- t(apply(mat, 1, "scale"))
  mat_cluster <- mat_scale[a$tree_row$order, a$tree_col$order]
}
if (scale == "none" && cluster_col == "TRUE" && cluster_row == "FALSE") {
  mat_cluster <- mat[, a$tree_col$order]
}
if (scale == "row" && cluster_col == "TRUE" && cluster_row == "FALSE") {
  mat_scale <- t(apply(mat, 1, "scale"))
  mat_cluster <- mat_scale[, a$tree_col$order]
}
if (scale == "none" && cluster_col == "FALSE" && cluster_row == "TRUE") {
  mat_cluster <- mat[a$tree_row$order, ]
}
if (scale == "row" && cluster_col == "FALSE" && cluster_row == "TRUE") {
  mat_scale <- t(apply(mat, 1, "scale"))
  mat_cluster <- mat_scale[a$tree_row$order, ]
}
if (scale == "none" && cluster_col == "FALSE" && cluster_row == "FALSE") {
  mat_cluster <- mat
}
if (scale == "row" && cluster_col == "FALSE" && cluster_row == "FALSE") {
  mat_scale <- t(apply(mat, 1, "scale"))
  mat_cluster <- mat_scale
}

if (color == 1) {
  my_color <- c("#005522", "#FFFACD", "firebrick3")
} else if (color == 2) {
  my_color <- c("navy", "white", "firebrick3")
}

col_fun <- colorRamp2(c(min(mat_cluster), median(mat_cluster), max(mat_cluster)), my_color)
mat <- t(mat_cluster)
factors <- rep("a", times = ncol(mat))
mat_list <- list(a = mat[, factors == "a"])
dend_list <- list(a = as.dendrogram(hclust(dist(t(mat_list[["a"]])))))
dend_list2 <- list(a = as.dendrogram(hclust(dist(mat_list[["a"]]))))

pdf(output)
circos.par(cell.padding = c(0, 0, 0, 0), track.margin = c(0, 0), gap.degree = 90, "track.height" = 0.55)
circos.initialize(factors, xlim = cbind(0, table(factors)))
circos.track(ylim = c(0, nrow(mat) + 6), bg.border = NA, panel.fun = function(x, y) {
  sector.index <- CELL_META$sector.index
  xlim <- CELL_META$xlim
  m <- mat_list[[sector.index]]
  dend <- dend_list[[sector.index]]
  m2 <- m[, order.dendrogram(dend)]
  if (cluster_row == "FALSE") {
    m2 <- m
  }
  col_mat <- col_fun(m2)
  nr <- nrow(m2)
  nc <- ncol(m2)
  col <- col_mat[1, ]
  for (i in 1:nr) {
    circos.rect(1:nc - 1, rep(nr - i, nc),
      1:nc, rep(nr - i + 1, nc),
      border = col_mat[i, ], col = col_mat[i, ]
    )
    circos.text(rep(xlim[2] + sgap, nr), 1:nr - 0.4, samplenames, adj = c(0, 0.5), cex = scex, niceFacing = TRUE)
  }
  for (i in 1:nc) {
    if (270 / nc * (i - 0.5) < 90 || 270 / nc * (i - 0.5) == 90) {
      circos.text(i - 0.5, CELL_META$cell.ylim[2] - 6 + uy(ggap, "mm"), genenames[i], facing = "clockwise", niceFacing = TRUE, cex = gcex, srt = -i * (270 / nc), adj = c(0, 0.5), font = 2)
    }
    if (270 / nc * (i - 0.5) > 90) {
      circos.text(i - 0.5, CELL_META$cell.ylim[2] - 6 + uy(ggap, "mm"), genenames[i], facing = "clockwise", niceFacing = TRUE, cex = gcex, srt = -i * (270 / nc) + 180, adj = c(0, 0.5), font = 2)
    }
  }
})

sector.index <- CELL_META$sector.index
m <- mat_list[[sector.index]]
dend <- dend_list[[sector.index]]
dend2 <- dend_list2[[sector.index]]
m2 <- m[, order.dendrogram(dend)]
col_mat <- col_fun(m2)
col <- col_mat[1, ]
genenames <- names(col)
max_height <- max(sapply(dend_list, function(x) attr(x, "height")))

if (cluster_row == "TRUE") {
  circos.track(
    ylim = c(0, max_height + 0.3), bg.border = NA, track.height = 0.3,
    panel.fun = function(x, y) {
      sector.index <- get.cell.meta.data("sector.index")
      dend <- dend_list[[sector.index]]
      circos.dendrogram(dend, max_height = max_height)
    }
  )
}

lgd <- Legend(col_fun = col_fun, title_position = "topcenter", grid_width = unit(0.45, "cm"), legend_height = unit(2.3, "cm"), labels_gp = gpar(cex = 0.6, font = 2))
draw(lgd, x = unit(0.8, "npc"), y = unit(0.75, "npc"))
if (cluster_col == "TRUE") {
  par(new = TRUE, pin = c(2, 1), plt = c(0.695, 0.815, 0.5, 0.56))
  plot(dend2, axes = FALSE, xaxs = "i", leaflab = "none", xlim = c(2.6, -1.5))
}
dev.off()
