
#' Title gene expression plot
#'
#' @param data the expression data of compounds
#' @param num the numbers of cluster
#' @param colors the color to show expression level
#'
#' @return
#' @export
#'
#' @examples geneplot(data=as.matrix(data_logFC[1:6,]))
#'
#' @tests
#'

geneplot <- function(data,num=3,colors = c("blue", "white", "red")) {
  lable1=row.names(data)
  lable2= substring(colnames(data),1,10)
  dend <- stats::as.dendrogram(stats::hclust(stats::dist(t(as.matrix(data)))))
  dend <-dend %>% dendextend::set("branches_k_color", k = num)
  col_fun = circlize::colorRamp2(breaks = c(min(data),median(data),max(data)),
                       colors = colors)
  col_mat = col_fun(data)
  mat = data[, order.dendrogram(dend)]
  nr = nrow(mat)
  nc = ncol(mat)
  par(mar=c(0,0,0,0))
  circlize::circos.clear()

  circlize::circos.par(canvas.xlim =c(-1.3,1.3),
                            canvas.ylim = c(-1.3,1.3),
                            cell.padding = c(0,0,0,0),
                            gap.degree =90)
  factors = "a"
  circlize::circos.initialize(factors, xlim = c(0, ncol(mat)))

  circlize::circos.track(ylim = c(0, nr),bg.border = NA,track.height = 0.1*nr,
               panel.fun = function(x, y) {
                 for(i in 1:nr) {
                   circlize::circos.rect(xleft = 1:nc - 1, ybottom = rep(nr - i, nc),
                               xright = 1:nc, ytop = rep(nr - i + 1, nc),
                               border = "white",
                               col = col_mat[i,])
                   circlize::circos.text(x = nc,
                               y = 6.4 -i,
                               labels = lable1[i],
                               facing = "downward", niceFacing = TRUE,
                               cex = 0.6,
                               adj = c(-0.2, 0))
                 }
               })

  suppressMessages(for(i in 1:nc){
    circlize::circos.text(x = i-0.4,
                y = 7,
                labels = lable2[i],
                facing = "clockwise", niceFacing = TRUE,
                cex = 0.6,adj = c(0, 0))
  })


  max_height <-max(attr(dend, "height"))
  circlize::circos.track(ylim = c(0, max_height),bg.border = NA,track.height = 0.3,
               panel.fun = function(x, y){
                 circlize::circos.dendrogram(dend = dend,
                                   max_height = max_height)
               })
  circlize::circos.clear()

  lgd <- ComplexHeatmap::Legend(at = c(-5,-2, 0, 2, 5), col_fun = col_fun,
                title_position = "topcenter",title = "Z-score")
  ComplexHeatmap::draw(lgd, x = unit(0.7, "npc"), y = unit(0.7, "npc"))
}

geneplot(data=as.matrix(data_logFC[1:6,]),colors = c("blue", "white", "red"))
geneplot(data=as.matrix(data_logFC[1:6,]),colors = c("skyblue", "white", "red"))
