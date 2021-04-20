library("circlize")
## 首先需要构建数据
set.seed(999)

mat <-
  matrix(rnorm(100 * 10), nrow = 10, ncol = 100)

mat <- as.matrix(data_logFC, ncol = 103)

mat <- mat[1:6, ]
## 构建颜色转变函数,数值将按照线性转变为对应的颜色.
col_fun <- colorRamp2(
  breaks = c(min(mat), median(mat), max(mat)),
  colors = c("skyblue", "white", "red")
)

col_fun <- colorRamp2(
  breaks = c(min(mat), median(mat), max(mat)),
  colors = c("#005522", "#FFFACD", "firebrick3")
)

col_fun <- colorRamp2(
  breaks = c(min(mat), median(mat), max(mat)),
  colors = c("navy", "white", "firebrick3")
)



## 因为我们只有一个热图,因此我们只需要一个sector(扇形)即可,所以factor只需要一个.
factors <- "a"
## 需要对数据进行聚类,从而得到其dendrogram对象,用于后续的绘图
dend <-
  as.dendrogram(hclust(dist(t(mat))))


### 对整体函数做一些设置
circos.par(
  cell.padding = c(0, 0, 0, 0),
  start.degree = 90
)
## 初始化一个circos,注意xlim
circos.initialize(factors, xlim = c(0, 100))
## 使用circos.track创建第一个track,并将热图使用circos.rect函数绘制在该轨道中
circos.track(
  ylim = c(0, 10),
  bg.border = NA,
  track.height = 0.6,
  panel.fun = function(x, y) {
    mat2 <- mat[, order.dendrogram(dend)]
    col_mat <- col_fun(mat2)
    nr <- nrow(mat2)
    nc <- ncol(mat2)
    for (i in 1:nr) {
      circos.rect(
        xleft = 1:nc - 1, ybottom = rep(nr - i, nc),
        xright = 1:nc, ytop = rep(nr - i + 1, nc),
        border = "white",
        col = col_mat[i, ]
      )
    }

    for (i in 1:nc) {
      circos.text(
        x = c(1:100) - 0.5,
        y = 10,
        labels = c(1:100),
        facing = "clockwise", niceFacing = TRUE,
        cex = 1,
        adj = c(0, 0.5)
      )
    }
  }
)
## 然后使用circos.dendrogram将circos.dendrogram画在内部

## 首先获得dendrogram对象的最高值
max_height <-
  max(attr(dend, "height"))

circos.track(
  ylim = c(0, max_height),
  bg.border = NA,
  track.height = 0.3,
  panel.fun = function(x, y) {
    circos.dendrogram(
      dend = dend,
      max_height = max_height
    )
  }
)

circos.clear()











library("circlize")
## 首先需要构建数据
set.seed(999)

mat <-
  matrix(rnorm(100 * 10), nrow = 10, ncol = 100)

mat <- as.matrix(data_logFC, ncol = 103)

mat <- mat[1:6, ]
## 构建颜色转变函数,数值将按照线性转变为对应的颜色.
col_fun <- colorRamp2(
  breaks = c(min(mat), median(mat), max(mat)),
  colors = c("skyblue", "white", "red")
)

col_fun <- colorRamp2(
  breaks = c(min(mat), median(mat), max(mat)),
  colors = c("#005522", "#FFFACD", "firebrick3")
)

col_fun <- colorRamp2(
  breaks = c(min(mat), median(mat), max(mat)),
  colors = c("navy", "white", "firebrick3")
)

## 因为我们只有一个热图,因此我们只需要一个sector(扇形)即可,所以factor只需要一个.
factors <- "a"
## 需要对数据进行聚类,从而得到其dendrogram对象,用于后续的绘图
dend <-
  as.dendrogram(hclust(dist(t(mat))))

## 需要对数据进行聚类,从而得到其dendrogram对象,用于后续的绘图
dend <-
  as.dendrogram(hclust(dist(t(mat))))
## 将数据分为三类,并标为不同的颜色
library(tidyverse)




### 对整体函数做一些设置
circos.par(
  cell.padding = c(0, 0, 0, 0),
  start.degree = 120
)
## 初始化一个circos,注意xlim
circos.initialize(factors, xlim = c(0, 104))

## 创建一个空track
circos.track(
  ylim = c(0, 10),
  bg.border = NA,
  track.height = 0.6
)

col_mat <- col_fun(mat)

nr <- nrow(mat)
nc <- ncol(mat)

for (i in 1:nr) {
  circos.rect(
    xleft = 1:nc - 1,
    ybottom = rep(nr - i, nc),
    xright = 1:nc,
    ytop = rep(nr - i + 1, nc),
    border = "white",
    sector.index = "b",
    track.index = 1,
    col = col_mat[i, ]
  )
}

## 添加文字
for (i in 1:nc) {
  circos.text(
    x = c(1:50) - 0.5,
    y = 10,
    labels = paste("Var", c(1:50), sep = ""),
    facing = "clockwise",
    niceFacing = TRUE,
    cex = 0.8,
    adj = c(0, 0.8),
    sector.index = "b",
    track.index = 1
  )
}

## 添加y坐标
circos.yaxis(
  side = "left", at = c(1:10),
  labels = paste("Sample", 1:10, sep = ""),
  tick = FALSE, sector.index = "b",
  labels.font = 0.8,
  col = "white"
)

## 首先获得dendrogram对象的最高值
max_height <-
  max(attr(dend, "height"))

## 创建新的空track
circos.track(
  ylim = c(0, max_height),
  bg.border = NA,
  track.height = 0.3
)

circos.dendrogram(
  dend = dend,
  max_height = max_height
)
circos.clear()
## 添加legend
library(ComplexHeatmap)

lgd <- Legend(
  at = c(-2, -1, 0, 1, 2), col_fun = col_fun,
  title_position = "topcenter",
  title = "Intensity"
)

draw(lgd,
  y = unit(100, "mm"),
  just = c("top", "center")
)
