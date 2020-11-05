## code to prepare `DATASET` dataset goes here
setwd("data-raw/")

# library(GEOquery)
# eSet <- getGEO("GSE85871", GSEMatrix=T, AnnotGPL=FALSE)

exprSet_simple <- exprs(eSet[[1]])
pdata <- pData(eSet[[1]])

exprSet_simple <- data.frame(exprSet_simple, check.rows = F, check.names = F)
ttt <- eSet$GSE85871_series_matrix.txt.gz@featureData@data
exprSet_simple$gsym <- eSet$GSE85871_series_matrix.txt.gz@featureData@data$`Gene Symbol`

exprdf <- exprSet_simple
# 删除没有gene symbol的探针组
exprdf <- exprdf[exprdf$gsym != "", ]
dim(exprdf)

# 有的探针组对应多个基因，用“///”分隔基因名，删掉这样的行
exprdf <- exprdf[!grepl("///", exprdf$gsym), ]
dim(exprdf)

# 有多个探针组对应同一个基因，取均值
exprdf_uniq <- aggregate(. ~ gsym, exprdf, mean)

dim(exprdf_uniq)

# 现在就可以用gene symbol作为行名了
rownames(exprdf_uniq) <- exprdf_uniq$gsym
# 删除gene symbol列
exprdf_uniq <- subset(exprdf_uniq, select = -gsym)
exprdf_uniq

# 保存所有基因的表达矩阵到文件

# save(exprdf_uniq,pdata,file ="GSE85871_simple_pdata.Rdata" )

library(dplyr)
expr <- exprdf_uniq %>% tibble::rownames_to_column("symbol")
# 减少一些磁盘占用
vroom::vroom_write(expr, path = "../inst/extdata/GSE85871_expr.tsv.gz")
vroom::vroom_write(pdata, path = "../inst/extdata/GSE85871_pdata.tsv.gz")
# xx = data.table::fread("../inst/extdata/GSE85871_expr.tsv.gz")
# xx = data.table::fread("../inst/extdata/GSE85871_pdata.tsv.gz")

setwd("../")

load("data-raw/pdata.rda")

pdata %>%
  dplyr::select(c(
    "title", "geo_accession", "platform_id",
    "cell line:ch1", "duration:ch1",
    "group:ch1", "perturbagen:ch1"
  )) %>%
  setNames(c(
    "title", "geo_accession", "platform_id", "cell_line",
    "duration", "group", "perturbagen"
  )) -> pdata
vroom::vroom_write(pdata, path = "inst/extdata/GSE85871_pdata.tsv.gz")

pdata <- data.table::fread(system.file(
  "extdata", "GSE85871_pdata.tsv.gz",
  package = "TCMR", mustWork = TRUE
), data.table = FALSE)

# Handle weird words
pdata$title <- sub("渭", "μ", pdata$title)
pdata$title <- gsub("灏戼㹥", "β", pdata$title)
pdata$perturbagen <- sub("渭", "μ", pdata$perturbagen)
pdata$title <- sub("尾", "β", pdata$title)
pdata$perturbagen <- sub("尾", "β", pdata$perturbagen)

vroom::vroom_write(pdata, path = "inst/extdata/GSE85871_pdata.tsv.gz")

