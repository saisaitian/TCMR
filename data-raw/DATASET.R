## code to prepare `DATASET` dataset goes here

library(GEOquery)
eSet <- getGEO("GSE85871", GSEMatrix=T, AnnotGPL=FALSE)

exprSet_simple=exprs(eSet[[1]])
pdata=pData(eSet[[1]])

exprSet_simple <- data.frame(exprSet_simple,check.rows = F,check.names = F)
ttt <- eSet$GSE85871_series_matrix.txt.gz@featureData@data
exprSet_simple$gsym<- eSet$GSE85871_series_matrix.txt.gz@featureData@data$`Gene Symbol`

exprdf <- exprSet_simple
#删除没有gene symbol的探针组
exprdf<-exprdf[exprdf$gsym!="",]
dim(exprdf)

#有的探针组对应多个基因，用“///”分隔基因名，删掉这样的行
exprdf<-exprdf[!grepl("///", exprdf$gsym),]
dim(exprdf)

#有多个探针组对应同一个基因，取均值
exprdf_uniq<-aggregate(.~gsym,exprdf,mean)

dim(exprdf_uniq)

#现在就可以用gene symbol作为行名了
rownames(exprdf_uniq)<-exprdf_uniq$gsym
#删除gene symbol列
exprdf_uniq<-subset(exprdf_uniq,select = -gsym)
exprdf_uniq

#保存所有基因的表达矩阵到文件

save(exprdf_uniq,pdata,file ="GSE85871_simple_pdata.Rdata" )

usethis::use_data(exprdf_uniq)
usethis::use_data(pdata)
