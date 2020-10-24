
options(stringsAsFactors = F)
library(GEOquery)
library(stringr)
library(limma)

#第一步：获取GEO数据
load("GSE85871_simple_pdata.Rdata")
#第二步：处理数据，获取想要的信息

pdata[c(61,62,129,130,211,212),1]

#第三步：pd里面能找到里面的分组信息，生成一个group_list向量。

group_list_all= str_split(as.character(pdata$title),'_',simplify = T)[,2]  #如何进行分组呢？
head(group_list_all)

#for循环前运行

mixdeg<-list(name="deg")


#31和65,106都是DMSO对照！
deg.cal <- function(i){

  if(i %in% 1:30){
    sel=c(2*i-1,2*i,61,62)
  }else if(i %in% 32:64 ){
    sel=c(2*i-1,2*i,129,130)
  }else if(i %in% 66:105 ){
    sel=c(2*i-1,2*i,211,212)
  }

  exp = exprdf_uniq[,sel]
  group_list = group_list_all[sel]

  #limma，获取deg
  {
    design=model.matrix(~factor(group_list,levels=group_list[c(3,1)]) ) #~表示以~后面连接的东西作为因变量,level表示factor的排序，因为我们希望DMSO为0，实验组为1，只有这样limma包才能正确的匹配出来
    fit=lmFit(exp,design)
    fit=eBayes(fit)   #EBayes统计
    deg=topTable(fit,coef=2,number = Inf)
  }

  #对deg添加内容
  {
    #1.加symbol列，把行名变成一列

    deg <- dplyr::mutate(deg,symbol=rownames(deg)) #相当于把列名单独成了一行
    #head(deg)

  }

  mixcache=list(deg)
  names(mixcache)=paste("deg",i,group_list[1],sep="_")
  mixcache
}
deg.cal(1)

deg.all =lapply(c(1:30,32:64,66:105), function(x)deg.cal(x))


