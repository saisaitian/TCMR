devtools::load_all()
data <- load_example_dataset()

group_list_all <- stringr::str_split(as.character(data$pdata$title),'_',simplify = T)[,2]
i=1

#31和65,106都是DMSO对照！
# 1. set the group -----------------------------------------------------------

num <- function(i){

  if(i %in% 1:30){
    num=c(2*i-1,2*i,61,62)
  }else if(i %in% 32:64 ){
    num=c(2*i-1,2*i,129,130)
  }else if(i %in% 66:105 ){
    num=c(2*i-1,2*i,211,212)
  }
  num
}

# 2. DE analysis -----------------------------------------------------------

deg.cal <- function(i){

  sel <- num(i)
  exp = data$expr[,sel]
  group_list = group_list_all[sel]

  #limma，获取deg
  {
    design=model.matrix(~factor(group_list,levels=group_list[c(3,1)]) ) #~表示以~后面连接的东西作为因变量,level表示factor的排序，因为我们希望DMSO为0，实验组为1，只有这样limma包才能正确的匹配出来
    fit=limma::lmFit(exp,design)
    fit=limma::eBayes(fit)   #EBayes统计
    deg=limma::topTable(fit,coef=2,number = Inf)
  }

  #对deg添加内容
  {
    #1.加symbol列，把行名变成一列

    deg <- dplyr::mutate(deg,symbol=rownames(deg)) #相当于把列名单独成了一行
    #head(deg)

  }
  mixcache=list(deg)
  names(mixcache)=group_list[1]
  mixcache
}
num(1)
data$expr[,c(1,2,61,62)]
aa=deg.cal(1)

deg.all =lapply(c(1:30,32:64,66:105), function(x)deg.cal(x))

#saveRDS(deg.all,file = './data-raw/deg.all.RDS')

#deg.all <- load_example_deg()
