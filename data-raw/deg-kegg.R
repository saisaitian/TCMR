
devtools::load_all()
library(purrr)
library(clusterProfiler)
library(ggplot2)
data("AnalyzedDEG")
one_report <- load_analyzedDEG(1)
names(one_report)

deg <- subset(one_report,abs(logFC)>1&P.Value<0.05)

getsigpath <- function(data,
                       p.cutoff = 0.05,
                       p.adj.cutoff = NULL,
                       sortby = 'pvalue',
                       n.path = 10){

  stopifnot(is.data.frame(data), !is.null(rownames(data)))
  identifiers = suppressWarnings(bitr(data$identifier, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db"))
  message("The number of gene is: ",paste0(nrow(identifiers), collapse = "  "),' for next analysis')
  kk <- suppressWarnings(enrichKEGG(gene= identifiers$ENTREZID,
                                    organism     = 'hsa',
                                    pvalueCutoff = 1,
                                    qvalueCutoff = 1))

  kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

  df <- as.data.frame(kk)

  if(!is.null(p.cutoff) ){
    df <- df[which(df$pvalue < p.cutoff),]
  }

  if(!is.null(p.adj.cutoff) ){
    df <- df[which(df$p.adjust < p.adj.cutoff),]
  }

  if(!is.null(p.cutoff)& !is.null(p.adj.cutoff)){
    df <- df[which(df$pvalue < p.cutoff & df$p.adjust < p.adj.cutoff),]
  }

  if(nrow(df) > n.path) {
    pathway<- df[1:n.path,]
  } else {
    pathway <- df
  }

  sortdf <- pathway[order(pathway[,sortby],pathway$GeneRatio),]

  sortdf$Description <- factor(sortdf$Description,levels = sortdf$Description)

  sortdf$GeneRatio <- DOSE::parse_ratio(sortdf$GeneRatio)

  return(sortdf)

}

aa=getsigpath(deg,p.cutoff = 0.05,p.adj.cutoff = 0.5,n.path = 10)
aaa=getsigpath(deg,p.cutoff = 0.05,p.adj.cutoff = 0.5,sortby ='p.adjust' ,n.path = 10)

names(aa)

getbarplot <- function(data,
                       color=c("blue","red"),
                       linetype='dashed',
                       linecol="chocolate",
                       pointcol="black",
                       filcol='red',
                       colby='pvalue'
                       ){

  if(!is.element(colby, c("pvalue", "p.adjust"))) {
    stop("the argument of distance should be one of pvalue, p.adjust.")
  }
  if(colby=='pvalue'){
    ylab <- expression('-log'[10]*'(P Value)')
  }
  if(colby=='p.adjust'){
    ylab <- expression('-log'[10]*'(P.adjust)')
  }

  mx <- max(-log10(data[,colby]))
  colbys <- data[,colby]
  ggplot(data=data, mapping=aes(x=Description,
                              y=-log10(colbys),
                              fill=-log10(colbys),
                              group = 1)) +
    geom_bar(stat='identity') +
    scale_fill_continuous(low=color[1], high=color[2])+
    geom_line(aes(x=Description, y=-log10(colbys)),stat="identity",
              linetype = linetype, color = linecol, size = 1.2)+
    geom_point(shape=21, color=pointcol, fill=filcol, size=5)+
    scale_x_discrete(limits=rev(data$Description)) +
    ylim(0,  mx) +
    labs(x='', y= ylab) + coord_flip() +
    geom_text(aes(label=data$Count), hjust=-1, size=5)+
    labs(fill=ylab) +
    theme_bw(15)+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour='white'),
          panel.background = element_blank())  +
    theme(title=element_text(family="myFont",size=12,color="black",
                             face="plain",hjust=0.2,lineheight=0.2),
          axis.title.x=element_text(size=12,face="plain",color="black",hjust=0.5),
          axis.title.y=element_text(size=14,color="black",hjust=0.5,angle=45),
          axis.text.x=element_text(family="myFont",size=8,color="black"),
          axis.text.y=element_text(family="myFont",size=12,face="plain",color="black"))+
    theme(legend.position=c(1,0.2),legend.justification = c(1,0.2))

}

getbarplot(aa)


getbarplot(aa,
           color=c("yellow","red"),
           linetype='dashed',
           linecol="black",
           pointcol="black",
           filcol='red',
           colby='p.adjust')


getdotplot <- function(data,
                       fill= 'pvalue',
                       color=c("blue","red")){

  if(!is.element(fill, c("pvalue", "p.adjust"))) {
    stop("the argument of distance should be one of pvalue, p.adjust.")
  }

  fills <- data[,fill]
  min_x <- min(data$GeneRatio)
  max_x <- max(data$GeneRatio)

  if(fill=='pvalue'){
    collab <- expression('-log'[10]*'(P Value)')
  }
  if(fill=='p.adjust'){
    collab <- expression('-log'[10]*'(P.adjust)')
  }

  ggplot(data,aes(GeneRatio,Description,colour=-log10(fills)))+
    geom_point(aes(size=Count))+
    scale_color_gradientn(colours=color)+
    scale_size_continuous(range = c(2,10))+
    scale_x_continuous(limits = c(min_x,max_x))+
    theme_bw(15)+
    ylab("")+
    labs(size="Count",colour=collab) +
    theme(legend.position=c(1,0.5),legend.justification = c(1,0.5))+
    theme(legend.background = element_blank())+
    theme(legend.key = element_blank())+
    theme(title=element_text(family="myFont",size=12,color="black",
                             face="plain",hjust=0.2,lineheight=0.2),
          axis.title.x=element_text(size=12,face="plain",color="black",hjust=0.5),
          axis.title.y=element_text(size=14,color="black",hjust=0.5,angle=45),
          axis.text.x=element_text(family="myFont",size=8,color="black"),
          axis.text.y=element_text(family="myFont",size=12,face="plain",color="black"))

}


getdotplot(aa,fill= 'p.adjust',
           color=c("blue","red"))


getdotplot(aa,fill= 'p.adjust',
           color=c("yellow","red"))


getbarplot(aa,
           color=c("yellow","red"),
           linetype='dashed',
           linecol="black",
           pointcol="black",
           filcol='red',
           colby='p.adjust')


getbarplot(aa,
           color=c("blue","red"),
           linetype='dashed',
           linecol="black",
           pointcol="black",
           filcol='red',
           colby='pvalue')
