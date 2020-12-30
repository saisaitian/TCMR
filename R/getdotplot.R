#' Dotplot of enrichResult
#'
#' @param data enrichResult object which is `data.frame` format.
#' @param fill dot color
#' @param color color to show
#'
#' @return ggplot object
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_color_gradientn
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @examples
#' data("Analyzedsigpath")
#' one_report <- load_analyzedsigpath(1)
#' getdotplot(one_report)
#'
#'
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
    theme(title=element_text(size=12,color="black",
                             face="plain",hjust=0.2,lineheight=0.2),
          axis.title.x=element_text(size=12,face="plain",color="black",hjust=0.5),
          axis.title.y=element_text(size=14,color="black",hjust=0.5,angle=45),
          axis.text.x=element_text(size=8,color="black"),
          axis.text.y=element_text(size=12,face="plain",color="black"))

}
