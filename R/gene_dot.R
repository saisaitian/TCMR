
#' Title Dotplot of special genes in special compounds of logFC matrix
#'
#' @param data the special genes expression data of special compounds
#' @param color the color to distinguish between different trends
#'
#' @return A `ggplot` object
#' @export
#' @importFrom reshape2 melt
#' @importFrom ggplot2 aes_
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
#' data=as.matrix(data_logFC[1:6,])
#' data <- data.frame(t(data),check.rows = F,check.names = F)
#' data$compound <- rownames(data)
#' data <- data[3:13,]
#' gene_dot(data)
gene_dot <- function(data,
                     color=c('red','blue')){

  tmp <- reshape2::melt(data = data,id.vars = 'compound')

  tmp$type <- ifelse(tmp$value > 0,'UP','Down')
  tmp$expression <- abs(tmp$value)

  ggplot2::ggplot(tmp,aes(x=variable,y=compound,colour=type))+
    geom_point(aes(size=expression))+
    theme_bw()+
    theme(panel.grid = element_blank(),
          axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
    scale_color_manual(values=c(UP='blue',Down='red'))+
    labs(x=NULL,y=NULL)+
    theme(
      title = element_text(
        size = 12, color = "black",
        face = "plain", hjust = 0.2, lineheight = 0.2
      ),
      axis.title.x = element_text(size = 12, face = "plain", color = "black", hjust = 0.5),
      axis.title.y = element_text(size = 14, color = "black", hjust = 0.5, angle = 45),
      axis.text.x = element_text(size = 8, color = "black"),
      axis.text.y = element_text(size = 12, face = "plain", color = "black")
    )

}









