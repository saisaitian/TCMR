#' A Function for Differential Expression Analysis Between two groups
#'
#' @param data expression dataset
#' @param group two groups
#' @param level two groups level
#' @param method the method for Differential Expression Analysis
#'
#' @return a `data.frame`.
#' @export
#'
#' @examples
#' data=load_example_dataset()
#' expr=data$expr[,c(1,2,61,62)]
#' group_list=group_list_all[c(1,2,61,62)]
#' tmp=deg.cal(expr,group=group_list,level=group_list[c(3,1)],method='limma')
#'
deg.cal <- function(data= expr,group,level, method='limma'){
  if(method == 'limma'){
    message("Now method is ", method, "...")
    design = model.matrix(~factor(group,levels=level) )
    fit = limma::lmFit(exp,design)
    fit = limma::eBayes(fit)
    deg = limma::topTable(fit,coef=2,number = Inf)
    deg = dplyr::mutate(symbol=rownames(deg),deg)
    deg = dplyr::select(deg,'symbol', dplyr::everything())
    message("Done")
  }
  return(deg)
}


# batchcal <- function(i,data=data){
#   group_list_all <- stringr::str_split(as.character(data$pdata$title), "_", simplify = T)[, 2]
#   if (i %in% 1:30) {
#     sel <- c(2 * i - 1, 2 * i, 61, 62)
#   } else if (i %in% 32:64) {
#     sel <- c(2 * i - 1, 2 * i, 129, 130)
#   } else if (i %in% 66:105) {
#     sel <- c(2 * i - 1, 2 * i, 211, 212)
#   }
#   exp <- data$expr[, sel]
#   group_list <- group_list_all[sel]
#   message('Now calculate is ',group_list[1],'  VS  ',group_list[3]  )
#   result = deg.cal(exp,group=group_list,level = group_list[c(3,1)],method = 'limma')
#   return(result)
# }

#aa=batchcal(1,data = data)
