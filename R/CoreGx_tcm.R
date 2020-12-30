
#' Title CoreGx method
#'
#' @param input a data.frame contains logFC value, which rownames should be geme symbols
#' @param data a data.frame contains logFC value data from 103 compounds
#'
#' @return a data.frame
#' @export
#'
#' @examples
#' data("data_logFC")
#' query2 <- data_logFC[1:60,1,drop=FALSE]
#' txp <- CoreGx_tcm(query2,data_logFC[1:3])
CoreGx_tcm <- function(input,data) {

  if(is(input, 'data.frame')){
    num <- sum(rownames(input) %in% rownames(data))
    if(num<=10){
      stop("the commom gene less than 10")
    }
    message(paste( num , "/", length(rownames(data)),
                  "genes in input share identifiers with reference database"))
    if(is.null(input)){
      stop(" Input is NULL !")
    }
  }else{
    stop(" Input is not data.frame !")
  }
  res_list <- NULL
  common <- intersect(rownames(input), rownames(data))
  input2 <- input[common,,drop=F]

  for (i in 1:ncol(data)) {
    tt <- suppressWarnings(CoreGx::connectivityScore(data[,i,drop=F],input2, method = 'fgsea', nperm=1e4,nthread = 1))
    direction=tt[1]
    direction[direction>=0]="up"
    direction[direction<0]="down"
    res <- data.frame(tcm=colnames(data)[i],
                      direction=direction,
                      raw_score = tt[1],
                      p=tt[2],
                      Nset=length(common),
                      stringsAsFactors = FALSE)
    res_list <- c(res_list, list(res))
    print(i)
  }
  result <- do.call(rbind, res_list)
  result$fdr <- stats::p.adjust(result$p, "fdr")
  result$scaled_score <- .S(result$raw_score)
  result <- result[order(abs(result$scaled_score), decreasing=TRUE), ]
  result<- result[, c("tcm", "direction", "raw_score",
                      "scaled_score","p","fdr",'Nset')]
  rownames(result) <- NULL
  return(result)

}


