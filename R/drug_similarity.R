
#' Title The similarity calculation between two drugs
#'
#' @param smile_a the smile of drug a
#' @param smile_b the smile of drug b
#' @param method the method for calculation
#' @param plot a logistical value
#'
#' @return a number
#' @export
#'
#' @examples
#' data("SMILES")
#' drug_similarity(smile_a = SMILES[1,],smile_b = SMILES[2,],method = 'fp',plot = T)

drug_similarity <- function(smile_a,
                            smile_b,
                            method='ap',
                            plot=T){

  message(paste0("\n", ">>> Calculating ", "drug_similarity"))

  sdf <- suppressWarnings(ChemmineR::smiles2sdf(c(smile_a,smile_b)))
  cid(sdf) <- c('Molecular A','Molecular B')
  ap <- sdf2ap(sdf)
  if(method=="ap"){
    index <- cmp.similarity(ap[1],ap[2])
  }
  if(method=="fp"){
  fp <- ChemmineR::desc2fp(ap)
  index <- fpSim(fp[1], fp[2])
  names(index) <- ""
  }
 if(plot==T){

   plot(sdf[1:2], regenCoords=TRUE,print=FALSE) # 'print=TRUE' returns SDF summaries
  }

  return(index)

}













