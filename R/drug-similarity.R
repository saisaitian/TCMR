#' Similarity calculation between two drugs
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
#' tcm.GetDrugSimilarity(smile_a = SMILES[1, ], smile_b = SMILES[2, ], method = "fp", plot = TRUE)
tcm.GetDrugSimilarity <- function(smile_a,
                                  smile_b,
                                  method = "fp",
                                  plot = T) {
  message(paste0("\n", ">>> Calculating ", "tcm.GetDrugSimilarity"))

  sdf <- suppressWarnings(ChemmineR::smiles2sdf(c(smile_a, smile_b)))
  ChemmineR::cid(sdf) <- c("Molecular A", "Molecular B")
  ap <- ChemmineR::sdf2ap(sdf)
  if (method == "ap") {
    index <- ChemmineR::cmp.similarity(ap[1], ap[2])
  }
  if (method == "fp") {
    # suppressWarnings(library(BioMedR))
    fp <- ChemmineR::desc2fp(ap)
    index <- ChemmineR::fpSim(fp[1], fp[2])
    names(index) <- ""
  }
  if (plot == T) {
    ChemmineR::plot(sdf[1:2], regenCoords = TRUE, print = FALSE) # 'print=TRUE' returns SDF summaries
  }

  return(index)
}
