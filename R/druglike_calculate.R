
#' Title druglike analysis
#'
#' @param smiles SMILES of compounds
#'
#' @return a data frame
#' @export
#'
#' @examples
#' data("SMILES")
#' druglike_calculate(SMILES)

druglike_calculate <- function(smiles){
  message(paste0("\n", ">>> Calculating ", "Lipinski Rule of Five - Druglikeness"))
  sdf <- apply(smiles, 1, ChemmineR::smiles2sdf) #1 = rows, convert smiles to sdfs
  openbabel <- lapply(sdf, ChemmineR::propOB)
  hba1 <- matrix(0L, nrow = length(smiles[,1]))
  formula <- hba1
  hbd <- hba1
  logp <- hba1
  mw <- hba1
  TPSA <- hba1
  for (i in 1:length(openbabel)){
    ob <- unlist(openbabel[i])
    formula[i] <- ob[3]
    hba1[i] <- ob[6]
    hbd[i] <- ob[8]
    logp[i] <- ob[9]
    mw[i] <- ob[11]
    TPSA[i] <- ob[13]
  }
  stats <- cbind(smiles,formula, hba1, hbd, logp, mw,TPSA)
  colnames(stats) <- c("SMILES",'Formula', "Hydrogen Bond Donors", "Hydrogen Bond Acceptors", "Octanol/Water Partition Coefficient (logP)", "Molecular Weight",'TPSA')
  return(stats)
}

