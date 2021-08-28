
#' #' Calculating drug score using z-score method
#'
#' @param eset transcriptomic data,please make sure a microarray profile or a normalized expression data was provided.
#' @param signature List of gene signatures
#' @return data frame with pdata and signature scores for gene sets; signatures in columns, samples in rows
#' @export

calculate_zscore<-function(eset,
                           signature){

  message(paste0("\n", ">>> Calculating signature score using z-score method"))

  eset<-scale(eset,center = T,scale = T)
  ###########################
  pdata<-data.frame("Index" = 1:length(colnames(eset)),"ID" = colnames(eset))

  #filter signatures
  signature<-signature[lapply(signature,function(x) sum(x%in%rownames(eset)==TRUE))>= 2]

  ##calculating signature score
  sigs <- names(signature)

  for (sig in sigs) {
    pdata[, sig] <- NA
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset)]
    tmp <- eset[genes, , drop=FALSE]
    pdata[, sig] <- colMeans(tmp)

  }
  if ("TMEscoreA_CIR"%in%sigs &"TMEscoreB_CIR" %in% sigs) {
    pdata[,"TMEscore_CIR"]<-pdata[,"TMEscoreA_CIR"]-pdata[,"TMEscoreB_CIR"]
  }
  if ("TMEscoreA_plus"%in%sigs &"TMEscoreB_plus" %in% sigs) {
    pdata[,"TMEscore_plus"]<-pdata[,"TMEscoreA_plus"]-pdata[,"TMEscoreB_plus"]
  }
  pdata<-tibble::as_tibble(pdata)
  return(pdata)
}

###################################################
#' Calculating drug score using ssGSEA method
#'
#' @param eset normalizated  transcriptomic data: normalized (CPM, TPM, RPKM, FPKM, etc.)
#' @param signature List of gene signatures
#' @return data frame with pdata and signature scores for gene sets; signatures in columns, samples in rows
#' @export
#' @import tibble

calculate_ssgsea<-function(eset,
                           signature){

  message(paste0("\n", ">>> Calculating signature score using ssGSEA method"))

  #filter signatures
  signature<-signature[lapply(signature,function(x) sum(x%in%rownames(eset)==TRUE))>= 5]

  ############################
  res <- GSVA::gsva(as.matrix(eset),
                     signature,
                     method="ssgsea",
                     kcdf="Gaussian",
                     ssgsea.norm=T)

  res <- as.data.frame(t(res))
  res <- tibble::rownames_to_column(res,var = "ID")

  if ("TMEscoreA_CIR"%in%colnames(res) & "TMEscoreB_CIR"%in%colnames(res)) {
    res[,"TMEscore_CIR"]<-res[,"TMEscoreA_CIR"] - res[,"TMEscoreB_CIR"]
  }

  if ("TMEscoreA_plus"%in%colnames(res) & "TMEscoreB_plus"%in%colnames(res)) {
    res[,"TMEscore_plus"]<-res[,"TMEscoreA_plus"] - res[,"TMEscoreB_plus"]
  }

  res<-tibble::as_tibble(res)
  return(res)
}
###################################################

#' Calculating drug score  on a gene expression dataset
#'
#' @param eset normalizaed  transcriptomic data: normalized (CPM, TPM, RPKM, FPKM, etc.)
#' @param signature List of gene signatures;
#' such as `signature_tme`
#' @param method the methods currently supported are `ssgsea`, `zscore`
#' default is zscore funciion
#' @importFrom stats cor
#' @importFrom stats na.omit
#' @importFrom stats prcomp
#' @return data frame with pdata and signature scores for gene sets; signatures in columns, samples in rows
#' @export
#' @references 1. Hänzelmann S, Castelo R, Guinney J (2013). “GSVA: gene set variation analysis for microarray and RNA-Seq data.” BMC Bioinformatics, 14, 7. doi: 10.1186/1471-2105-14-7
#' 2. Mariathasan, S., Turley, S., Nickles, D. et al. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature 554, 544–548 (2018).
#'
drug_score<-function(eset,
                     signature = signature_collection,
                     method = "zscore"){

  if(class(signature)=="list"){
    signature<-lapply(signature,function(x) stats::na.omit(x))
    signature<-lapply(signature,function(x) as.character(x))
    signature<-lapply(signature,function(x) unique(x))
    signature<-lapply(signature,function(x) x[!x==""])
  }

  ##########################################

  method<-tolower(method)

  if(!method%in%c("zscore","ssgsea")) stop("At present, we only provide two methods to calculate the score: zscore, ssGSEA")

  # run selected method
  res = switch(method,
               ssgsea = calculate_ssgsea(eset,
                                         signature = signature),
               zscore = calculate_zscore(eset,
                                         signature = signature))
  return(res)
}
