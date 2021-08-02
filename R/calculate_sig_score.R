
#' Calculating signature score using z-score method
#'
#' @param eset transcriptomic data,please make sure a microarray profile or a normalized expression data [FPKM or TPM without log2 transformation is recommended] was provided.
#' @param signature List of gene signatures
#' @param mini_gene_count filter out signatures with genes less than minimal gene in expression set;
#' @return data frame with pdata and signature scores for gene sets; signatures in columns, samples in rows
#' @export
#' @examples There is no example and please refer to vignette.
#'


calculate_sig_score_zscore<-function(eset,
                                     signature,
                                     mini_gene_count){

  message(paste0("\n", ">>> Calculating signature score using z-score method"))

  eset<-scale(eset,center = T,scale = T)
  ###########################
  pdata<-data.frame("Index" = 1:length(colnames(eset)),"ID" = colnames(eset))

  signature<-signature[lapply(signature,function(x) sum(x%in%rownames(eset)==TRUE))>= mini_gene_count]

  ###########################
  #calculating signature score
  goi <- names(signature)
  for (sig in goi) {
    pdata[, sig] <- NA
    genes <- signature[[sig]]
    genes <- genes[genes %in% rownames(eset)]
    tmp <- eset[genes, , drop=FALSE]
    pdata[, sig] <- sigScore(tmp,methods = "zscore")
  }
  if ("TMEscoreA_CIR"%in%goi &"TMEscoreB_CIR" %in% goi) {
    pdata[,"TMEscore_CIR"]<-pdata[,"TMEscoreA_CIR"]-pdata[,"TMEscoreB_CIR"]
  }
  if ("TMEscoreA_plus"%in%goi &"TMEscoreB_plus" %in% goi) {
    pdata[,"TMEscore_plus"]<-pdata[,"TMEscoreA_plus"]-pdata[,"TMEscoreB_plus"]
  }
  pdata<-tibble::as_tibble(pdata)
  return(pdata)
}
###################################################


#' Calculating signature score using ssGSEA method
#' @param eset normalizated  transcriptomic data: normalized (CPM, TPM, RPKM, FPKM, etc.)
#' @param signature List of gene signatures
#' @param mini_gene_count filter out signatures with genes less than minimal gene in expression set; default is 5;
#' the minimal gene count for ssGSEA methods should larger than 5 for the robustness of the calculation
#' @return data frame with pdata and signature scores for gene sets; signatures in columns, samples in rows
#' @export
#' @import tibble
#' @examples There is no example and please refer to vignette.

calculate_sig_score_ssgsea<-function(eset,
                                     signature,
                                     mini_gene_count){

  message(paste0("\n", ">>> Calculating signature score using ssGSEA method"))

  signature<-signature[lapply(signature,function(x) sum(x%in%rownames(eset)==TRUE))>= mini_gene_count]

  if(mini_gene_count<=5) mini_gene_count <- 5
  ############################

  ##############################
  res <- GSVA:: gsva(as.matrix(eset),
                     signature,
                     method="ssgsea",
                     kcdf="Gaussian",
                     min.sz= mini_gene_count,
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

#' Calculating signature score  on a gene expression dataset
#'
#' @param eset normalizaed  transcriptomic data: normalized (CPM, TPM, RPKM, FPKM, etc.)
#' @param signature List of gene signatures;
#' such as `signature_tme`
#' @param method he methods currently supported are `ssgsea`, `zscore`
#' @param mini_gene_count filter out signatures with genes less than minimal gene in expression set;
#' default is zscore funciion
#' @return data frame with pdata and signature scores for gene sets; signatures in columns, samples in rows
#' @export
#' @examples There is no example and please refer to vignette.
#' @references 1. Hänzelmann S, Castelo R, Guinney J (2013). “GSVA: gene set variation analysis for microarray and RNA-Seq data.” BMC Bioinformatics, 14, 7. doi: 10.1186/1471-2105-14-7
#' 2. Mariathasan, S., Turley, S., Nickles, D. et al. TGFβ attenuates tumour response to PD-L1 blockade by contributing to exclusion of T cells. Nature 554, 544–548 (2018).
#'
calculate_sig_score<-function(eset,
                              signature = signature_collection,
                              method = "zscore",
                              mini_gene_count = 3){


  if(class(signature)=="list"){
    signature<-lapply(signature,function(x) na.omit(x))
    signature<-lapply(signature,function(x) as.character(x))
    signature<-lapply(signature,function(x) unique(x))
    signature<-lapply(signature,function(x) x[!x==""])
  }

  ##########################################

  method<-tolower(method)

  if(!method%in%c("zscore","ssgsea")) stop("At present, we only provide two methods to calculate the score: zscore, ssGSEA")
  # run selected method
  res = switch(method,
               ssgsea = calculate_sig_score_ssgsea(eset,
                                                   signature = signature,
                                                   mini_gene_count = mini_gene_count),
               zscore = calculate_sig_score_zscore(eset,
                                                   signature = signature,
                                                   mini_gene_count = mini_gene_count))
  return(res)
}

