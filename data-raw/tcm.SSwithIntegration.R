
data <- data_logFC
input <- data_logFC[, 1, drop = FALSE]
tcm.SSwithIntegration(input,data)

tcm.SSwithIntegration <- function(input, data, method = "spearman", threshold=1, cores = 1L){

  stopifnot(is.data.frame(input))

  message(paste0("\n", ">>> Running ", "tcm.SSwithIntegration"))

  if(is.data.frame(input)){
    message(paste0("try to convert profile to genesets by logFC threshold ±",threshold))
    upset <- rownames(input)[ input>=threshold]
    downset <- rownames(input)[ input<=-threshold]
    input2 <- list(upset = upset, downset = downset)

    message("convert done and calculation begining")

    cmap_kk <- suppressMessages(tcm.SSwithCMAP(input = input2, data = data,cores= cores))
    lincs_kk <- suppressMessages(tcm.SSwithLINCS(input2, data_logFC,cores= cores))
    zs_kk <- suppressMessages(tcm.SSwithZS(input2, data_logFC,cores= cores))
    xsum_kk <- suppressMessages(tcm.SSwithXSum(input2,data,cores= cores))
  }

  cor_kk <- suppressMessages(tcm.SSwithCorrelation(input, data_logFC, method = method))
  xcos_kk <- suppressMessages(tcm.SSwithXCos(input,data))
  gcmap_kk <- suppressMessages(tcm.SSwithGCMAP(input, data))

  cmap_kk <- dplyr::arrange(cmap_kk,scaled_score)
  lincs_kk <- dplyr::arrange(lincs_kk,scaled_score)
  zs_kk <- dplyr::arrange(zs_kk,scaled_score)
  cor_kk <- dplyr::arrange(cor_kk,scaled_score)
  gcmap_kk <- dplyr::arrange(gcmap_kk,scaled_score)
  xsum_kk  <- dplyr::arrange(xsum_kk ,scaled_score)
  xcos_kk  <- dplyr::arrange(xcos_kk ,scaled_score)

  C1 <- data.frame(
      name=  cmap_kk$tcm,
      method = "CMAP"
    )  %>% na.omit()
  C2 <- data.frame(
      name= lincs_kk$tcm,
      method = "LINCS"
    )  %>% na.omit()
  C3 <- data.frame(
      name= zs_kk$tcm,
      method = "ZS"
    )  %>% na.omit()
  C4 <- data.frame(
      name= cor_kk$tcm,
      method = "Cor"
    )  %>% na.omit()
  C5 <- data.frame(
      name = gcmap_kk$tcm,
      method = "gCMAP"
    )  %>% na.omit()
  C6 <- data.frame(
      name= xsum_kk$tcm,
      method = "XSUM"
    )  %>% na.omit()

  C7 <- data.frame(
    name= xcos_kk$tcm,
    method = "XCOS"
  )  %>% na.omit()

  CA <- rbind(C1,C2,C3,C4,C5,C6,C7)
  summary <- CA %>% dplyr::group_by(name) %>% dplyr::summarise(Freq=dplyr::n(), method = paste(method, collapse = ", ")) %>% na.omit() %>% dplyr::arrange(desc(Freq))

  star_glist <- list(C1$name,C2$name,C3$name,C4$name,C5$name,C6$name,C7$name)
  star_rra <- RobustRankAggreg::aggregateRanks(star_glist)
  star_rra <- cbind(star_rra, "rra_rank" = seq(1,nrow(star_rra)))

  colnames(star_rra) <- tolower(colnames(star_rra))

  summary <- dplyr::left_join(star_rra, summary, by="name")

  return(summary)

}





















