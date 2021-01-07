#' Lincs method for Signature Search
#'
#' @inherit ss_cmap
#'
#' @export
#' @importFrom ExperimentHub ExperimentHub
#' @examples
#' data("data_logFC")
#' upset <- rownames(data_logFC)[1:100]
#' downset <- rownames(data_logFC)[400:550]
#' input <- list(upset = upset, downset = downset)
#' lincs_kk <- ss_lincsscore(input, data_logFC)
ss_lincsscore <- function(input, data) {
  stopifnot(all(c("upset", "downset") %in% names(input)), is.list(input))
  data <- as.matrix(data)

  upset <- input$upset
  downset <- input$downset
  if (!is.null(upset)) {
    num <- sum(upset %in% rownames(data))
    if (num <= 10) {
      stop("the up gene less than 10")
    }
    message(paste(
      sum(upset %in% rownames(data)), "/", length(upset),
      "genes in up set share identifiers with reference database"
    ))

    upset <- upset[upset %in% rownames(data)]
  }
  if (!is.null(downset)) {
    num <- sum(downset %in% rownames(data))
    if (num <= 10) {
      stop("the down gene less than 10")
    }
    message(paste(
      sum(downset %in% rownames(data)), "/", length(downset),
      "genes in down set share identifiers with reference database"
    ))
    downset <- downset[downset %in% rownames(data)]
  }
  if (is.null(upset) & is.null(downset)) {
    stop("Both upset and downset share zero identifiers with reference database,
          please make sure that at least one share identifiers!")
  }


  lincsup <- apply(data, 2, function(x) {
    calES(
      sigvec = sort(x, decreasing = TRUE),
      Q = upset
    )
  })

  lincsdown <- apply(data, 2, function(x) {
    calES(
      sigvec = sort(x, decreasing = TRUE),
      Q = downset
    )
  })

  lincsout <- ifelse(sign(lincsup) != sign(lincsdown), (lincsup - lincsdown) / 2, 0)

  eh <- suppressMessages(ExperimentHub::ExperimentHub())
  CSnull <- suppressMessages(eh[["EH3234"]])
  CSnull[CSnull[, "Freq"] == 0, "Freq"] <- 1
  myrounding <- max(nchar(as.character(CSnull[, "WTCS"]))) - 3
  lincs_round <- round(as.numeric(lincsout), myrounding)

  CS_pval <- vapply(lincs_round, function(x) {
    sum(CSnull[abs(CSnull[, "WTCS"]) > abs(x), "Freq"]) / sum(CSnull[, "Freq"])
  }, FUN.VALUE = numeric(1))
  CS_fdr <- stats::p.adjust(CS_pval, "fdr")
  # grouping <- paste(gsub("^.*?_", "", names(lincsout)),
  #                   as.character(ifelse(ESout > 0, "up", "down")), sep="_")

  # lincsout_na <- as.numeric(lincsout)
  # lincsout_na[lincsout_na == 0] <- NA

  # groupmean <- tapply(lincsout_na, grouping, mean, na.rm=TRUE)
  # groupmean[is.na(groupmean)] <- 0
  #
  # groupmean[groupmean==0] <- 10^-12
  #
  # ncs <- as.numeric(lincsout) / abs(groupmean[grouping])
  score <- .S(lincsout)

  result <- data.frame(
    name = names(lincsout),
    trend = as.character(ifelse(lincsout > 0, "up", "down")),
    raw_score = as.numeric(lincsout),
    scaled_score = score,
    Pval = CS_pval,
    FDR = CS_fdr,
    N_upset = length(upset),
    N_downset = length(downset), stringsAsFactors = FALSE
  )
  result <- result[order(abs(result$scaled_score), decreasing = TRUE), ]
  names(result) <- c("tcm", "direction", "raw_score","scaled_score", "pvalue", "fdr", "Nupset",'Ndownset')
  row.names(result) <- NULL
  return(result)
}

calES <- function(sigvec, Q) {
  L <- names(sigvec)
  N <- length(L)
  NH <- length(Q)
  Ns <- N - NH
  hit_index <- as.numeric(L %in% Q)
  miss_index <- 1 - hit_index
  R <- abs(as.numeric(sigvec))
  NR <- sum(R[hit_index == 1])
  if (NR == 0) {
    return(0)
  }
  ESvec <- cumsum((hit_index * R * 1 / NR) - (miss_index * 1 / Ns))
  ES <- ESvec[which.max(abs(ESvec))]
  return(ES)
}
