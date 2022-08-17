#' CMAP method for Signature Search
#'
#' @param input A list containing UP (in `upset` element) and DOWN (in `downset` element) genes.
#' @param data A `data.frame` containing `logFC` value data from 103 compounds.
#' @param cores Number of cores to run the task
#'
#' @return A data.frame.
#' @export
#' @examples
#' data("data_logFC")
#' upset <- rownames(data_logFC)[1:100]
#' downset <- rownames(data_logFC)[400:550]
#' input <- list(upset = upset, downset = downset)
#' cmap_kk <- tcm.SSwithCMAP(input = input, data = data_logFC)
tcm.SSwithCMAP <- function(input, data, cores = 1L) {
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
      "up_genes in up set share identifiers with reference database"
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
      "down_genes in down set share identifiers with reference database"
    ))
    downset <- downset[downset %in% rownames(data)]
  }
  if (is.null(upset) & is.null(downset)) {
    stop("Both upset and downset share zero identifiers with reference database,
          please make sure that at least one share identifiers!")
  }

  message(paste0("\n", ">>> Running ", "tcm.SSwithCMAP"))
  rankLup <- parallel::mclapply(
    colnames(data),
    function(x) sort(rank(-1 * data[, x])[upset]),
    mc.cores = cores
  )
  rankLdown <- parallel::mclapply(
    colnames(data),
    function(x) sort(rank(-1 * data[, x])[downset]),
    mc.cores = cores
  )
  raw.score <- vapply(seq_along(rankLup), function(x) {
    .s(rankLup[[x]], rankLdown[[x]], n = nrow(data))
  },
  FUN.VALUE = numeric(1)
  )
  score <- .S(raw.score)

  library(furrr)
  plan(multisession, workers = cores) # availableCores()-1
  system.time(
    result <- furrr::future_map(1:1000, function(i) {

      ## Prepare the random query signatures
      upset <- sample(rownames(data), size = length(upset))
      downset <- sample(rownames(data), size = length(downset))
      ## Compute the random scores for each sample in the reference lists
      rankLup <- parallel::mclapply(
        colnames(data),
        function(x) sort(rank(-1 * data[, x])[upset]),
        mc.cores = cores
      )
      rankLdown <- parallel::mclapply(
        colnames(data),
        function(x) sort(rank(-1 * data[, x])[downset]),
        mc.cores = cores
      )
      raw.score <- vapply(seq_along(rankLup), function(x) {
        .s(rankLup[[x]], rankLdown[[x]], n = nrow(data))
      },
      FUN.VALUE = numeric(1)
      )
      raw.score
    }, .progress = TRUE)
  )
  permuteScore <- do.call(cbind,result)
  permuteScore[is.na(permuteScore)] <- 0
  ## Compute the p-value based on bootstrap method
  pValue <- rowSums(abs(permuteScore) >= abs(raw.score)) / 1000
  ## Compute the adjusted p-value. The adjusting method can be reseted
  ## (Refer to p.adjust()).
  pAdjust <- stats::p.adjust(pValue, "fdr")
  result <- data.frame(
    set = colnames(data),
    trend = ifelse(score > 0.3, 'up', ifelse(score > -0.3, 'none', 'down')),
    raw_score = raw.score,
    scaled_score = score,
    Pval = pValue,
    FDR = pAdjust,
    stringsAsFactors = FALSE
  )
  result <- result[order(abs(result$scaled_score), decreasing = TRUE), ]
  names(result) <- c("tcm", "direction", "raw_score", "scaled_score", "pvalue", "fdr")
  rownames(result) <- NULL
  return(result)
}

## Function to compute a and b
.ks <- function(V, n) {
  t <- length(V)
  if (t == 0) {
    return(0)
  } else {
    if (is.unsorted(V)) V <- sort(V)
    d <- seq_len(t) / t - V / n
    a <- max(d)
    b <- -min(d) + 1 / t
    ifelse(a > b, a, -b)
  }
}


.s <- function(V_up, V_down, n) {
  ks_up <- .ks(V_up, n)
  ks_down <- .ks(V_down, n)
  ifelse(sign(ks_up) == sign(ks_down), 0, ks_up - ks_down)
}

# Function to scale scores
.S <- function(scores) {
  p <- max(scores)
  q <- min(scores)
  ifelse(scores == 0, 0, ifelse(scores > 0, scores / p, -scores / q))
}


#' Cor method for for Signature Search
#'
#' @param input A `data.frame` containing `logFC` value, which `rownames` should be gene symbols.
#' @param data A `data.frame` contains `logFC` value data from 103 compounds.
#' @param method A string specifying correlation calculation method, should be one of "pearson", "spearman", "kendall".
#'
#' @return A `data.frame`.
#' @export
#'
#' @examples
#' data("data_logFC")
#' query2 <- data_logFC[1:60, 1, drop = FALSE]
#' aa <- tcm.SSwithCorrelation(query2, data_logFC, method = "pearson")
#' aaa <- tcm.SSwithCorrelation(query2, data_logFC, method = "spearman")
#' aaaa <- tcm.SSwithCorrelation(query2, data_logFC, method = "kendall")
tcm.SSwithCorrelation <- function(input, data, method = c("pearson", "spearman", "kendall")) {
  stopifnot(is.data.frame(input))
  method <- match.arg(method)
  num <- sum(rownames(input) %in% rownames(data))
  if (num <= 10) {
    stop("the commom gene less than 10")
  }

  message(paste(
    num, "/", length(rownames(data)),
    "genes in input share identifiers with reference database"
  ))
  if (is.null(input)) {
    stop(" Input is NULL !")
  }
  message(paste0("\n", ">>> Running ", "tcm.SSwithCorrelation"))
  res_list <- NULL
  common <- intersect(rownames(input), rownames(data))
  input2 <- input[common, , drop = F]
  data2 <- data[common, ]

  for (i in 1:ncol(data2)) {
    tt <- stats::cor.test(input2[, 1], data2[, i], method = method)
    p <- tt$p.value
    cor <- tt$estimate
    res <- data.frame(
      tcm = colnames(data2)[i],
      raw_score = cor,
      p = p,
      Nset = length(common),
      stringsAsFactors = FALSE
    )
    res_list <- c(res_list, list(res))
  }
  result <- do.call(rbind, res_list)
  result$fdr <- stats::p.adjust(result$p, "fdr")
  result$scaled_score <- .S(result$raw_score)
  result$direction = ifelse(scaled_score > 0.3, 'up', ifelse(scaled_score > -0.3, 'none', 'down'))
  result <- result[order(abs(result$scaled_score), decreasing = TRUE), ]
  result <- result[, c(
    "tcm", "direction", "raw_score",
    "scaled_score", "p", "fdr"
  )]
  names(result) <- c("tcm", "direction", "raw_score", "scaled_score", "pvalue", "fdr")
  rownames(result) <- NULL
  return(result)
}

#' gcmap Method for Signature Search
#'
#' @inherit tcm.SSwithCorrelation
#' @inheritParams tcm.SSwithCMAP
#' @param higher A cutoff value. If `logFC` larger than or equal to `higher` will be included
#' in the gene set with  `+1`.
#' @param lower A cutoff value. If `logFC` smaller than or equal to `higher` will be included
#' in the gene set to `-1`, others should be set to `0`.
#'
#' @export
#' @examples
#' data("data_logFC")
#' query2 <- data_logFC[1:60, 1, drop = FALSE]
#' gcmap_kk <- tcm.SSwithGCMAP(input = query2, data = data_logFC)
tcm.SSwithGCMAP <- function(input, data, higher = 1, lower = -1, cores = 1L) {
  stopifnot(is.data.frame(input))
  num <- sum(rownames(input) %in% rownames(data))
  if (num <= 10) {
    stop("the commom gene less than 10")
  }

  message(paste(
    num, "/", length(rownames(data)),
    "genes in input share identifiers with reference database"
  ))

  message(paste0("\n", ">>> Running ", "tcm.SSwithGCMAP"))
  data2 <- ifelse(data > higher, 1, ifelse(data > -lower, 0, -1))
  ## subset objects to shared genes
  matched.features <- match(rownames(input), rownames(data2))
  matched.sets <- data2[stats::na.omit(matched.features), ]

  ## extract scores for each gene set
  sets.up <- parallel::mclapply(
    seq(ncol(matched.sets)),
    function(x) which(matched.sets[, x] == 1),
    mc.cores = cores
  )

  sets.down <- parallel::mclapply(
    seq(ncol(matched.sets)),
    function(x) which(matched.sets[, x] == -1),
    mc.cores = cores
  )

  ## transform experiment to (reverse) ranks
  rank.matrix <- apply(input, 2, function(x) {
    length(x) - rank(x) + 1
  })

  ## calculate connectivity score
  raw.score <- apply(rank.matrix, 2, function(r) {
    vapply(seq_along(sets.up), function(n) {
      .s(r[sets.up[[n]]], r[sets.down[[n]]], length(r))
    }, FUN.VALUE = numeric(1))
  })

  raw.score <- matrix(raw.score, ncol = ncol(input))

  library(furrr)
  plan(multisession, workers = cores)
  system.time(
    result <- furrr::future_map(1:1000, function(i) {

      sets <- sample(rownames(data), size = nrow(input))

      matched.features <- match(sets, rownames(data2))

      matched.sets <- data2[stats::na.omit(matched.features), ]

      ## extract scores for each gene set
      sets.up <- parallel::mclapply(
        seq(ncol(matched.sets)),
        function(x) which(matched.sets[, x] == 1),
        mc.cores = cores
      )

      sets.down <- parallel::mclapply(
        seq(ncol(matched.sets)),
        function(x) which(matched.sets[, x] == -1),
        mc.cores = cores
      )
      ## Compute the random scores for each sample in the reference lists
      ## calculate connectivity score
      raw.score <- apply(rank.matrix, 2, function(r) {
        vapply(seq_along(sets.up), function(n) {
          .s(r[sets.up[[n]]], r[sets.down[[n]]], length(r))
        }, FUN.VALUE = numeric(1))
      })

      raw.score <- matrix(raw.score, ncol = ncol(input))
      raw.score

    }, .progress = TRUE)
  )

  permuteScore <- as.matrix(do.call(cbind,result))

  permuteScore[is.na(permuteScore)] <- 0
  ## Compute the p-value based on bootstrap method
  pValue <- rowSums(abs(permuteScore) >= abs(raw.score[,1])) / 1000
  ## Compute the adjusted p-value. The adjusting method can be reseted
  ## (Refer to p.adjust()).
  pAdjust <- stats::p.adjust(pValue, "fdr")

  scaled_score <- .S(raw.score[,1])

  ## store results
  results <- data.frame(
    tcm = colnames(data2),
    direction = ifelse(scaled_score > 0.3, 'up', ifelse(scaled_score > -0.3, 'none', 'down')),
    raw_score = raw.score[, 1],
    Pval = pValue,
    FDR = pAdjust
  )


  ## Apply scaling of scores to full data set
  results[, "scaled_score"] <- .S(results$raw_score)

  results <- results[order(abs(results$scaled_score), decreasing = TRUE), ]
  results <- results[, c(
    "tcm", "direction", "raw_score",
    "scaled_score", "Pval", "FDR"
  )]

  results <- results[order(abs(results$scaled_score), decreasing = TRUE), ]
  names(results) <- c("tcm", "direction", "raw_score", "scaled_score", "pvalue", "fdr")

  return(results)
}

#' Lincs method for Signature Search
#'
#' @inherit tcm.SSwithCMAP
#' @export
#' @examples
#' data("data_logFC")
#' upset <- rownames(data_logFC)[1:100]
#' downset <- rownames(data_logFC)[400:550]
#' input <- list(upset = upset, downset = downset)
#' lincs_kk <- tcm.SSwithLINCS(input, data_logFC)

tcm.SSwithLINCS <- function(input, data) {
  stopifnot(all(c("upset", "downset") %in% names(input)), is.list(input))
  data <- as.matrix(data)
  upset <- input$upset
  downset <- input$downset

  if (is.null(colnames(data)) || is.null(rownames(data))) {
    stop("Warning: data should have both rownames and colnames!")
  }

  if (!is.null(upset)) {
    num <- sum(upset %in% rownames(data))
    if (num <= 10) {
      stop("the up gene less than 10")
    }

    message(paste(
      sum(upset %in% rownames(data)), "/", length(rownames(data)),
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
      sum(downset %in% rownames(data)), "/", length(rownames(data)),
      "genes in down set share identifiers with reference database"
    ))
    downset <- downset[downset %in% rownames(data)]
  }
  if (is.null(upset) & is.null(downset)) {
    stop("Both upset and downset share zero identifiers with reference database,
          please make sure that at least one share identifiers!")
  }


  if (!is.null(upset)) {
    num <- sum(union(upset,downset) %in% rownames(data))
    if (num <= 10) {
      stop("the up gene less than 10")
    }
  }
  message(paste0("\n", ">>> Running ", "tcm.SSwithLINCS"))
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
  score <- .S(lincsout)

  library(furrr)
  plan(multisession, workers = cores) # availableCores()-1
  system.time(
    result <- furrr::future_map(1:1000, function(i) {

      ## Prepare the random query signatures
      upset <- sample(rownames(data), size = length(upset))
      downset <- sample(rownames(data), size = length(downset))
      ## Compute the random scores for each sample in the reference lists
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
      as.vector(lincsout)


    }, .progress = TRUE)
  )

  permuteScore <- do.call(cbind,result)

  permuteScore[is.na(permuteScore)] <- 0

  ## Compute the p-value based on bootstrap method
  pValue <- rowSums(abs(permuteScore) >= abs(lincsout)) / 1000
  ## Compute the adjusted p-value. The adjusting method can be reseted
  pAdjust <- stats::p.adjust(pValue, "fdr")

  result <- data.frame(
    name = names(lincsout),
    trend = ifelse(score > 0.3, 'up', ifelse(score > -0.3, 'none', 'down')),
    raw_score = as.numeric(lincsout),
    scaled_score = score,
    Pval = pValue,
    FDR = pAdjust,
    N_upset = length(upset),
    stringsAsFactors = FALSE
  )
  result <- result[order(abs(result$scaled_score), decreasing = TRUE), ]
  names(result) <- c("tcm", "direction", "raw_score", "scaled_score", "pvalue", "fdr")
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


#' ZS method for Signature Search
#'
#' @inherit tcm.SSwithCMAP
#'
#' @export
#' @examples
#' data("data_logFC")
#' upset <- rownames(data_logFC)[1:100]
#' downset <- rownames(data_logFC)[400:550]
#' input <- list(upset = upset, downset = downset)
#' lincs_kk <- tcm.SSwithZS(input, data_logFC)

tcm.SSwithZS <- function(input, data, cores = 1L) {

  stopifnot(all(c("upset", "downset") %in% names(input)), is.list(input))
  upset <- input$upset
  downset <- input$downset
  if (is.data.frame(data)) {data <- as.matrix(data)}

  if (is.null(colnames(data)) || is.null(rownames(data))) {
    stop("Warning: data should have both rownames and colnames!")
  }

  if (!is.null(upset)) {
    num <- sum(upset %in% rownames(data))
    if (num <= 10) {
      stop("the up gene less than 10")
    }
    message(paste(
      sum(upset %in% rownames(data)), "/", length(rownames(data)),
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
      sum(downset %in% rownames(data)), "/", length(rownames(data)),
      "genes in down set share identifiers with reference database"
    ))
    downset <- downset[downset %in% rownames(data)]
  }
  if (is.null(upset) & is.null(downset)) {
    stop("Both upset and downset share zero identifiers with reference database,
          please make sure that at least one share identifiers!")
  }


  if (!is.null(upset)) {
    num <- sum(union(upset,downset) %in% rownames(data))
    if (num <= 10) {
      stop("the up gene less than 10")
    }
  }
  message(paste(
    num, "/", length(rownames(data)),
    "genes in input share identifiers with reference database"
  ))

  message(paste0("\n", ">>> Running ", "tcm.SSwithZS"))

  ## Convert the gene expression matrix to ranked list
  matrixToRankedList <- function(refMatrix) {
    ## Allocate memory for the refList
    refList <- vector("list", ncol(refMatrix))
    for(i in 1:ncol(refMatrix)) {
      ## Sort the reference gene lists based on the absolute value
      refSort <- refMatrix[order(abs(refMatrix[, i]), decreasing=TRUE), i]
      ## Make ranks for the reference gene lists with sign
      ## if two variables are equal, their ranks will be averaged
      refRank <- rank(abs(refSort)) * sign(refSort)
      refList[[i]] <- refRank
    }
    return(refList)
  }

  ## The core part for computing the zhangscore
  computeScore <- function(refRank, queryRank) {
    if(length(intersect(names(refRank), names(queryRank))) > 0) {
      ## Compute the maximal theoretical score
      maxTheoreticalScore <- sum(abs(refRank)[1:length(queryRank)] *
                                   abs(queryRank))
      ## The final score
      score <- sum(queryRank * refRank[names(queryRank)], na.rm=TRUE)/
        maxTheoreticalScore
    }
    else {
      score <- NA
    }
    return(score)
  }

  ## Prepare the ranked reference lists
  refList <- matrixToRankedList(data)
  ## Prepare the query signatures
  queryVector <- c(rep(1, length(upset)), rep(-1, length(downset)))
  names(queryVector) <- c(upset,downset)
  ## Compute the scores for each sample in the reference lists. mcCore is the
  ## number of cores to use for parallel computing. Set it based on your computer.
  raw.score <- parallel::mclapply(refList, computeScore, queryRank = queryVector,
                                  mc.cores = cores)
  raw.score <- as.vector(do.call(rbind, raw.score))

  library(furrr)
  plan(multisession, workers = cores) # availableCores()-1
  system.time(
    result <- furrr::future_map(1:1000, function(i) {

      bootSample <- sample(c(-1,1), replace = TRUE,
                           size = length(upset) + length(downset))
      names(bootSample) <- sample(rownames(data), replace = FALSE,
                                  size = length(upset) + length(downset))

      bootScore <- parallel::mclapply(refList, computeScore, queryRank = bootSample,
                                      mc.cores = cores)
    }, .progress = TRUE)
  )

  permuteScore <- do.call(cbind,result)

  permuteScore[is.na(permuteScore)] <- 0

  ## Compute the p-value based on bootstrap method
  pValue <- rowSums(abs(permuteScore) >= abs(score)) / 1000
  ## Compute the adjusted p-value. The adjusting method can be reseted
  pAdjust <- stats::p.adjust(pValue, "fdr")
  score <- .S(raw.score)
  results <- data.frame(
    tcm = colnames(data),
    direction = ifelse(score > 0.3, 'up', ifelse(score > -0.3, 'none', 'down')),
    raw_score = raw.score,
    scaled_score = score,
    Pval = pValue,
    FDR = pAdjust,
    stringsAsFactors = FALSE
  )

  results <- results[order(abs(results$scaled_score), decreasing = TRUE), ]
  results <- results[, c(
    "tcm", "direction", "raw_score",
    "scaled_score", "Pval", "FDR"
  )]

  names(results) <- c("tcm", "direction", "raw_score", "scaled_score", "pvalue", "fdr")

  return(results)
}


#' tcm.SSwithXCos method for Signature Search
#'
#' @param input A `data.frame` containing `logFC` value, which `rownames` should be gene symbols.
#' @param data A `data.frame` contains `logFC` value data from 103 compounds.
#' @param cores Number of cores to run the task
#'
#' @return A data.frame.
#' @export
#'
#' @examples

tcm.SSwithXCos <- function(input, data, cores = 1L) {

  stopifnot(is.data.frame(input))

  num <- sum(rownames(input) %in% rownames(data))
  if (num <= 10) {
    stop("the commom gene less than 10")
  }
  res_list <- NULL
  common <- intersect(rownames(input), rownames(data))
  input2 <- input[common, , drop = F]
  input3 <- input2[,1]
  names(input3) <- rownames(input2)

  if (is.data.frame(data)) {data <- as.matrix(data)}

  if (is.null(colnames(data)) || is.null(rownames(data))) {
    stop("Warning: data must have both rownames and colnames!")
  }

  message(paste(
    num, "/", length(rownames(data)),
    "genes in input share identifiers with reference database"
  ))

  message(paste0("\n", ">>> Running ", "tcm.SSwithXCos"))
  # matrixToRankedList(data)
  ## Convert the gene expression matrix to ranked list.
  matrixToRankedList <- function(data) {
    ## Allocate memory for the refList
    refList <- vector("list", ncol(data))
    for(i in 1:ncol(data)) {
      ## Sort the reference gene lists based on fold change of gene expression.
      refList[[i]] <- c(head(data[order(data[, i], decreasing=TRUE), i],
                             n = 500),
                        tail(data[order(data[, i], decreasing=TRUE), i],
                             n = 500))}
    return(refList)
  }
  ## The core part for computing the XCos score
  XCos <- function(refList, query) {
    reservedRef <- refList[match(intersect(names(refList), names(query)),
                                 names(refList))]
    reservedRef[order(names(reservedRef))]
    reservedQuery <- query[match(intersect(names(refList), names(query)),
                                 names(query))]
    reservedQuery[order(names(reservedQuery))]
    ## Compute the cosine similarity
    if (length(reservedRef) == 0) {
      return(NA)
    }
    else {
      return((crossprod(reservedRef, reservedQuery) /
                sqrt(crossprod(reservedRef) * crossprod(reservedQuery)))[1, 1])
    }
  }

  ## Prepare the ranked reference lists
  refList <- matrixToRankedList(data)
  ## Compute the scores for each sample in the reference lists. mcCore is the
  ## number of cores to use for parallel computing. Set it based on your computer.
  raw.score <- parallel::mclapply(refList, XCos, query = input3, mc.cores = cores)
  raw.score <- as.vector(do.call(rbind, raw.score))
  ## Allocate memory for the permuteScore that are used to compute the p-value.
  ## The permuteNum can be reseted. Notice large permuteNum means low speed.

  library(furrr)
  plan(multisession, workers = 16L) # availableCores()-1
  system.time(
    result <- furrr::future_map(1:1000, function(i) {

      names(input3) <- sample(rownames(data), size = length(input3))
      ## Compute the random scores for each sample in the reference lists
      bootScore <- parallel::mclapply(refList, XCos, query = input3, mc.cores = cores)
    }, .progress = TRUE)
  )

  permuteScore <- do.call(cbind,result)

  permuteScore[is.na(permuteScore)] <- 0

  ## Compute the p-value based on bootstrap method
  pValue <- rowSums(abs(permuteScore) >= abs(raw.score)) / 1000
  ## Compute the adjusted p-value. The adjusting method can be reseted
  ## (Refer to p.adjust()).
  pAdjust <- stats::p.adjust(pValue, "fdr")

  score <- .S(raw.score)
  results <- data.frame(
    tcm = colnames(data),
    direction = ifelse(score > 0.3, 'up', ifelse(score > -0.3, 'none', 'down')),
    raw_score = raw.score,
    scaled_score = score,
    Pval = pValue,
    FDR = pAdjust,
    stringsAsFactors = FALSE
  )

  results <- results[order(abs(results$scaled_score), decreasing = TRUE), ]
  results <- results[, c(
    "tcm", "direction", "raw_score",
    "scaled_score", "Pval", "FDR"
  )]

  names(results) <- c("tcm", "direction", "raw_score", "scaled_score", "pvalue", "fdr")

  return(results)

}
################################################################################

#' tcm.SSwithXSum method for Signature Search
#'
#' @inherit tcm.SSwithXCos
#' @export
#'
#' @examples
#' data <- data_logFC
#' upset <- rownames(data_logFC)[1:100]
#' downset <- rownames(data_logFC)[400:550]
#' input <- list(upset = upset, downset = downset)
#' tcm.SSwithXSum(input,data)

tcm.SSwithXSum <- function(input, data, cores = 1L) {

  stopifnot(all(c("upset", "downset") %in% names(input)), is.list(input))
  if (is.data.frame(data)) {data <- as.matrix(data)}
  data <- as.matrix(data)
  upset <- input$upset
  downset <- input$downset

  if (is.null(colnames(data)) || is.null(rownames(data))) {
    stop("Warning: data should have both rownames and colnames!")
  }

  if (!is.null(upset)) {
    num <- sum(upset %in% rownames(data))
    if (num <= 10) {
      stop("the up gene less than 10")
    }
    message(paste(
      sum(upset %in% rownames(data)), "/", length(rownames(data)),
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
      sum(downset %in% rownames(data)), "/", length(rownames(data)),
      "genes in down set share identifiers with reference database"
    ))
    downset <- downset[downset %in% rownames(data)]
  }
  if (is.null(upset) & is.null(downset)) {
    stop("Both upset and downset share zero identifiers with reference database,
          please make sure that at least one share identifiers!")
  }


  if (!is.null(upset)) {
    num <- sum(union(upset,downset) %in% rownames(data))
    if (num <= 10) {
      stop("the matched gene less than 10")
    }
  }

  if (!is.character(upset)) {upset <- as.character(upset)}
  if (!is.character(downset)) {queryUp <- as.character(downset)}
  message(paste0("\n", ">>> Running ", "tcm.SSwithXCos"))
  ## Convert the gene expression matrix to ranked list.
  matrixToRankedList <- function(data) {
    ## Allocate memory for the refList
    refList <- vector("list", ncol(data))
    for(i in 1:ncol(data)) {
      ## Sort the reference gene lists based on fold change of gene expression.
      refList[[i]] <- c(head(data[order(data[, i], decreasing=TRUE), i],
                             n = 500),
                        tail(data[order(data[, i], decreasing=TRUE), i],
                             n = 500))
    }
    return(refList)
  }
  ## The core part for computing the XSum score
  XSum <- function(refList, upset, downset) {
    scoreUp <- sum(refList[match(upset, names(refList))], na.rm = T)
    scoreDown <- sum(refList[match(downset, names(refList))], na.rm = T)
    return(scoreUp - scoreDown)
  }

  ## Prepare the ranked reference lists
  refList <- matrixToRankedList(refMatrix)
  ## Compute the scores for each sample in the reference lists. mcCore is the
  ## number of cores to use for parallel computing. Set it based on your computer.
  raw.score <- mclapply(refList, XSum, upset = upset,
                    downset = downset, mc.cores = cores)
  raw.score <- as.vector(do.call(rbind, raw.score))
  score <- .S(raw.score)

  ## Allocate memory for the permuteScore that are used to compute the p-value.
  ## The permuteNum can be reseted. Notice large permuteNum means low speed.

  library(furrr)
  plan(multisession, workers = 16L) # availableCores()-1
  system.time(
    result <- furrr::future_map(1:1000, function(i) {
      ## Prepare the random query signatures
      bootUp <- sample(rownames(data), size = length(upset))
      bootDown <- sample(rownames(data), size = length(downset))
      ## Compute the random scores for each sample in the reference lists
      bootScore <- mclapply(refList, XSum, upset = bootUp,
                            downset = bootDown, mc.cores = cores)
    }, .progress = TRUE)
  )

  permuteScore <- do.call(cbind,result)

  permuteScore[is.na(permuteScore)] <- 0

  ## Compute the p-value based on bootstrap method
  pValue <- rowSums(abs(permuteScore) >= abs(score)) / permuteNum
  ## Compute the adjusted p-value. The adjusting method can be reseted

  pAdjust <- stats::p.adjust(pValue, "fdr")

  results <- data.frame(
    tcm = colnames(data),
    direction = ifelse(score > 0.3, 'up', ifelse(score > -0.3, 'none', 'down')),
    raw_score = raw.score,
    scaled_score = score,
    Pval = pValue,
    FDR = pAdjust,
    stringsAsFactors = FALSE
  )

  results <- results[order(abs(results$scaled_score), decreasing = TRUE), ]
  results <- results[, c(
    "tcm", "direction", "raw_score",
    "scaled_score", "Pval", "FDR"
  )]

  names(results) <- c("tcm", "direction", "raw_score", "scaled_score", "pvalue", "fdr")

  return(results)
}
################################################################################











