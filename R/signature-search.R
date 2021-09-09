#' CMAP method for Signature Search
#'
#' @param input A list containing UP (in `upset` element) and DOWN (in `downset` element) genes.
#' @param data A `data.frame` containing `logFC` value data from 103 compounds.
#' @param cores Number of cores to run the task
#'
#' @return A data.frame.
#' @export
#' @importFrom ExperimentHub ExperimentHub
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
    message(paste0("\n", ">>> Running ", "tcm.SSwithCMAP"))
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
  input$upset <- upset
  input$downset <- downset

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
  eh <- suppressMessages(ExperimentHub::ExperimentHub())
  CSnull <- suppressMessages(eh[["EH3234"]])
  CSnull[CSnull[, "Freq"] == 0, "Freq"] <- 1
  myrounding <- max(nchar(as.character(CSnull[, "WTCS"]))) - 3
  camp_round <- round(as.numeric(raw.score), myrounding)

  CS_pval <- vapply(camp_round, function(x) {
    sum(CSnull[abs(CSnull[, "WTCS"]) > abs(x), "Freq"]) / sum(CSnull[, "Freq"])
  }, FUN.VALUE = numeric(1))
  CS_fdr <- stats::p.adjust(CS_pval, "fdr")
  result <- data.frame(
    set = colnames(data),
    trend = ifelse(score >= 0, "up", "down"),
    raw_score = raw.score,
    scaled_score = score,
    Pval = CS_pval,
    FDR = CS_fdr,
    N_upset = length(upset),
    N_downset = length(downset), stringsAsFactors = FALSE
  )
  result <- result[order(abs(result$scaled_score), decreasing = TRUE), ]
  names(result) <- c("tcm", "direction", "raw_score", "scaled_score", "pvalue", "fdr", "Nupset", "Ndownset")
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
  message(paste0("\n", ">>> Running ", "tcm.SSwithCorrelation"))
  message(paste(
    num, "/", length(rownames(data)),
    "genes in input share identifiers with reference database"
  ))
  if (is.null(input)) {
    stop(" Input is NULL !")
  }

  res_list <- NULL
  common <- intersect(rownames(input), rownames(data))
  input2 <- input[common, , drop = F]
  data2 <- data[common, ]

  for (i in 1:ncol(data2)) {
    tt <- stats::cor.test(input2[, 1], data2[, i], method = method)
    p <- tt$p.value
    cor <- tt$estimate
    direction <- cor
    direction[cor >= 0] <- "up"
    direction[cor < 0] <- "down"
    res <- data.frame(
      tcm = colnames(data2)[i],
      direction = direction,
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

  result <- result[order(abs(result$scaled_score), decreasing = TRUE), ]
  result <- result[, c(
    "tcm", "direction", "raw_score",
    "scaled_score", "p", "fdr", "Nset"
  )]
  names(result) <- c("tcm", "direction", "raw_score", "scaled_score", "pvalue", "fdr", "Nset")
  rownames(result) <- NULL
  return(result)
}

#' CoreGx method for Signature Search
#' @inherit tcm.SSwithCorrelation
#' @inheritParams tcm.SSwithCMAP
#' @export
#'
#' @examples
#' data("data_logFC")
#' query2 <- data_logFC[1:60, 1, drop = FALSE]
#' txp <- tcm.SSwithCoreGx(query2, data_logFC[1:10])
#' @testexamples
#' expect_is(txp, "data.frame")

tcm.SSwithCoreGx <- function(input, data, cores = 1L) {
  stopifnot(is.data.frame(input))
  num <- sum(rownames(input) %in% rownames(data))
  if (num <= 10) {
    stop("the commom gene less than 10")
  }
  message(paste0("\n", ">>> Running ", "tcm.SSwithCoreGx"))
  message(paste(
    num, "/", length(rownames(data)),
    "genes in input share identifiers with reference database"
  ))

  res_list <- NULL
  common <- intersect(rownames(input), rownames(data))
  input2 <- input[common, , drop = F]

  res_list <- parallel::mclapply(
    seq_len(ncol(data)),
    function(i) {
      tt <- suppressWarnings(CoreGx::connectivityScore(data[, i, drop = F], input2, method = "fgsea", nperm = 1e4, nthread = 1))
      # tt <- suppressWarnings(PharmacoGx:::connectivityScore(data[, i, drop = F], input2, method = "fgsea", nperm = 1e4, nthread = 1))
      direction <- tt[1]
      direction[direction >= 0] <- "up"
      direction[direction < 0] <- "down"
      res <- data.frame(
        tcm = colnames(data)[i],
        direction = direction,
        raw_score = tt[1],
        p = tt[2],
        Nset = length(common),
        stringsAsFactors = FALSE
      )
    },
    mc.cores = cores
  )

  result <- do.call(rbind, res_list)
  result$fdr <- stats::p.adjust(result$p, "fdr")
  result$scaled_score <- .S(result$raw_score)
  result <- result[order(abs(result$scaled_score), decreasing = TRUE), ]
  result <- result[, c(
    "tcm", "direction", "raw_score",
    "scaled_score", "p", "fdr", "Nset"
  )]
  names(result) <- c("tcm", "direction", "raw_score", "scaled_score", "pvalue", "fdr", "Nset")
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
#' @importFrom ExperimentHub ExperimentHub
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

  message(paste0("\n", ">>> Running ", "tcm.SSwithGCMAP"))
  message(paste(
    num, "/", length(rownames(data)),
    "genes in input share identifiers with reference database"
  ))
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

  eh <- suppressMessages(ExperimentHub::ExperimentHub())
  CSnull <- suppressMessages(eh[["EH3234"]])
  CSnull[CSnull[, "Freq"] == 0, "Freq"] <- 1
  myrounding <- max(nchar(as.character(CSnull[, "WTCS"]))) - 3
  gcmap_round <- round(as.numeric(raw.score), myrounding)

  CS_pval <- vapply(gcmap_round, function(x) {
    sum(CSnull[abs(CSnull[, "WTCS"]) > abs(x), "Freq"]) / sum(CSnull[, "Freq"])
  }, FUN.VALUE = numeric(1))
  CS_fdr <- stats::p.adjust(CS_pval, "fdr")

  score <- matrix(raw.score, ncol = ncol(input))
  ## store results
  results <- data.frame(
    tcm = colnames(data2),
    direction = ifelse(score[, 1] >= 0, "up", "down"),
    raw_score = score[, 1],
    Pval = CS_pval,
    FDR = CS_fdr,
    Nset = colSums(as.matrix(abs(matched.sets)))
  )


  ## Apply scaling of scores to full data set
  results[, "scaled_score"] <- .S(results$raw_score)

  results <- results[order(abs(results$scaled_score), decreasing = TRUE), ]
  results <- results[, c(
    "tcm", "direction", "raw_score",
    "scaled_score", "Pval", "FDR", "Nset"
  )]

  results <- results[order(abs(results$scaled_score), decreasing = TRUE), ]
  names(results) <- c("tcm", "direction", "raw_score", "scaled_score", "pvalue", "fdr", "Nset")

  return(results)
}

#' Lincs method for Signature Search
#'
#' @inherit tcm.SSwithCMAP
#'
#' @export
#' @importFrom ExperimentHub ExperimentHub
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
  if (!is.null(upset)) {
    num <- sum(upset %in% rownames(data))
    if (num <= 10) {
      stop("the up gene less than 10")
    }
    message(paste0("\n", ">>> Running ", "tcm.SSwithLINCS"))
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
  names(result) <- c("tcm", "direction", "raw_score", "scaled_score", "pvalue", "fdr", "Nupset", "Ndownset")
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
