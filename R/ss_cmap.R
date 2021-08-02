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
#' cmap_kk <- ss_cmap(input = input, data = data_logFC)
ss_cmap <- function(input, data, cores = 1L) {
  stopifnot(all(c("upset", "downset") %in% names(input)), is.list(input))
  data <- as.matrix(data)
  upset <- input$upset
  downset <- input$downset
  if (!is.null(upset)) {
    num <- sum(upset %in% rownames(data))
    if (num <= 10) {
      stop("the up gene less than 10")
    }
    message(paste0("\n", ">>> Running ", "ss_cmap"))
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
