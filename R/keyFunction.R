
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


