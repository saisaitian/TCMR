

universe <- rownames(data)
data.matrix <- query2

data2 <- ifelse(data > 0.585, 1, ifelse(data > -0.585, 0, -1))

names(data2)

common.genes <- intersect(rownames(data.matrix), rownames(data))

if (length(common.genes) == 0) {
  stop("None of the query gene identifiers could be found in the reference dataset.",
    call. = FALSE
  )
}

query.common <- abs(matrix(data.matrix[common.genes, ], nrow = length(common.genes)))

sets.common <- abs(matrix(data2[common.genes, ], nrow = length(common.genes)))

coincidence <- Matrix::crossprod(query.common, sets.common)

query.and.sets <- t(matrix(coincidence, ncol = ncol(data2), dimnames = list(names(data.matrix), names(data))))

dimnames(query.and.sets)

query.all <- matrix(rep(Matrix::colSums(abs(data.matrix)), ncol(data2)),
  nrow = ncol(data2), ncol = ncol(data.matrix), byrow = TRUE, dimnames = dimnames(query.and.sets)
)

sets.all <- matrix(rep(Matrix::colSums(abs(data2)), ncol(data.matrix)),
  nrow = ncol(data2), ncol = ncol(data.matrix), byrow = FALSE, dimnames = dimnames(query.and.sets)
)

neither <- length(universe) - sets.all - query.all + query.and.sets

sets.not.query <- length(universe) - query.all - neither

query.not.sets <- query.all - query.and.sets

query.colsum <- sets.not.query + neither

## p-value calculation
p.values <- matrix(
  unlist(
    mclapply(row.names(query.and.sets), function(k) {
      k <- 1
      fisher_p(query.and.sets[k, ], sets.not.query[k, ], query.all[k, ], query.colsum[k, ])
    })
  ),
  ncol = ncol(data.matrix), byrow = TRUE,
  dimnames = list(names(data2), names(data.matrix))
)

lor <- log((query.and.sets * neither) / (query.not.sets * sets.not.query))
lor[query.not.sets == 0] <- Inf
lor[sets.not.query == 0] <- Inf
lor[query.and.sets == 0] <- 0



fisher_p <- function(x, y, m, n, relErr = 1 + 1e-7) {
  ## 'x' and 'y' are entries in the top two cells; 'm' and 'n' are column totals.
  ## Code is excerpted from fisher.test, for efficiency. Note that 'support'
  ## varies in length with the input variables, so vectorization is only possible
  ## via an mapply.
  mapply(
    function(x, y, m, n) {
      k <- x + y
      lo <- max(0, k - n)
      hi <- min(k, m)
      support <- lo:hi
      d <- dhyper(support, m, n, k, log = TRUE)
      d <- exp(d - max(d))
      d <- d / sum(d)
      sum(d[d <= d[x - lo + 1] * relErr])
    },
    x, y, m, n
  )
}
