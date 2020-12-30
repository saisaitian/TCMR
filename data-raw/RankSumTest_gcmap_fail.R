

data.matrix <- query2

data2 <- ifelse(data>1,1,ifelse(data>-1,0,-1))
## subset objects to shared genes
matched.features <- match(rownames(data.matrix), rownames(data2))
matched.sets <- data2[na.omit(matched.features),]

## extract scores for each gene set
sets.up <- lapply(seq(ncol(matched.sets)),
                  function(x) which(matched.sets[ ,x ] == 1))

sets.down <- lapply(seq(ncol(matched.sets)),
                    function(x) which(matched.sets[ ,x] == -1))


## transform experiment to ranks
experiment <- apply(data.matrix, 2, rank)

## perform wilcox rank sum test
stats <- apply( experiment, 2, function( m ) {
  sapply( seq_along( sets.up ), function( n ) {
    res <- signedRankSumTest(m, sets.up[[n]], sets.down[[n]],
                             adjust.ties=T, input.is.ranks=TRUE)
    list(p.value=res["p"],
         U=res["U"],
         effect=res["z"]
    )
  })
})

## store results
results <- data.frame(
  set = dimnames(data)[2],
  trend = ifelse( unlist(stats[[1]]["effect",]) >= 0, "correlated", "anticorrelated"),
  pval  = unlist(stats[[1]]["p.value",]),
  padj  = p.adjust( unlist( stats[[1]]["p.value",]), method="BH"),
  effect = unlist(stats[[1]]["effect",]),
  nSet = colSums(as.matrix(abs(data))),
  nFound = colSums(as.matrix(abs(matched.sets)))
)

## Apply scaling of scores to full data set
results[,"effect"] <-.connnectivity_scale(results$effect)
results <- results[order(abs(results$effect), decreasing=TRUE), ]

signedRankSumTest <- function(
  statistics,
  index.up,
  index.down=NULL,
  input.is.ranks = FALSE,
  correlation = 0,
  df = Inf,
  adjust.ties=TRUE) {

  if( length( index.up   ) == 0 ) { index.up <- NULL }
  if( length( index.down ) == 0 ) { index.down <- NULL }
  if( is.null(index.up) & is.null(index.down) ) {
    return( c(U=NA, p=1, z=0) )
  } else {
    n <- length(statistics)
    if ( input.is.ranks == FALSE ) {
      r <- rank(statistics)
    } else {
      r <- statistics
    }
    rev.r <- n-r+1
    r1 <- c( r[index.up], rev.r[index.down] )
    n1 <- length(r1)
    n2 <- n - n1
    U <- n1 * n2 + n1 * (n1 + 1)/2 - sum(r1)
    mu <- n1 * n2/2
    if (correlation == 0 || n1 == 1) {
      sigma2 <- n1 * n2 * (n + 1)/12
    }
    else {
      sigma2 <- asin(1) * n1 * n2 + asin(0.5) * n1 * n2 * (n2 - 1)
      + asin(correlation/2) * n1 * (n1 - 1) * n2 * (n2 - 1)
      + asin((correlation + 1)/2) * n1 * (n1 - 1) * n2
      sigma2 <- sigma2/2/pi
    }

    if( adjust.ties == TRUE) {
      TIES <- any( duplicated( r ))
      if (TIES) {
        NTIES <- table(r)
        adjustment <- sum(NTIES * (NTIES + 1) * (NTIES - 1))/(n *
                                                                (n + 1) * (n - 1))
        sigma2 <- sigma2 * (1 - adjustment)
      }
    }

    zlowertail <- (U + 0.5 - mu)/sqrt(sigma2)
    zuppertail <- (U - 0.5 - mu)/sqrt(sigma2)
    p.less <- pt(zuppertail, df = df, lower.tail = FALSE)
    p.greater <- pt(zlowertail, df = df)
    direction <- ifelse( p.less <= p.greater, "less", "greater")

    c(U=U,
      p=ifelse( direction == "less", p.less, p.greater),
      z=ifelse( direction == "less", -zuppertail, -zlowertail)
    )
  }
}
