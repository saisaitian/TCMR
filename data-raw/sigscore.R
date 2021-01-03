
sigscore <- function(a, b,
                     method = c("gsea"),
                     nperm = 10000,
                     nthread = 1,
                     seed = 10, ...) {
  if (!is(a, "matrix")) {
    a <- as.matrix(a)
  }
  if (!is(b, "matrix")) {
    b <- as.matrix(b)
  }

  if (nrow(b) >= nrow(a)) {
    stop("Signature (b) should not be larger than signature (a)")
  }

  if (is.null(rownames(a)) || is.null(rownames(b)) ||
    !length(intersect(rownames(b), rownames(a)))) {
    stop("Row names of a and b are either missing or have no intersection")
  }

  if (nperm < 100) {
    stop("The minimum number of permutations for permutation testing is 100")
  }
  b <- b[!is.na(b[, 1]), , drop = FALSE]
  a <- a[!is.na(a[, 1]), , drop = FALSE]
  gset <- cbind(gene = rownames(b), set = ifelse(as.numeric(b[, 1]) >= 0, "UP", "DOWN"))
  gset <- piano::loadGSC(gset)
  set.seed(seed)
  nes <- piano::runGSA(
    geneLevelStats = a[, 1], geneSetStat = "gsea",
    gsc = gset, nPerm = nperm + (nperm %% nthread), ncpus = nthread,
    verbose = FALSE, adjMethod = "none"
  )

  nes$Dir <- nes$pDistinctDirUp
  nes$Dir[is.na(nes$pDistinctDirUp), 1] <- nes$pDistinctDirDn[is.na(nes$pDistinctDirUp), 1]
  nes.up <- c(
    nes$statDistinctDir[which(names(nes$gsc) == "UP"), 1],
    nes$Dir[which(names(nes$gsc) == "UP"), 1]
  )
  nes.down <- c(
    nes$statDistinctDir[which(names(nes$gsc) == "DOWN"), 1],
    nes$Dir[which(names(nes$gsc) == "DOWN"), 1]
  )
  if (length(nes.up) == 0) {
    score <- c(sigscore = -nes.down[1], p = nes.down[2])
  } else if (length(nes.down) == 0) {
    score <- c(sigscore = nes.up[1], p = nes.up[2])
  } else if (complete.cases(cbind(nes.up, nes.down)) &&
    sign(nes.up[1]) != sign(nes.down[1])) {
    score <- c(sigscore = (nes.up[1] - nes.down[1]) / 2, p = survcomp::combine.test(p = c(
      nes.up[2],
      nes.down[2]
    ), method = "fisher", na.rm = TRUE))
  } else {
    score <- c(score = 0, p = 1)
  }
  return(score)
}




sigscore(data[2], data[1:10, 3, drop = F])

a <- data[2]
b <- data[1:10, 3, drop = F]
