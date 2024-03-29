% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/signature-search.R
\name{tcm.SSwithGCMAP}
\alias{tcm.SSwithGCMAP}
\title{gcmap Method for Signature Search}
\usage{
tcm.SSwithGCMAP(input, data, higher = 1, lower = -1, cores = 1L)
}
\arguments{
\item{input}{A \code{data.frame} containing \code{logFC} value, which \code{rownames} should be gene symbols.}

\item{data}{A \code{data.frame} contains \code{logFC} value data from 103 compounds.}

\item{higher}{A cutoff value. If \code{logFC} larger than or equal to \code{higher} will be included
in the gene set with  \code{+1}.}

\item{lower}{A cutoff value. If \code{logFC} smaller than or equal to \code{higher} will be included
in the gene set to \code{-1}, others should be set to \code{0}.}

\item{cores}{Number of cores to run the task}
}
\value{
A \code{data.frame}.
}
\description{
gcmap Method for Signature Search
}
\examples{
data("data_logFC")
query2 <- data_logFC[1:60, 1, drop = FALSE]
gcmap_kk <- tcm.SSwithGCMAP(input = query2, data = data_logFC)
}
