
.onLoad <- function(libname, pkgname) {

  invisible(suppressPackageStartupMessages(library("ComplexHeatmap")))
  invisible(suppressPackageStartupMessages(library("tidyHeatmap")))

  invisible(suppressPackageStartupMessages(
    sapply(c("tibble", "tidyverse", "survival", "survminer", "ggplot2",
             "ggpubr","limma","limSolve","preprocessCore","e1071","GSVA"),
           requireNamespace, quietly = TRUE)
  ))
}



##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0("==========================================================================\n",
                " ", pkgname, " v", pkgVersion, "  ",

                "  For help: https://github.com/XXXX", "\n\n")

  citation <-paste0(" If you use ", pkgname, " in published research, please cite:\n",
                    " SSSSS*.\n",
                    " TCMR: XXXXXXXXXXXXX \n ",
                    " XXXXXXXXXXXXXX. XXXXXXXXXXXX. XXXXXXX. \n",
                    # " XXXX, 2020", "\n",
                    " DOI: XXXXXX\n",
                    # " PMID:  ","\n",
                    "==========================================================================")

  packageStartupMessage(paste0(msg, citation))
}

