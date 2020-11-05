
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TCMR

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
<!-- badges: end -->

The goal of TCMR is to …

## Installation

You can install the released version of TCMR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("TCMR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("saisaitian/TCM_Microarray")
```

## Example

This is a basic example which shows you how to solve a common problem.

### Load example dataset

``` r
library(TCMR)
data <- load_example_dataset()
#> Loading dataset GSE85871...
#> Done.
```

### Load analyzed DEG data

``` r
data("AnalyzedDEG")
head(AnalyzedDEG)
#>   id  batch                           vs           filename
#> 1  1 batch1       Glycyrrhizic acid:DMSO GSE85871-DEG-1.rds
#> 2  2 batch1 Hydroxysafflor yellow A:DMSO GSE85871-DEG-2.rds
#> 3  3 batch1         Anhydroicaritin:DMSO GSE85871-DEG-3.rds
#> 4  4 batch1              Hyperoside:DMSO GSE85871-DEG-4.rds
#> 5  5 batch1              Hesperidin:DMSO GSE85871-DEG-5.rds
#> 6  6 batch1                Puerarin:DMSO GSE85871-DEG-6.rds
```

Use a subset of `AnalyzedDEG` to select corresponding DEG results.

``` r
head5_reports <- head(AnalyzedDEG) %>% 
  load_analyzedDEG()

str(head5_reports, max.level = 1)
#> List of 6
#>  $ :Classes 'data.table' and 'data.frame':   12548 obs. of  7 variables:
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ :Classes 'data.table' and 'data.frame':   12548 obs. of  7 variables:
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ :Classes 'data.table' and 'data.frame':   12548 obs. of  7 variables:
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ :Classes 'data.table' and 'data.frame':   12548 obs. of  7 variables:
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ :Classes 'data.table' and 'data.frame':   12548 obs. of  7 variables:
#>   ..- attr(*, ".internal.selfref")=<externalptr> 
#>  $ :Classes 'data.table' and 'data.frame':   12548 obs. of  7 variables:
#>   ..- attr(*, ".internal.selfref")=<externalptr>
```

``` r
one_report <- load_analyzedDEG(2)
head(one_report)
#>    identifier     logFC  AveExpr         t      P.Value adj.P.Val         B
#> 1:     ADAM30 -3.319485 3.497113 -34.74913 4.466436e-05 0.3022539 0.9008748
#> 2:    SYNPO2L  3.293256 1.646628  27.62886 8.994106e-05 0.3022539 0.7949715
#> 3:     IL36RN  3.229622 2.808680  25.72108 1.118756e-04 0.3022539 0.7523525
#> 4:       GYPE -3.338898 2.053063 -23.01748 1.569436e-04 0.3022539 0.6753000
#> 5:       SCTR  2.543554 2.602965  21.24630 2.002828e-04 0.3022539 0.6105549
#> 6:        C1S -2.803254 2.748042 -21.14723 2.031525e-04 0.3022539 0.6065191
```

### Run DEG analysis

``` r
ix <- c(1, 2, 61, 62)
expr <- data$expr[, ix]
group <- data$pdata$perturbagen[ix]

# Run DEG analysis
report <- deg_caller(expr, group = group, level = group[c(3, 1)])
#> Info: Glycyrrhizic acid vs DMSO (reference group)
#> N: DMSO:#2  Glycyrrhizic acid:#2
#> Constructing design matrix...
#> Running DEG analysis with limma...
#> Warning: Zero sample variances detected, have been offset away from zero
#> Reporting results...
#> Done.
head(report)
#>    identifier     logFC  AveExpr         t      P.Value adj.P.Val         B
#> 1:      CDH12  3.411526 1.981608  30.80083 5.441761e-05 0.2925089 1.1921795
#> 2:      NLGN1 -3.204683 1.645785 -26.29158 8.915954e-05 0.2925089 1.0903322
#> 3:       HCRT  3.181742 3.144149  22.99756 1.353037e-04 0.2925089 0.9811453
#> 4:    SOSTDC1  3.582065 1.791032  22.73374 1.402533e-04 0.2925089 0.9706183
#> 5:   SERPINB2 -2.191310 2.182794 -20.95775 1.806605e-04 0.2925089 0.8909550
#> 6:   ATP6V1G2  2.698387 3.463790  20.22094 2.019340e-04 0.2925089 0.8527498
```

See `?deg_caller` for more details.
