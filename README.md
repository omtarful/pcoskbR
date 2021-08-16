
<!-- README.md is generated from README.Rmd. Please edit that file -->
# pcoskbR

pcoskbR provides an R interface to acess datasets related to PCOS in PCOSKBR2. As well as, data mining tools for comorbidity prediction and identification of enriched pathways and hub genes.

## Installation

You can install the released version of pcoskbR from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("pcoskbR")
```

## Searching for information about PCOS

Here we select the page and the filter to get a dataframe of the result from PCOSKBR2.

``` r
library(pcoskbR)
#> Loading required package: rvest
#> Loading required package: stringr
manually_curated_genes = browse(page = "Genes", filter = "Manually curated")
head(manually_curated_genes)
#>       Gene Symbol Entrez ID                                           Aliases
#> row         ABCA1        19 ABC-1, ABC1, CERP, HDLCQTL13, HDLDT1, HPALP1, TGD
#> row.1       ACACA        31                     ACAC, ACACAD, ACC, ACC1, ACCA
#> row.2         ACE      1636                            ACE1, CD143, DCP, DCP1
#> row.3       ACTA2        59                                             ACTSA
#> row.4        ACTB        60                                  BRWS1, PS1TP5BP1
#> row.5       ACTG1        71                ACT, ACTG, DFNA20, DFNA26, HEL-176
#>                                       Gene Name Chromosomal Location
#> row   ATP binding cassette subfamily A member 1               9q31.1
#> row.1              Acetyl-CoA carboxylase alpha                17q12
#> row.2           Angiotensin I converting enzyme              17q23.3
#> row.3              Actin alpha 2, smooth muscle             10q23.31
#> row.4                                Actin beta               7p22.1
#> row.5                             Actin gamma 1              17q25.3
#>            Record type
#> row   Manually curated
#> row.1 Manually curated
#> row.2 Manually curated
#> row.3 Manually curated
#> row.4 Manually curated
#> row.5 Manually curated
```
