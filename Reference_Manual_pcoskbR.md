<!-- toc -->

December 13, 2021

# DESCRIPTION

```
Package: pcoskbR
Type: Package
Title: What the Package Does (Title Case)
Version: 0.1.0
Author: Omar Hassoun
Maintainer: Omar Hassoun <omtarful@gmail.com>
Description: More about what it does (maybe more than one line)
    Use four spaces when indenting paragraphs within the Description.
License: GPL-2
Encoding: UTF-8
LazyData: true
imports:
  rvest (>= 1.0.0),
  stringr (>= 1.4.0),
  ggplot2 (>= 3.3.5),
  reshape (>= 0.8.8)
  dplyr (>= 1.0.5),
  STRINGdb (>= 2.2.2),
  igraph (>= 1.2.6),
  hpar (>= 1.32.1),
  biomaRt (>= 2.46.3),
  visNetwork (>= 2.1.0),
  plotly (>= 4.10.0),
  readxl (>= 1.3.1)
Suggets:
  knitr,
  rmarkdown
RoxygenNote: 7.1.2
Depends: 
    R (>= 2.10)
Suggests: 
    rmarkdown,
    knitr
VignetteBuilder: knitr```


# `browse`

Browse: Retrieves information from the PCOSKB database


## Description

This function scrapes data from PCOSKB. Given a set of filters and corresponding values, it retrieves information about genes, miRNAs, SNPs, diseases, pathways, and ontology associated with PCOS.


## Usage

```r
`browse(page = c("Genes", "miRNA", "Diseases", "miRNA Diseases", "SNPs", "Pathways", "Ontologies"),`
```


## Arguments

Argument      |Description
------------- |----------------
`page`     |     Page from which data is retrieved. A possible list of pages can be retreived using the [listPage](#listpage) function.
`filter`     |     Filter (one) that should be used in the query. A possible list of pages can be retieved using the [listFilters](#listfilters) function.
`search_qry`     |     Searches a keyword in the query.


## Value

A `data.frame` . Size and attributes depend on query.


## Author

Omar Hassoun


## Examples

```r
`browse(page = "Genes", filter = "Manually curated", search_qry = "Fem")`
```


