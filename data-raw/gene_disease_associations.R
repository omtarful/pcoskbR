gene_disease_associations = read.csv("inst/extdata/disease_gene_pcoskbr2.csv")
sysdata_filenames <- load("R/sysdata.rda")
save(list = c(sysdata_filenames, "gene_disease_associations"), file = "R/sysdata.rda")
