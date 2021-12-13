human_data = readxl::read_xlsx("inst/extdata/NCBI2Reactome2.txt")
sysdata_filenames <- load("R/sysdata.rda")
save(list = c(sysdata_filenames, "ncbi2reactome"), file = "R/sysdata.rda")

