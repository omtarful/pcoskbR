human_data = readxl::read_xlsx("inst/extdata/hsa_MTI.xlsx")
sysdata_filenames <- load("R/sysdata.rda")
save(list = c(sysdata_filenames, "human_data"), file = "R/sysdata.rda")
