library("hpar")
data("hpaNormalTissue")
sysdata_filenames <- load("R/sysdata.rda")
save(list = c(sysdata_filenames, "hpaNormalTissue"), file = "R/sysdata.rda")
