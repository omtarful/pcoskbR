## code to prepare `DATASET` dataset goes here
string_db <- STRINGdb::STRINGdb$new( version="11", species=9606, input_directory="inst/extdata")

usethis::use_data(string_db,internal = TRUE, overwrite = TRUE)
