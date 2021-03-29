
#' @title Imports the table from the PCOSKB website
#' @description Imports the table from the PCOSKB website -in the form of a dataframe- into the R workspace.
#' @return the table in the PCOSKB website but in the form of a dataframe
#' @export
#' @import dplyr rvest
#' @examples
getPCOS = function(){
  url = "./data/Genes associated with PCOS.html"

  pcoskb = read_html(url)

  name = pcoskb %>%
    html_nodes("form table tr") %>%
    html_text()
  column_names = as.vector(strsplit(name[1], split = "\\n\\s+"))[[1]][2:7]
  mati <- matrix(nrow = 1, ncol = 6)
  mati = rbind(mati, column_names)

  for (n in 2:length(name)) {
    row = as.vector(strsplit(name[n], split = "\\n\\s+"))[[1]][2:7]
    mati = rbind(mati, row)
  }

  mati =mati[3:nrow(mati),]

  PRG.df = as.data.frame(mati)
  names(PRG.df) = column_names
  return(PRG.df)
}

