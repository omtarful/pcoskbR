source("./R/browse.R")

#' lisFilters: Prints the filters for the page inputted.
#'
#' @param page It's a string for the page you want to obtain filters to.
#'
#' @return Prints all the filters for the page inputted.
#' @export
#'
#' @examples
listFilters = function(page)
{
  filters = NULL
  if(page == "Genes" | page == "miRNA")
  {
    filters = getText(url = "http://pcoskb.bicnirrh.res.in/gene.php", css_selector = "#limitt")
  }
  else if(page == "Diseases" | page == "miRNA Diseases")
  {
    filters = getText(url = "http://pcoskb.bicnirrh.res.in/disease.php", css_selector = "#par")
  }
  else if(page == "Pathways")
  {
    filters = getText(url = "http://pcoskb.bicnirrh.res.in/pathway_pcos.php", css_selector = "#limitt")
  }
  else if(page == "Ontologies")
  {
    filters = getText(url = "http://pcoskb.bicnirrh.res.in/go_d.php", css_selector = "#limitt")
  }
  else if(page == "SNPs")
  {
    print("This page has no filters")
  }
  else
  {
    return("Page doesn't exist")
  }

  if(!is.null(filters))
  {
    result = str_extract_all(filters, pattern = "[A-Z][a-z ]+")[[1]]
    if(page == "Pathways")
    {
      result = append(result, "KEGG", after = 1)
    }
  }
  else
  {
    result = NULL
  }
  cat(result, sep = "\n")
}
