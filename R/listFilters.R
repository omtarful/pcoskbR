#' lisFilters
#' @description Lists or search filters available for a given page in PCOSKBR.
#' @usage
#' \code{listFilters(page)}
#' @param page Page. A possible list of algorithms can be retrieved using the \code{listPages} function.
#'
#' @return A \code{vector} containing functions.
#' @export
#' @author Omar Hassoun
#' @examples
#' \code{listFilters(page = "Genes")}
listFilters = function(page)
{
  filters = NULL
  if(page == "Genes" | page == "miRNA")
  {
    url = "http://pcoskb.bicnirrh.res.in/gene.php"
    tag = "#limitt option"
  }
  else if(page == "Diseases" | page == "miRNA Diseases")
  {
    url = "http://pcoskb.bicnirrh.res.in/disease.php"
    tag = "#par option"
  }
  else if(page == "Pathways")
  {
    url = "http://pcoskb.bicnirrh.res.in/pathway_pcos.php"
    tag = "#limitt option"
  }
  else if(page == "Ontologies")
  {
    url = "http://pcoskb.bicnirrh.res.in/go_d.php"
    tag = "#limitt option"
  }
  else if(page == "SNPs")
  {
    print("This page has no filters")
  }
  else
  {
    stop("Argument 'Page' is not valid.")
  }
  #result
  pcoskb = rvest::session(url)
  pcoskb = rvest::read_html(pcoskb)
  result = pcoskb %>%
    rvest::html_elements(tag) %>%
    rvest::html_attrs()
  result = unname(unlist(result))
  if(!is.null(result))
  {
    result[which(result == "")] = "all"
  }
  else
  {
    result = NULL
  }
  return(result)
}
