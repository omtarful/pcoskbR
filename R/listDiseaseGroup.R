#' List Disease Groups
#' @description Lists disease groups given a dataset.
#' @usage
#' \code{listDiseaseGroup(dataset = c("Genes", "miRNAs"))}
#' @param dataset Dataset. Either "miRNA" or "Genes".
#'
#' @return A \code{vector} containing disease groups.
#' @export
#' @author Omar Hassoun
#' @examples
#' \code{listDiseaseGroup(dataset = "Genes")}
listDiseaseGroup = function(dataset = "Genes")
{
  if(!(dataset != "Genes" | dataset != "miRNAs"))
    stop("Argument 'dataset' must be specified.")
  url = NULL
  if(dataset == "Genes")
    url = "http://pcoskb.bicnirrh.res.in/disease.php"
  if(dataset == "miRNAs")
    url = "http://pcoskb.bicnirrh.res.in/diseasem.php"
  pcoskb = rvest::session(url)
  pcoskb = rvest::read_html(pcoskb)
  name = pcoskb %>%
    rvest::html_elements("#par option") %>% # returns multiple nodes that contain the css_selector
    rvest::html_attr("value")
  return(name)
}

