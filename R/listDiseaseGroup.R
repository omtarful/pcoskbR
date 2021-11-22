#' List Disease Group
#'
#' @param dataset String containing either the miRNAs or Genes dataset. Genes is the default.
#'
#' @return character vector of disease groups in input dataset.
#' @export
#'
#' @examples
listDiseaseGroup = function(dataset = "Genes")
{
  pcoskb = rvest::session("http://pcoskb.bicnirrh.res.in/disease.php")
  pcoskb = rvest::read_html(pcoskb)
  name = pcoskb %>%
    rvest::html_elements("option") %>% # returns multiple nodes that contain the css_selector
    rvest::html_attr("value")
  disease_group = name[2:18]
  if(dataset == "miRNAs")
  {
    #scrape the disease group  from the RNAs
    return(disease_group[c(-5,-7)])
  }
  else if(dataset == "Genes")
  {
   return(disease_group)
  }
  else
  {
    #throw an error
  }

}

