

#' listPages: Prints the pages in PCOSKBR2
#'
#' @return  returns vector with pages you can search for in PCCOSKB.
#' @export
#'
#' @examples
listPages = function()
{
  pages = c("Genes", "miRNA", "Diseases", "miRNA Diseases", "SNPs", "Pathways", "Ontologies")
  return(pages)
}
