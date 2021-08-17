

#' listPages: Prints the pages in PCOSKBR2
#'
#' @return  Prints pages you can search for in PCCOSKBR2.
#' @export
#'
#' @examples
listPages = function()
{
  cat(c("Genes", "miRNA", "Diseases", "miRNA Diseases", "SNPs", "Pathways", "Ontologies"), sep = "\n")
}
