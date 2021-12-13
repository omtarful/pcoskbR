#' listPages
#' @description Lists pages available in PCOSKBR2.
#' @usage \code{listPages()}
#' @return A \code{vector} of pages.
#' @author Omar Hassoun
#' @export
#'
#' @examples
#' \code{listPages()}
listPages = function()
{
  pages = c("Genes", "miRNA", "Diseases", "miRNA Diseases", "SNPs", "Pathways", "Ontologies")
  return(pages)
}
