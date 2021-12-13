#' Algorithm for Comorbidity Analysis
#' @description
#' Lists algorithm to be inputed in \code{generateComorbidityAnalysis()} function.
#' @usage
#' \code{listAlgorithms()}
#' @return
#' @export
#' @author Omar Hassoun
#' @examples
#' \code{listAlgorithms()}
listCAlgorithms = function()
{
  return(c("Shared genes", "Uniqueness of shared genes", "Shared ontologies",  "Network-based seperation"))
}
