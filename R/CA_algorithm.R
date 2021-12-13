
#' Comorbidity Analysis Tool
#'
#' @description Given a list of diseases used as input, this tool can be used to predict comorbidity according to an algorithm. For more details of
#'  the methodology refer to \href{https://www.nature.com/articles/s41598-020-71418-8}{PCOSKBR2: a database of genes, diseases, pathways, and networks associated with polycystic ovary syndrome}.
#' @usage
#' \code{generateComorbidityAnalysis(disease_list =  c("Psychomotor Disorders", "Psychosexual Disorders" , "Pubertal Disorder"),
#'  algorithm = c("Shared genes", "Uniqueness of shared genes", "Shared ontologies", "Network-based seperation"))}
#' @param disease_list Diseases for which comorbidity with PCOS will be analyzed. A possible list of diseases can be retrieved using the \code{listDiseases} function.
#' @param algorithm One of four algorithms used to perform the analysis. A possible list of algorithms can be retrieved using the \code{listAlgorithms} function.
#'
#' @author Omar Hassoun
#' @return A \code{matrix} with scores in each cell indicating the risk of comorbidity between two diseases and generates a heatmap with the risks.
#' @export
#'
#' @examples
#' \code{generateComorbidtiyAnalysis(disease_list = c("Psychomotor Disorders", "Psychosexual Disorders" , "Pubertal Disorder"),
#'                                    algorithm = "Shared genes")}
generateComorbidityAnalysis = function(disease_list, algorithm)
{
  dataset = "Genes"
  if(missing(disease_list))
    stop("Argument 'disease_list' must be specified")
  if(!isInDataset(dataset = dataset, disease_list = disease_list))
    stop("Argument 'disease_list' contains one or more diseases not found in dataset")
  if(!algorithm %in% listCAlgorithms())
    stop("Argument 'algorithm' is not valid.")

  disease_list = append(disease_list, "PCOS", after = 0)
  score_matrix = matrix(nrow = (length(disease_list)-1), ncol = (length(disease_list)-1))
  rownames(score_matrix) = disease_list[-length(disease_list)]
  colnames(score_matrix) = disease_list[-1]
  if(algorithm == "Shared genes")
  {
    for (row_name in disease_list[-length(disease_list)])
    {
      for (col_name in disease_list[-1])
      {
        score_matrix[row_name,col_name] = CA_sharedGenes(row_name, col_name)
      }
    }

  }
  else if(algorithm == "Uniqueness of shared genes")
  {
    for (row_name in disease_list[-length(disease_list)])
    {
      for (col_name in disease_list[-1])
      {
        score_matrix[row_name,col_name] = getUniquenessCScore(row_name, col_name)
      }
    }
  }
  else if(algorithm == "Shared ontologies")
  {
    for (row_name in disease_list[-length(disease_list)])
    {
      for (col_name in disease_list[-1])
      {
        score_matrix[row_name,col_name] = getCA_GOScore(row_name, col_name)
      }
    }
  }
  else if(algorithm == "Network-based seperation")
  {
    for (row_name in disease_list[-length(disease_list)])
    {
      for (col_name in disease_list[-1])
      {
        score_matrix[row_name,col_name] = getCA_interactomeScore(row_name,col_name)
      }
    }

  }
  else {
    #throw error
  }

  print(generateHeatmap(score_matrix))
  return(score_matrix)
}


generateHeatmap = function(score_matrix)
{
  score_matrix_mod = reshape::melt(score_matrix)

  ggp = ggplot2::ggplot(data = score_matrix_mod, ggplot2::aes(X1, X2)) +
    ggplot2::geom_tile(ggplot2::aes(fill = value)) +
    ggplot2::scale_fill_gradient(low = "yellow", high = "dark orange", na.value = "white") +
    ggplot2::coord_flip() +
    ggplot2::geom_text(ggplot2::aes(label = round(value, 1))) +
    ggplot2::theme(axis.title = ggplot2::element_blank())
  return(ggp)
}


