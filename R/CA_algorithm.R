
#' Comorbidity Analysis Tool
#'
#' @description This tool can be used to predict comorbidity based on shared genes and their ontologies,
#' uniqueness of the shared genes, and the network-based separation between them. For more details of
#'  the methodology refer to \href{https://www.nature.com/articles/s41598-020-71418-8}{PCOSKBR2: a database of genes, diseases, pathways, and networks associated with polycystic ovary syndrome}.
#'
#' @param disease_list character vector containing diseases to be analyzed.
#' @param algorithm character vector containing algorithm to be perform the analysis.
#'
#' @return prints a heatmap containing comorbidity between input diseases.
#' @export
#'
#' @examples
generateComorbidityAnalysis = function(disease_list, algorithm)
{
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


