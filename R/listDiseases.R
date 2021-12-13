#' List Diseases: list diseases according to each dataset and disease group.
#'
#' @param dataset character containing the dataset. Can be either miRNA or Genes.
#' @param disease_group character containing disease group. Use \link[pcoskbR]{listDiseaseGroup()}for more info on diease groups.
#'
#' @return A \code{vector} with the diseases belonging to the disease group and dataset.
#' @export
#'
listDiseases = function(dataset = "Genes", disease_group = "all")
{
  diseases =NULL
  if(checkDiseasesGroup(disease_group, dataset))
    {
    if(dataset == "Genes")
    {
      diseases = browse(page ="Diseases", filter = disease_group)
    }
    else if(dataset == "miRNAs")
    {
      diseases = browse(page = "miRNA Diseases", filter = disease_group)
    }
    else
    {
      #it means they had the wrong dataset
      #returns errors
      stop("Argument 'dataset' is not valid.")
    }
  }
  else
  {
    #it means the diseases group is not  in the dataset
    #throws error
  }
  diseases = diseases[,"Disease"]
  return(as.vector(diseases))
}
#checks if disease group is in dataset
checkDiseasesGroup = function(disease_group, dataset)
{
  if(is.null(disease_group))
  {
    return(TRUE)
  }
  disease_group_list = listDiseaseGroup(dataset)
  if(disease_group %in% disease_group_list)
  {
    return(TRUE)
  }
  else
  {
    return(FALSE)
  }
}

isInDataset = function(disease_list, dataset)
{
  if(length(disease_list) == 1)
    if(disease_list == "PCOS")
      stop("Analysis can't be done solely with PCOS")

  all_diseases = NULL
  if(dataset == "miRNAs")
  {
    all_diseases = gene_disease_associations %>% dplyr::filter(grepl(pattern = "MIR", x = geneSymbol))
    all_diseases = all_diseases$disease_merge
  }
  else if(dataset == "Genes")
  {
    all_diseases = gene_disease_associations %>% dplyr::filter(!grepl(pattern = "MIR", x = geneSymbol))
    all_diseases = all_diseases$disease_merge
  }
  else
  {
   #returns error if it's not in dataset
    stop("Argument 'dataset' must be either \"miRNAs\" or \"Genes\"")
  }
  #returns true if all diseases are in dataset, false otherwise
  return(sum(disease_list %in% all_diseases) == length(disease_list))
}
