#' List Diseases: list diseases according to each dataset and disease group.
#'
#' @param dataset character containing the dataset. Can be either miRNA or Genes.
#' @param disease_group character containing disease group. Use \link[pcoskbR]{listDiseaseGroup()}for more info on diease groups.
#'
#' @return character vector with the diseases belonging to the disease group and dataset.
#' @export
#'
#' @examples
listDiseases = function(dataset, disease_group)
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
