getDiseaseTable =  function(disease_list, dataset = "Genes") {
  #add function to check if every disease is part of the dataset
  disease_table = NULL #initialize table of diseases u need to do the analysis
  page = NULL
  if (dataset == "miRNAs") {
    page =  "miRNA Diseases"
  } else if (dataset == "Genes") {
    page = "Diseases"
  } else {
    #return error
  }

  for (disease in disease_list) {
    row = browse(page, search_qry = disease )
    if (nrow(row) > 1) {
      row = row[which(row$Disease  == disease),]
    }
    disease_table = rbind(disease_table, row)
  }
  return(disease_table)
}

getNetworkAnalysis = function(disease_list, dataset)
{
  disease_table = getDiseaseTable(disease_list, dataset)
  ##this creates the node table
  nodes = data.frame(id = disease_table$Disease, #this is an id
                     label = disease_table$Disease, #this labels the node with the disease name
                     value = rep(NA, nrow(disease_table))) #this gives the size to the node
  #this assigns a size to each node depending on the amount of genes the disease has
  for (row in 1:nrow(nodes))
  {
    genes = stringr::str_split(disease_table[row, 2], pattern = ", ")[[1]]
    nodes[row, "value"] = length(genes)
  }
  #this creates all pairs of diseases  and edge dataframe
  edges =  t(combn(disease_table$Disease, 2))
  colnames(edges) = c("from", "to")
  #this manipulates the width of each edge
  value = rep(NA, nrow(edges))
  edges = cbind(edges, value)
  for (row in 1:nrow(edges))
  {
    disease_A = edges[row,1]
    disease_B = edges[row,2]
    genes_disease_A = disease_table[which(disease_table$Disease == disease_A),2]
    genes_disease_A = strsplit(genes_disease_A, split = ", ")[[1]]
    genes_disease_B = disease_table[which(disease_table$Disease == disease_B),2]
    genes_disease_B = strsplit(genes_disease_B, split = ", ")[[1]]
    shared_genes = intersect(genes_disease_A, genes_disease_B)
    edges[row, "value"] = length(shared_genes)
  }

  #figure out how to work on edges data.frame
  edges <- as.data.frame(edges)
  edges$value= as.numeric(as.character(edges$value))

  valid = which(edges$value == 0)
  if(length(valid) > 0)
  {
    edges = edges[-valid,]
  }
  print(visNetwork::visNetwork(nodes, edges, width = "100%"))
  return(disease_table)
}


