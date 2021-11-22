
#' Pathway enrichment Analysis:  For selected diseases, this tool gives a list of enriched pathways.
#'
#' @param disease_list character vector containing diseases to be analyzed.
#' @param dataset character containing the dataset. Can be either miRNA or Genes.
#' @param database character vector containing the database from which pathways will be extracted. It can either be "KEGG" or "Reactome".
#'
#' @return
#' @export
#'
#' @examples
generatePathwayTable = function(disease_list, dataset, database = "KEGG")
{
  disease_table = getDiseaseTable(disease_list, dataset)
  gene_population = sapply(disease_table$`Gene Symbol`, strsplit, split = ", ")
  gene_population = unname(unlist(gene_population))
  pathways = NULL
  if(dataset == "miRNAs")
  {
    split_mirna = strsplit(x= gene_population, split = "R")[[1]]
    pattern = paste0("hsa-", split_mirna[1], "R-", split_mirna[2], "")
    filtered_human_data = human_data %>% dplyr::filter(grepl(pattern, x = miRNA, ignore.case = TRUE, perl = TRUE))
    gene_population = unique(filtered_human_data$`Target Gene`)
    #result of querying the miRNA and then downloading
    pathways = getPathways(gene_population, database)
  }
  else if(dataset == "Genes")
  {
    pathways = getPathways(gene_population, database)
  }
  else
  {

  }

  #this is the gene population, it's N in the formula


  pcos_associated_genes = gene_disease_associations %>% dplyr::filter(disease_merge == "PCOS")
  comb = getComb(pathways, gene_population) #pcos genes
  pathways = pathways[,-c(2,4)]
  pathways = cbind(pathways, comb)
  pathways = cbind(pathways, rep(NA, nrow(pathways)))
  colnames(pathways) = c("Pathway name",  "Pathway genes associated with PCOS", "Pathway genes associated with PCOS and the selected diseases", "Hypergeometric probability")
  N = union(countGenes(pathways$`Gene (Gene Symbol)`), gene_disease_associations$geneSymbol) #population of genes in pcos and pathway dataset
  n = gene_population #sample size = number of gnes used in input

  for (i in 1:nrow(pathways))
  {
    K = unique(countGenes(pathways$`Pathway genes associated with PCOS`[i])) #number of pcos genes in ith pathway
    k = intersect(K, gene_population) #number of inpiut genes in sample overlapping iwth genes in ith pathay
    pathways[i,"Hypergeometric probability"] = dhyper(x = length(k), m = length(K), n = length(setdiff(N, K)), k = length(n))
  }
  pathways = pathways[,c(1,3,2,4)]
  return(pathways)
}
#gets patthways for each disease
getPathways = function(gene_population, database)
{
  pathway_table <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(pathway_table) <- c("Pathway", "Pathway ID", "Gene (Gene Symbol)", "Source" )


  pathway_table = pathway_table %>%
    dplyr::mutate(across(everything(), as.character))
  for (gene in gene_population)
  {
    url = paste0("http://pcoskb.bicnirrh.res.in/result1.php?opt=6&qury=", gene)
    gene_pathway_table = getTableSection(url, n_col = 4, tag = "td td td td", n_page = 0)
    colnames(gene_pathway_table) = gene_pathway_table[1,]
    gene_pathway_table = gene_pathway_table[-1,]
    if(nrow(gene_pathway_table) == 0)
    {
      next
    }
    gene_pathway_table = gene_pathway_table %>% dplyr::filter(Source == database)
    #add an aditional grep check
    if(nrow(gene_pathway_table) > 0 )
    {
      gene_pathway_table = gene_pathway_table[checkIfCorrect(gene, gene_pathway_table),]
    }
    else
    {
      next
    }

    pathway_table = dplyr::union(pathway_table, gene_pathway_table)
  }
  return(pathway_table)
}

checkIfCorrect = function(gene, gene_pathway_table)
{
  #check if all the strings in the gene pathway table contain the exact given genE
  boolean_vec = sapply(gene_pathway_table$`Gene (Gene Symbol)`, function(x)
  {
    genes = countGenes(x)
    if(length(which(gene == genes)) == 0)
    {
      return(FALSE)
    }
    else
    {
      return(TRUE)
    }
  })
  return(unname(boolean_vec))
}

getComb = function(pathways, gene_population)
{
  path_pcos = sapply(pathways$`Gene (Gene Symbol)`, function(x){
    genes  = strsplit(x, split = ", ")[[1]]
    genes = unname(sapply(genes, str_replace_all, pattern = "[\r\n\\s]", replacement = ""))
    genes = intersect(gene_population, genes)
    return(paste(genes, collapse = ", "))
  })
  return(unname(path_pcos))
}


countGenes = function(gene_column)
{
  genes = sapply(gene_column, strsplit, split = ", | ")
  genes = unname(unlist(genes))
  genes = stringr::str_replace_all(genes, "[\r\n\\s]" , "")
  genes = unique(genes[which(genes != "")])
  return(genes)
}
#takes the output of the pathway analysis function and returns a network where the nodes are the pathways and
#
viewNetwork = function(pathways)
{
  #create nodes
  nodes = data.frame(id = pathways$`Pathway name`, #this is an id
                     label = pathways$`Pathway name`, #this labels the node with the disease name
                     value = rep(NA, nrow(pathways))) #this gives the size to the node

  #this assigns a size to each node depending on the amount of genes the disease has
  for (row in 1:nrow(nodes))
  {
    genes = countGenes(pathways$`Pathway genes associated with PCOS`[row])
    nodes[row, "value"] = length(genes)
  }
  #create edges
  #this creates all pairs of pathways  and edge dataframe
  edges =  t(combn(pathways$`Pathway name`, 2))
  colnames(edges) = c("from", "to")
  #this manipulates the width of each edge
  value = rep(NA, nrow(edges))
  edges = cbind(edges, value)
  for (row in 1:nrow(edges))
  {
    pathway_A = edges[row,1]
    pathway_B = edges[row,2]
    genes_pathway_A = pathways[which(pathways$`Pathway name` == pathway_A),3]
    genes_pathway_A = countGenes(genes_pathway_A)
    genes_pathway_B = pathways[which(pathways$`Pathway name` == pathway_B),3]
    genes_pathway_B = countGenes(genes_pathway_B)
    shared_genes = intersect(genes_pathway_A, genes_pathway_B)
    edges[row, "value"] = length(shared_genes)
  }
  #feed nodes and edges to the function
  edges <- as.data.frame(edges)
  edges$value= as.numeric(as.character(edges$value))

  valid = which(edges$value == 0)
  if(length(valid) > 0)
  {
    edges = edges[-valid,]
  }

  print(visNetwork::visNetwork(nodes, edges, width = "100%") %>%  visNetwork::visIgraphLayout())
}




