#' Gene Network Analysis:  It provides an output of a list of genes with their network properties and tissue specific interactions for each gene.  tab.
#'
#' @param disease_list character vector containing diseases to be analyzed.
#' @param database character vector containing the database from which pathways will be extracted. It can either be "KEGG" or "Reactome".
#'
#' @return a dataframe containing a list of genes with network properties and tissue specific interactions.
#' @export
#'
#' @examples
getGeneNetworkAnalysis = function(disease_list, database)
{
  dataset = "Genes"
  disease_table = getDiseaseTable(disease_list, dataset)
  gene_population = sapply(disease_table$`Gene Symbol`, strsplit, split = ", ")
  gene_population = unname(unlist(gene_population))
  #this loop gets the pathway for each gene in gene_population
  pathways = getPathways(gene_population, database)
  gene_population  = countGenes(pathways$`Gene (Gene Symbol)`)
  #this is what gets the values for
  gene_network_table = generateGeneNetworkTable(gene_population)
  #use the new gene population
  pathway_table = advancedSearch(gene_network_table$gene_population) #pathway table
  pathway_table = pathway_table %>% dplyr::filter(Source == database)
  #adds pathways to genes
  gene_network_table$STRING_id = ""
  for (row in 1:length(pathway_table$`Gene (Gene Symbol)`))
  {
    #thee are the gene contained in the row
    genes_in_row = countGenes(pathway_table$`Gene (Gene Symbol)`[row])
    genes = intersect(genes_in_row, gene_network_table$gene_population)
    positions = match(genes, gene_network_table$gene_population)
    for (i in positions)
    {
      gene_network_table[i,2] = paste0(gene_network_table[i,2], pathway_table$Pathway[row], ", ")
    }
  }
  #generate the network for the function
  colnames(gene_network_table) = c("Gene symbol", "KEGG pathways", "Degree", "Closeness centrality", "Betweeness centrality")
  gene_network_table = cbind(gene_network_table, isHubGene = rep(NA, nrow(gene_network_table)))
  gene_network_table = cbind(gene_network_table, isBottleNeckGene = rep(NA, nrow(gene_network_table)))
  #this determines which genes are hub and which are bottlneck
  for (row in 1:nrow(gene_network_table)) {
    #: Degree>(Mean of Degree+(2* Standard Deviation)) OR Closeness centrality>(Mean of closeness centrality+(2* Standard Deviation))
    gene_network_table[row,6] = if ((gene_network_table[row,3] > (mean(gene_network_table$Degree)+(2*sd(gene_network_table$Degree)) )) | (gene_network_table[row,4] > (mean(gene_network_table$`Closeness centrality`)+(2*sd(gene_network_table$`Closeness centrality`)) )) ) TRUE else FALSE
    #Degree>(Mean of Degree+(2* Standard Deviation)) OR Closeness centrality>(Mean of closeness
    gene_network_table[row,7] = if ((gene_network_table[row,3] < mean(gene_network_table$Degree)) & (gene_network_table[row,5] > mean(gene_network_table$`Betweeness centrality`))) TRUE else FALSE
  }
  return(gene_network_table)
}


generateGeneNetworkTable = function(gene_population)
{
  #maps the genes to stringDB ids
  gene_network_table <- string_db$map(as.data.frame(gene_population), "gene_population", removeUnmappedRows = TRUE )

  subgraph = string_db$get_subnetwork(gene_network_table$STRING_id) #generates an igraph object of the gene population  #sorts the gene network table so u can bind it to the order graph charasteristics
  subgraph = removeSmallComponent(subgraph)
  gene_network_table = gene_network_table[order(gene_network_table$STRING_id),]
  #this is degree
  obj_degree = igraph::degree(subgraph)
  obj_degree = obj_degree[order(names(obj_degree))]
  gene_network_table = cbind(gene_network_table, obj_degree) ##this is the one you will manipulate
  #this is closeness centrality
  obj_closeness = igraph::closeness(subgraph)
  obj_closeness = obj_closeness[order(names(obj_closeness))]
  gene_network_table = cbind(gene_network_table, obj_closeness)
  #this is betweeness centrality
  obj_betweeness = igraph::betweenness(subgraph, directed = FALSE, normalized = TRUE )
  obj_betweeness = obj_betweeness[order(names(obj_betweeness))]
  gene_network_table = cbind(gene_network_table, obj_betweeness)
  return(gene_network_table)

}

#' View Interactions: Shows a network of interactions for an input gene.
#'
#' @param gene A character containing the gene name (symbol)
#' @param tissue character containing tissue type.
#'
#' @return
#' @export
#'
#' @examples
viewInteractions = function(gene, tissue = "all")
{
  #maps gene as string id
  single_gene <- string_db$map(as.data.frame(gene), "gene", removeUnmappedRows = TRUE )
  neighbors = string_db$get_neighbors(single_gene$STRING_id)
  #neighbors of gene
  #check interaction between gene and neuighbors
  gene_and_neighbors = append(neighbors, single_gene$STRING_id)
  #check interaction between genes and neighbors
  gene_interactions = string_db$get_interactions(gene_and_neighbors)

  filtered_interactions = gene_interactions %>% dplyr::filter(from == single_gene$STRING_id | to == single_gene$STRING_id)
  all_genes = union(filtered_interactions$to, filtered_interactions$from) ##all genes INS interacts with

  relations <- data.frame(from=rep(gene, length(all_genes)), to=all_genes)
  #converts to ensembl ids usable in biomart
  all_genes_ids = sapply(all_genes, function(x)
  {
    return(substring(x, 6))
  })
  all_genes_ids = unname(all_genes_ids)
  #uses ensembl to comnvert ensembl ids to Gene names
  ensembl <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  G_list <- biomaRt::getBM(filters = "ensembl_peptide_id",
                           attributes = c("ensembl_peptide_id", "hgnc_symbol"),
                           values = all_genes_ids, mart = ensembl)
  tissues = c("adipose tissue", "adrenal gland", "appendix",
              "bone marrow", "breast", "bronchus", "caudate", "cerebellum", "cerebral cortex",  "cervix, uterine",
              "colon", "duodenum", "endometrium 1", "endometrium 2",  "epididymis",
              "esophagus", "fallopian tube","gallbladder", "heart muscle", "hippocampus",
              "kidney", "liver", "lung", "lymph node", "nasopharynx","oral mucosa",
              "ovary", "pancreas", "parathyroid gland", "placenta", "prostate", "rectum",
              "salivary gland", "seminal vesicle", "skeletal muscle", "skin 1", "skin 2",
              "small intestine",   "smooth muscle", "soft tissue 1", "soft tissue 2",
              "spleen" ,  "stomach 1", "stomach 2",  "testis", "thyroid gland", "tonsil",
              "urinary bladder", "vagina", "N/A","hypothalamus", "retina",
              "pituitary gland", "choroid plexus")
  #cell types u don't use
  cell_type_diff = c("cells in cortex/medulla",       "cells in cuticle" ,
                     "cells in external root sheath",  "cells in internal root sheath",
                     "lactating glandular cells",  "sebaceous cells", "secretory cells",
                     "sweat ducts", "corneal epithelial cells",  "hyaloid membrane",
                     "lens epithelial cells", "lens fiber cells",  "cortical cells",
                     "ductal cells", "cells in dentate nucleus")

  filtered_hpa = hpaNormalTissue %>% dplyr::filter(Tissue %in% tissues &
                                                     !(Cell.type %in% cell_type_diff) &
                                                 Gene.name %in% G_list$hgnc_symbol)

  #data for the graph
  nodes = data.frame(id=unique((filtered_hpa$Gene.name)),
                     label = unique(filtered_hpa$Gene.name),
                     color = rep("yellow", length(unique(filtered_hpa$Gene.name))),
                     shadow = rep(TRUE, length(unique(filtered_hpa$Gene.name))),
                     color.border = "black",
                     shape = "circle")
  #change the color of the nodes in x or y tissue
  if(tissue != "all")
  {
    if(tissue %in% tissues)
    {
      #modify color in nodes so it colors the nodes the way u want
      #filters to nodes that contain tissue info
      filtered_hpa_color = filtered_hpa %>% dplyr::filter(Tissue == tissue)
      #l
      nodes$color[nodes$label %in% unique(filtered_hpa_color$Gene.name)] = "pink"
    }
    else
    {
     #throw an error telling that tissue does not exist
    }
  }

  edges <- data.frame(from= rep(gene, length(unique(filtered_hpa$Gene.name))),
                      to= unique(filtered_hpa$Gene.name),
                      color = "orange")
  #graph
  #vizualize graph
  print(visNetwork::visNetwork(nodes = nodes, edges = edges, width = "100%") %>%
                                                                      visNetwork::visIgraphLayout())
  #add option to change the color

}

removeSmallComponent = function(g)
{
  sub_gs <- igraph::components(g)$membership

  small_sub <- names(which(table(sub_gs) == 2))


  (rm_nodes <- names(which(sub_gs == small_sub)))

  g2 <- igraph::delete_vertices(g, rm_nodes)
  return(g2)
}





