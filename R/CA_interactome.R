getCA_interactomeScore = function(D_i, D_j)
{
  Dij_genes = (gene_disease_associations %>% dplyr::filter(D_j == disease_merge | D_j == disease_merge))$geneSymbol
  Dij_network_table <- STRINGdb::string_db$map(as.data.frame(Dij_genes), "Dij_genes", removeUnmappedRows = TRUE )
  Dij_subgraph = STRINGdb::string_db$get_subnetwork( Dij_network_table$STRING_id) #generates an igraph object of the gene population
  Dij = igraph::mean_distance(igraph::simplify(Dij_subgraph))
  #anemia
  Dii_genes = (gene_disease_associations %>% dplyr::filter(disease_merge == D_i))$geneSymbol
  Dii_network_table <- STRINGdb::string_db$map(as.data.frame(Dii_genes), "Dii_genes", removeUnmappedRows = TRUE )
  Dii_subgraph = STRINGdb::string_db$get_subnetwork(Dii_network_table$STRING_id) #generates
  Dii = igraph::mean_distance(igraph::simplify(Dii_subgraph))
  #Kawasaki disease
  Djj_genes = (gene_disease_associations %>% dplyr::filter(disease_merge == D_j))$geneSymbol
  Djj_network_table <- STRINGdb::string_db$map(as.data.frame(Djj_genes), "Djj_genes", removeUnmappedRows = TRUE )
  Djj_subgraph = STRINGdb::string_db$get_subnetwork(Djj_network_table$STRING_id) #
  Djj = igraph::mean_distance(igraph::simplify(Djj_subgraph))


  return((Dij - (Dii+Dij)/2))
}

