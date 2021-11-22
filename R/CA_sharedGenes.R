CA_sharedGenes = function(D_i, D_j)
{
  if(D_i == D_j)
  {
    return(NA)
  }
  #these are the genes associated with PCOS
  pcos_associated_genes = gene_disease_associations %>% dplyr::filter(disease_merge == "PCOS")
  #genes associated with D_j
  j_associated_genes = gene_disease_associations %>% dplyr::filter(disease_merge == D_j)
  #genes associated
  i_associated_genes =gene_disease_associations %>% dplyr::filter(disease_merge == D_i)
  #genes associated with pcos and anemia
  G_i = intersect(pcos_associated_genes$geneSymbol, i_associated_genes$geneSymbol)
  #genes associated with PCOs and autoimmune
  G_j = intersect(pcos_associated_genes$geneSymbol, j_associated_genes$geneSymbol)
  #comorbidity score when PCOS,
  score = (length(intersect(G_i, G_j))/max(length(G_i), length(G_j)))*100
  return(score)
}
