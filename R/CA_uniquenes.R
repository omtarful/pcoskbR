
uniqueness = function(g_i)
{
  #D_T: total number of diseases in the gene-disease dataset
  D_T = length(unique(gene_disease_associations$disease_merge))
  #D_g_i: number of diseases associated with the ith gene
  D_g_i = gene_disease_associations %>% dplyr::filter(g_i == geneSymbol)
  D_g_i = length(D_g_i$geneSymbol)
  uniqueness_score = 1 - sqrt(D_g_i/D_T)
  return(uniqueness_score)
}

#if there are genes in common between D_i and D_j

getNCommonGenes = function(D_i, D_j)
{
  D_j_genes = (gene_disease_associations[which(D_j == gene_disease_associations$disease_merge),])$geneSymbol
  D_i_genes = (gene_disease_associations[which(D_i == gene_disease_associations$disease_merge),])$geneSymbol
  PCOS_genes = (gene_disease_associations[which("PCOS" == gene_disease_associations$disease_merge),])$geneSymbol
  ngenes = intersect(intersect(PCOS_genes, D_j_genes), intersect(PCOS_genes, D_i_genes))
  return(ngenes)
}

getUniquenessCScore = function(D_i, D_j)
{
  if(D_i == D_j){
    return(NA)
  }
  ngenes = getNCommonGenes(D_i, D_j)
  sum = 0
  for (g_i in ngenes)
  {
    sum = sum + uniqueness(g_i)
  }

  score = (sum/length(ngenes))*100
  return(score)
}


