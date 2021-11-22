mirna_diseases = browse(page = "miRNA Diseases")
mirna_disease_associations = NULL
for (row in 1:nrow(mirna_diseases))
{
  disease_row = NULL
  disease_mirnas = countGenes(mirna_diseases[row,2])
  if(length(disease_mirnas) > 1)
  {
    disease_row = data.frame(Disease = rep(mirna_diseases[row,1], length(disease_mirnas)), `Gene Symbol`= disease_mirnas)
    colnames(disease_row) = c("Disease", "Gene Symbol")
  }
  else
  {
    disease_row = mirna_diseases[row,c(1,2)]
  }
  mirna_disease_associations = rbind(mirna_disease_associations, disease_row)
}
sysdata_filenames <- load("R/sysdata.rda")
save(list = c(sysdata_filenames, "mirna_disease_associations"), file = "R/sysdata.rda")
