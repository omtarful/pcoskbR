
#scrapes the GO:IDs from website for a given gene
getGeneGO = function(gene)
{
  css_selector = "table+ table div+ table td:nth-child(2) , table+ table div+ table td:nth-child(1)"
  #url for the given gene to scrape
  url = paste0("http://pcoskb.bicnirrh.res.in/gene_pcos.php?gene=", gene)
  pcoskb = rvest::session(url)
  pcoskb = rvest::read_html(pcoskb)
  name = pcoskb %>%
    rvest::html_nodes(css_selector) %>% # returns multiple nodes that contain the css_selector
    rvest::html_text()
  name = trimws(name) #trims leading and trailing spaces
  GO_df = NULL
  indices <- seq(1,by=2, len=length(name))
  for (i in indices)
  {
    GO_df = rbind(GO_df, name[i:(i+1)])
  }
  colnames(GO_df) = as.vector(GO_df[1,])
  GO_df = as.data.frame(na.omit(GO_df[-1,]))
  return(GO_df)
}



getCA_GOScore = function(D_i, D_j)
{
  if(D_i == D_j)
  {
    return(NA)
  }
  GO_i_df = getDiseaseGO(D_i)
  GO_j_df = getDiseaseGO(D_j)
  #this is the score using just Cellular component
  GO_i = GO_i_df$`GO ID`
  GO_j = GO_j_df$`GO ID`
  #this is the score using not Celullar component
  GO_i = (GO_i_df %>% dplyr::filter(Ontology == "Biological process" | Ontology == "Molecular function"))$`GO ID`
  GO_j = (GO_j_df %>% dplyr::filter(Ontology == "Biological process" | Ontology == "Molecular function"))$`GO ID`
  score_1 = (length(intersect(GO_i, GO_j))/length(union(GO_i, GO_j)))*100

  return(score_1)
}

#takes a vector of genes (of a given disease) and returns all GO terms
#fix it so it doesn't use cellular components
getDiseaseGO = function(disease_name)
{
  if(disease_name == "PCOS")
  {
    return(pcos_ontology)
  }
  GO = gene_disease_associations %>% dplyr::filter(disease_name == disease_merge)
  gene_vector = GO$geneSymbol
  GO = data.frame(GO_id = NA, Ontology=NA)
  colnames(GO) = c("GO ID", "Ontology")
  for (gene in gene_vector)
  {
    GO = dplyr::union(getGeneGO(gene), GO) #this loops though every gene and gets the list
  }
  return(unique(na.omit(GO)))
}





