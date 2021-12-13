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
#takes the output of the pathway analysis function and returns a network where the nodes are the pathways
#function works regardless of database or dataset
viewNetwork = function(pathways)
{
  if(missing(pathways))
    stop("Argument 'disease_list' must be specified")
  #check if pathways have the correct format
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



getPathway_mirna = function(gene_pathway_table, pathway_id, pcos_associated_mirna, dataset)
{
  gene_associated_w_pathways = NULL
  if("Reactome" == dataset)
  {
    hsa_ncbi2reactome = ncbi2reactome %>% dplyr::filter(grepl(pattern = "R-HSA", x = V2))
    filtered_hsa_ncbi2reactome = hsa_ncbi2reactome %>% dplyr::filter(V2 == pathwayId)
    gene_associated_w_pathways = filtered_hsa_ncbi2reactome$V1
  }
  else if("KEGG" == dataset)
  {
    gene_associated_w_pathways = gene_pathway_table %>% dplyr::filter(pathway_id == `associated_pathway`)
    gene_associated_w_pathways = gene_associated_w_pathways$entrez_id
  }
  else
  {

  }

  #search for them in mirTar
  miRNA_data = human_data %>% dplyr::filter(`Target Gene (Entrez Gene ID)` %in% gene_associated_w_pathways)
  #search for the miRNA inside pcosk
  miRNA_names = miRNA_data$miRNA
  #convert miRNA to precursrs
  mature_to_precursor = miRBaseConverter::miRNA_MatureToPrecursor(miRNA_names)
  miRNA_names = mature_to_precursor$Precursor
  result = pcos_associated_mirna[pcos_associated_mirna$mirbase_id %in% miRNA_names,1]
  return(result)
}

#gets pathway using
createMirnaPathway = function(disease_target_genes, database, disease_mirnas)
{
  #convert geneSymbol to mirBase id
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  pcos_mirna_table =biomaRt::getBM(attributes = c("hgnc_symbol", "mirbase_id"), filters = "hgnc_symbol", values = unique(mirna_disease_associations$`Gene Symbol`),  mart = mart)
  pcos_mirna_table[16,2] = "hsa-mir-320a"
  pthwys = NULL
  if(database == "KEGG")
  {
    kegg_gene_pathway_table = NULL
    for (id in disease_target_genes)
    {
      geneInfo = tryCatch( KEGGREST::keggGet(paste0("hsa:", id)), error = function(e) NULL)
      if(!is.null(geneInfo))
      {
        entrez_id = geneInfo[[1]][["ENTRY"]]
        gene_pathways = names(geneInfo[[1]][["PATHWAY"]])
        gene_pathway_table = data.frame(entrez_id=rep(entrez_id, length(gene_pathways)),
                                        associated_pathway = gene_pathways)
        kegg_gene_pathway_table = rbind(kegg_gene_pathway_table, gene_pathway_table)
      }
    }
    pthwys = getPathwayUsingId(kegg_gene_pathway_table$associated_pathway)

    pathway_table = data.frame(`Pathway` = pthwys$Pathway,
                               `Pathway genes associated with PCOS` = pthwys$`Pathway ID`,
                               `Pathway genes associated with PCOS and the selected diseases` = rep(NA, length( pthwys$Pathway)),
                               `Hypergeometric probability` = rep(NA, length( pthwys$Pathway)))
    #gwt pathway mirna
    for (row in 1:nrow(pathway_table))
    {
      pathway_id = pathway_table[row,2]
      pathway_table[row,2] = paste0(getPathway_mirna(gene_pathway_table = kegg_gene_pathway_table, pathway_id, pcos_associated_mirna = pcos_mirna_table, dataset = "KEGG"), collapse = ", ")
      pathway_table[row,3] = intersect(countGenes(pathway_table[row,2]), disease_mirnas)
    }
    #this is the gene population, it's N in the formula

  }
  else if(database == "Reactome")
  {
    #filter human pathways
    hsa_ncbi2reactome = ncbi2reactome %>% dplyr::filter(grepl(pattern = "R-HSA", x = V2))
    #filter to target genes
    filtered_hsa_ncbi2reactome = hsa_ncbi2reactome %>% dplyr::filter(V1 %in% disease_target_genes)
    #filter to pathways uniquely in PCOSKB
    pthwys = getPathwayUsingId(filtered_hsa_ncbi2reactome$V2)
    for (row in 1:nrow(pthwys)) {
      #find entrez gene associated with each pathway
      pathwayId = pthwys[row,2]
      pathway_genes = filtered_hsa_ncbi2reactome %>% dplyr::filter(V2 == pathwayId)
      pathway_gene_ids = pathway_genes$V1
      #find all miRNA associated with the pathway
      pathway_mirna = human_data %>% dplyr::filter(`Target Gene (Entrez Gene ID)` %in% pathway_gene_ids)
      #convert it to precursor
      miRNAs = unique(pathway_mirna$miRNA)
      #the extra empty string is added to the vector is for when there is a single miRNA,
      #to avoid error
      if(length(miRNAs) == 1)
      {
        miRNAs= append(miRNAs, "")
      }
      mature_to_precursor = miRBaseConverter::miRNA_MatureToPrecursor(miRNAs)
      miRNA_names = mature_to_precursor$Precursor
      result = pcos_mirna_table[pcos_mirna_table$mirbase_id %in% miRNA_names,1]
      pthwys[row,2] = paste0(result[!is.na(result)], collapse = ", ")

    }
    return(pthwys)

  }
  else
  {

  }

  #filters kegg pathways that are solely in pcoskbr
  out_pathways = pthwys
  #get miRNA for each pathway


  gene_population = disease_mirnas
  pcos_associated_mirna = mirna_disease_associations %>% dplyr::filter(Disease == "PCOS")
  colnames(pathway_table) = c("Pathway name",  "Pathway genes associated with PCOS", "Pathway genes associated with PCOS and the selected diseases", "Hypergeometric probability")
  N = union(countGenes(pathway_table$`Gene (Gene Symbol)`), gene_disease_associations$geneSymbol) #population of genes in pcos and pathway dataset
  n = gene_population #sample size = number of gnes used in input

  for (i in 1:nrow(pathway_table))
  {
    K = unique(countGenes(pathway_table$`Pathway genes associated with PCOS`[i])) #number of pcos genes in ith pathway
    k = intersect(K, gene_population) #number of inpiut genes in sample overlapping iwth genes in ith pathay
    pathway_table[i,"Hypergeometric probability"] =dhyper(x = length(k), m = length(K), n = length(setdiff(N, K)), k = length(n))
  }
  pathway_table = pathway_table[,c(1,3,2,4)]

}

getPathwayUsingId = function(kegg_pathway_ids)
{
  query = ""
  for (pathway_id in kegg_pathway_ids) {
    query = paste0(query, " ", paste0("[ID]{", pathway_id, "}"))
  }
  url = "http://pcoskb.bicnirrh.res.in/advsrh.php"
  advcd_search  = rvest::session(url)
  html_node = rvest::read_html(advcd_search)

  search_form = html_node %>% html_form()
  search_form = search_form[[1]]
  search = search_form %>% html_form_set(qry = query, make = "Associated Pathways", type = "ID")
  resp <- html_form_submit(search)
  tables = read_html(resp) %>% html_nodes("p+ table td") %>% html_text()
  name = tables
  n_col = 4
  n_row = length(name)/n_col #calculates the number of rows in the table
  table_section = NULL
  #loop splits row of text into a vector
  for (i in 1:n_row)
  {
    start_index = 1 + n_col*(i - 1) # calculates start index  of the ith row in the name vector
    end_index = start_index + n_col - 1 #calculates the end_index of the ith row in the name vector
    name = trimws(name)
    row = name[start_index:end_index] #slices the elements to make the ith row
    table_section = rbind(table_section, row) #binds the (i-1)th row to the ith row
  }
  table_section = as.data.frame(table_section)
  col_names = table_section[1,]
  colnames(table_section) = col_names
  table_section = table_section[-1,]
}

#' Pathway enrichment Analysis
#' @description
#' For selected diseases, this tool gives a list of enriched pathways.
#' @usage
#' \code{PathwayAnalysis(disease_list, dataset = c("miRNAs", "Genes"), database = c("KEGG", "Reactome"))}
#' @param disease_list Diseases for which comorbidity with PCOS will be analyzed. A possible list of diseases can be retrieved using the \code{listDiseases} function.
#' @param dataset Dataset. It can either be "miRNAs" or "Genes".
#' @param database Database from which pathways will be extracted. It can either be "KEGG" or "Reactome".
#'
#' @author Omar Hassoun
#'
#' @return \code{data.frame} of enriched pathways.
#' @export
#'
#' @examples
#' \code{PathwayAnalysis(disease_list = c("Tangier Disease", "Hypoalphalipoproteinemia", "Neuropathy"),
#'          dataset = "Genes",
#'          database = "KEGG")}

PathwayAnalysis = function(disease_list, dataset, database)
{
  #check if dataset argument and disease list argument is correct

  if(missing(disease_list))
    stop("Argument 'disease_list' must be specified")
  if(missing(dataset))
    stop("Argument 'dataset' must be specified")
  if(database != "KEGG" & database != "Reactome")
    stop("Argument database can either be \"KEGG\" or \"Reactome\"")
  if(!isInDataset(dataset = dataset, disease_list))
    stop("Argument 'disease_list' contains one or more diseases not found in dataset")
  if(dataset != "miRNAs" & dataset != "Genes")
    stop("Argument 'dataset' must be specified")
  #getting gene population
  disease_table = getDiseaseTable(disease_list, dataset = dataset)
  gene_population = countGenes(disease_table$`Gene Symbol`) #outputs gene symbols
  #map miRNAs to genes using mirTar
  disease_mirnas = NULL #miRNAs to be used as input in n variable
  if(dataset == "miRNAs")
  {
    disease_mirnas = gene_population
    #convert miRNAs to mirBase ids using biomaRt
    mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
    disease_mirna_table = biomaRt::getBM(attributes = c("hgnc_symbol", "mirbase_id"), filters = "hgnc_symbol", values = unique(disease_mirnas),  mart = mart)
    disease_mirnas_mirTar = disease_mirna_table$mirbase_id
    disease_mirnas_mirTar = paste0(disease_mirnas_mirTar, "-?")
    disease_data = human_data %>% dplyr::filter(grepl(paste(disease_mirnas_mirTar, collapse = "|"), x = miRNA, ignore.case = TRUE))
    #get pathways
    disease_target_genes = unique(disease_data$`Target Gene (Entrez Gene ID)`)
    gene_population = disease_target_genes #outputs entrez gene ids
  }
  #step 2: Map genes to pathways depending on database
  if(dataset == "Genes")
    #outputs a data.frame of mapped pathways with their KEGG id, Source and Database
    #pathwa-gene association table from pcoskbR
    pathway = getPathways(gene_population, database = database)
  else
  {
    #makes the whole analysis
    pathway = createMirnaPathway(disease_target_genes = gene_population, database = "KEGG", disease_mirnas)
    return(pathway)
  }
  #assign name for each column and calculate score
  colnames(pathway) =  c("Pathway name",  "Pathway genes associated with PCOS", "Pathway genes associated with PCOS and the selected diseases", "Hypergeometric probability")
  for (row in 1:nrow(pathway)) {
    pathway[row,2] = paste0(countGenes(pathway[row,3]) , collapse = ", ") #assigns genes to pathways
    pathway[row,3] = paste0(intersect(countGenes(pathway[row,2]), gene_population), collapse = ", ") #assigns pcos and selected diseaes
  }
  pathway = pathway[,c(1,3,2,4)]

  #vakues u need to calculate p value
  N = unique(gene_disease_associations$geneSymbol) #population of genes
  n = gene_population #sample size = number of gnes used in input
  #assign proper value to each column
  for (i in 1:nrow(pathway))
  {
    K = unique(countGenes(pathway$`Pathway genes associated with PCOS`[i])) #number of pcos genes in ith pathway
    k = intersect(K, gene_population) #number of inpiut genes in sample overlapping iwth genes in ith pathay
    pathway[i,"Hypergeometric probability"] =dhyper(x = length(k), m = length(K), n = length(setdiff(N, K)), k = length(n))
  }
  return(pathway)
}
