#' Venn Analsis Tool
#' @description
#' Illustrates shared genes, ontologies or pathways between diseases as 6-way Venn Diagram. It takes a vector of up to 6 diseases and options.
#' @param disease_list Diseases for which comorbidity with PCOS will be analyzed. A possible list of diseases can be retrieved using the \code{listDiseases} function.
#' @param option It can either be "Genes", "Pathways" or "Ontologies".
#'
#' @return A Venn diagram.
#' @export
#'
#' @examples
#' \code{getVennAnalysis(disease_list = c("Tangier Disease", "Hypercholesterolemia", "Neuropathy" ), option = "Genes")}
getVennAnalysis = function(disease_list, option  = "Genes")
{
  #check if Venn Analysis
  if(missing(disease_list))
    stop("Argument 'disease_list' must be specified")
  if(!isInDataset(dataset = "Genes", disease_list = disease_list))
    stop("Argument 'disease_list' contains one or more diseases not found in dataset")
  if(length(disease_list) > 6)
    stop("Argument 'disease_list' accepts a vector of up to 6 diseases" )
  if(option != "Genes" & option != "Pathways" & option != "Ontologies")
    stop("Argument 'option' can only be \"Genes\", \"Pathways\" or \"Ontologies\".")
  #u need a list of every set
  optionList <- generateDiseaseList(disease_list, option)
  names(optionList) = disease_list
  #this is how u create a plot
  #this is how uturn them into triangles
  svg_url = "https://vnote-1251564393.cos.ap-chengdu.myqcloud.com/typora-img/triangles.svg"
  knitr::include_graphics(svg_url)

  vertex_coordinates <- list(c(-69277,-32868,135580,121186, 70900,199427),
                             c( 81988,-44426, 38444,206222,121044,165111),
                             c(203271,  9619, 39604, 82683, 84652,206669),
                             c(333561,225349, 61764, 76805, 38980,182461),
                             c(131886,385785, 38136,111491, 94208, 24690),
                             c(-60184,274046,142476, 39903,103276,183962))

  triangles <- lapply(vertex_coordinates, ggVennDiagram::triangle)


  position <- tibble::tribble(
    ~x,       ~y,
    -50000,     50000,
    60000,          0,
    160000,     20000,
    280000,    170000,
    140000,    300000,
    -20000,   270000
  )
  label_position = ggVennDiagram::label_position(position)

  shape = ggVennDiagram::VennPlotData(setEdge = triangles,
                       setLabel = label_position)

  venn =  ggVennDiagram::Venn(optionList)

  data =  ggVennDiagram::plotData_add_venn(plotData = shape, venn = venn)
  items <- ggVennDiagram::venn_region(data) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(text = stringr::str_wrap(paste0(.data$item, collapse = " "),
                                           width = 40)) %>%
    sf::st_as_sf()
  label_coord = sf::st_centroid(items$geometry) %>% sf::st_coordinates()
  p <- ggplot2::ggplot(items) +
    ggplot2::geom_sf(ggplot2::aes(fill=id), show.legend = FALSE) +
    ggplot2::geom_sf_text(ggplot2::aes_string(label = "name"),
                 data = data@setLabel,
                 inherit.aes = F) +
    ggplot2::geom_text(ggplot2::aes_string(label = "count", text = "text"),
              x = label_coord[,1],
              y = label_coord[,2],
              show.legend = FALSE) +
    ggplot2::theme_void()

  ax <- list(
    showline = FALSE
  )
  plotly::ggplotly(p, tooltip = c("text")) %>%
    plotly::layout(xaxis = ax, yaxis = ax)
}

#generates a list that
generateDiseaseList = function(disease_list, option = "Genes")
{
  optionList = list()
  if(option == "Genes")
  {
    for (disease in disease_list)
    {
      gene_vector = (gene_disease_associations %>% dplyr::filter(disease_merge == disease))$geneSymbol
      optionList[[disease]] = gene_vector
    }
  }
  else if(option == "Pathways")
  {
      disease_table = getDiseaseTable(disease_list)
      for (row in 1:nrow(disease_table))
      {
        genes = countGenes(disease_table[row,2])
        pathways = union(getPathways(genes, database = "KEGG")$Pathway, getPathways(genes, database = "Reactome")$Pathway)
        optionList[[disease_table[row,1]]] = pathways
      }
  }
  else if(option == "Ontologies")
  {
    for (disease in disease_list)
    {
      optionList[[disease]] = getDiseaseGO(disease)$`GO ID`
    }
  }
  else
  {
    #this throws an error
  }
  return(optionList)
}

getIntersections = function(disease_list, disease_sets)
{
  optionList <- generateDiseaseList(disease_list, option)
  names(optionList) = disease_list
  #this is how u create a plot
  #this is how uturn them into triangles
  svg_url = "https://vnote-1251564393.cos.ap-chengdu.myqcloud.com/typora-img/triangles.svg"
  knitr::include_graphics(svg_url)

  vertex_coordinates <- list(c(-69277,-32868,135580,121186, 70900,199427),
                             c( 81988,-44426, 38444,206222,121044,165111),
                             c(203271,  9619, 39604, 82683, 84652,206669),
                             c(333561,225349, 61764, 76805, 38980,182461),
                             c(131886,385785, 38136,111491, 94208, 24690),
                             c(-60184,274046,142476, 39903,103276,183962))

  triangles <- lapply(vertex_coordinates, ggVennDiagram::triangle)


  position <- tibble::tribble(
    ~x,       ~y,
    -50000,     50000,
    60000,          0,
    160000,     20000,
    280000,    170000,
    140000,    300000,
    -20000,   270000
  )
  label_position = ggVennDiagram::label_position(position)

  shape = ggVennDiagram::VennPlotData(setEdge = triangles,
                                      setLabel = label_position)

  venn =  ggVennDiagram::Venn(optionList)

  data =  ggVennDiagram::plotData_add_venn(plotData = shape, venn = venn)
  items <- gVennDiagram::venn_region(data) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(text = stringr::str_wrap(paste0(.data$item, collapse = " "),
                                           width = 40)) %>%
    sf::st_as_sf()
  label_coord = sf::st_centroid(items$geometry) %>% sf::st_coordinates()
  n_elem = NULL
  genes = NULL
  for (row in 1:length(items$name))
  {
    diseases = strsplit(items$name[row], split = "..", fixed = TRUE)[[1]]
    if(setequal(diseases, disease_sets))
    {
      n_elem = items$count[row]
      genes = strsplit(items$text[row], split = "~")[[1]]
      break
    }
  }
  cat(paste0("Count: ", n_elem))
  print()
  cat(paste0("Genes: ", genes))

}

