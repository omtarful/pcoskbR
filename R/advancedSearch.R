advancedSearch = function(gene_list)
{
  query = ""
  for (gene in gene_list) {
    query = paste0(query, " ", paste0("[Gene]{", gene, "}"))
  }
  url = "http://pcoskb.bicnirrh.res.in/advsrh.php"
  advcd_search  = rvest::session(url)
  html_node = rvest::read_html(advcd_search)

  search_form = html_node %>% rvest::html_form()
  search_form = search_form[[1]]
  search = search_form %>% html_form_set(qry = query, make = "Associated Pathways")

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
  return(table_section)
}
