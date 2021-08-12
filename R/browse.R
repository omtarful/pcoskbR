#How does it do it?
#it uses the rvest package to scrape and submit the filters in the website

#' Browse: enables users to query the page for genes, miRNAs, SNPs, diseases, pathways, and ontology associated with PCOS.
#'
#' @param tab A string representing which page will it retrieve the page from (i.e genes, miRNA, SNPs...).
#' @param filter A string to filter tables using certain subsets.
#' @param search_qry A string used to filter the tables according to a keyword.
#'
#' @return It returns a data frame with results fuiltered according to parameters.
#' @export
#'
#' @examples
browse = function(page, filter = NULL, search_qry = "")
{
  #this is the url of the page u scrape data from
  url = "http://pcoskb.bicnirrh.res.in/"
  #Initialize start_page and end_page variables
  start_page = 0 #it always starts at 0
  end_page = 0 #it changes according to the filter
  n_col = NULL #this is the number of columns in the website
  tag = NULL #this is the tag u need to get the table rows
  end_tag = NULL #this is the tag u need to get the end_page
  title_tag = NULL #ihis is the tag u need to get the column names
  #it first checks which page u selected and uses that to construct the url eleme
  if(page == "Genes" | page == "miRNA")
  {
    n_col = 6
    if(page == "Genes")
    {
      tab_choice = "gene"
    } else
    {
      tab_choice = "mirna"
    }
    par = "&limitt="
    tag = "#form6 tr+ tr td+ td"
    end_tag = "table+ table strong"
    title_tag =  "#form6 strong"
  }
  else if(page == "Diseases" | page == "miRNA Diseases")
  {
    n_col = 3
    if(page == "Diseases")
    {
      tab_choice = "disease"
    }
    else
    {
      tab_choice = "diseasem"
    }
    par = "&par="
    tag = ".expander td"
    end_tag = "p strong"
    title_tag = "table~ table+ table strong"
  }
  else if(page == "SNPs")
  {
    n_col = 7
    tab_choice = "snp"
    par = ""
    filter = ""
    tag = "td td td tr+ tr td"
    end_tag = "table+ table strong"
    title_tag = ".style2 td :nth-child(1)"
  }
  else if(page == "Pathways")
  {
    n_col = 4
    tab_choice = "pathway_pcos"
    par =  "&limitt="
    tag = "table+ table tr+ tr td"
    end_tag = "table:nth-child(3) strong"
    title_tag = "table~ table+ table strong"
  }
  else if(page == "Ontologies")
  {
    n_col = 4
    tab_choice = "go_d"
    par  = "&limitt="
    tag = ".expander td"
    end_tag = "tr+ tr strong"
    title_tag = "table+ table strong"
  }
  else
  {
    #this is where u handle pages that do not exist
  }
  #creates the url u need to scrape the page
  url = construct_url(tab_choice = tab_choice, last_part_of_par = par, filter = filter, search_qry = search_qry)
  end_page = getEnd_page(url, tag = end_tag) #get the last page in the website
  table = NULL #initializes table
  #it goes through every page in the website, scrapes the rows
  #and concatinates the the current page with the last till end_page is reached
  for (i in start_page:(end_page))
  {
    section = unique(getTableSection(url = url, n_col = n_col, tag = tag, n_page = i)) #scrapes the table in the i page
    table = rbind(table, section) #binds the i page with the i - 1
  }
  colnames(table) = getColumnNames(url = url, title_tag = title_tag) #names the columns in the table
  return(table)
}
#gets the section inside the table
# it takes the already constructed url, the number of columns, the tag to scrape the rows and the page
getTableSection = function(url, n_col, tag, n_page)
{
  url = sub(0, n_page, url) #replaces the 0 in the url for whatever the value of page number is
  table_section = NULL #initilizes the table section u need to construct data frame
  name = getText(url = url, css_selector = tag) #this gets the text u need to use to create the rows
  a <- grep("snp", url) #this catches if the url is in the snp page
  if(length(a) != 0 && a == 1) #if the url is using the snp page it enters the if statement and modifies name
  {
    del = seq(0, length(name), by = 8) #creates a sequence of indexes to delete from name
    del = del[-1] #skips the first index
    name = name[-del] #deletes those indexes from name
    name = str_replace_all(name, "[\r\n]" , "") #removes carriage return
    # the reason wh the indexes are deleetd is that snp is not really working properly
  }
  n_row = length(name)/n_col #calculates the number of rows in the table
  #loop splits row of text into a vector
  for (i in 1:n_row)
  {
    start_index = 1 + n_col*(i - 1) # calculates start index  of the ith row in the name vector
    end_index = start_index + n_col - 1 #calculates the end_index of the ith row in the name vector
    row = name[start_index:end_index] #slices the elements to make the ith row
    table_section = rbind(table_section, row) #binds the (i-1)th row to the ith row
  }
  table_section = as.data.frame(table_section) #converts table into data frame
  return(table_section)
}
#get column names
getColumnNames = function(url, title_tag)
{
  #add function to remove duplicates
  return(as.vector(unique(getText(url = url, css_selector = title_tag))))
}
#This functions returns the end of a page
getEnd_page = function(url, tag)
{
  name = getText(url = url, css_selector = tag) #gets the text to get end_page
  a = grep("go_d", url)
  if(a == 1 && length(a) != 0) #checks if it's the ontologies page and modifies it to work
  {
    name = name[10]
  }
  end_page = strsplit(name, split = "\\s+")[[1]] #gets end_page
  n_elements = length(strsplit(name, split = "\\s+")[[1]]) #gets the number of eleemnts in the vector
  if(a == 1 && length(a) != 0)
  {
    n_elements = n_elements - 1 #this is another mod so the ontologies work
  }
  end_page = as.numeric(end_page[n_elements-1])
  return(ceiling(end_page/10)-1)
}

#getText takes 3 arguments and returns either a string vector or a string
# url: this is the already concatinated url
# css_selector: this is the css_selector that get the part u need
# number (optional) : this is needed if u want to get the end_page,
#it's the index of the vector name,it contains the text that u need to get the page
getText = function(url, css_selector, number = NULL)
{
  pcoskb = session(url)
  pcoskb = read_html(pcoskb)
  name = pcoskb %>%
    html_nodes(css_selector) %>% # returns multiple nodes that contain the css_selector
    html_text()
  #this is the text u need to select a position in the vector, u didn't use it but u won't delete it cz it doesn't harm
  if(!is.null(number))
  {
    name = name[number]
  }
  return(name)
}

#this function takes the elements in an url and appends them
# tab_choice: this is the last part of the page
# last_part_of_par: it's the last parameter
# page: it's the page number
construct_url = function(tab_choice, last_part_of_par, filter, n_page = 0, search_qry = "")
{
  filter = sub(" ", "%20", filter) #replaces the pages in par for %20 so it works in url
  url = "http://pcoskb.bicnirrh.res.in/"
  par =  paste("&char=&qry=", search_qry, last_part_of_par, filter, sep = "")
  div = paste(tab_choice, ".php", sep ="") #concatinates the page with the .php
  #this concatinates an url depending on the table
  url = paste(url, div, "?page=", n_page, par, sep = "") #concatinates everything
  return(url)
}


