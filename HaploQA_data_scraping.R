library(rvest)
library(httr)
library(dplyr)
library(tidyr)
library(rstudioapi)
library(logr)

#### Set toggles
generate_summary <- TRUE # set to TRUE to generate summary table
output_summary <- TRUE # set to TRUE to output summary table to data directory
generate_sample_individuals <- TRUE # set to TRUE to generate individual sample data
output_samples <- TRUE # set to TRUE to output sample data to directory


#### Environment set-up
# get the directory this file is stored in
root <- dirname(getSourceEditorContext()$path)
# data output directory
data_dir <- file.path(root, 'haploqa_collab_cross') # name of desired output folder
# create if folder not exist
dir.create(data_dir, showWarnings = FALSE) 
# set main URL domain
url_domain <- 'https://haploqa.jax.org' 
# target sample URL
sample_url <- 'https://haploqa.jax.org/tag/Collaborative%20Cross.html' # change to the URL of sample

#### html file prep
# read html file
haploqa_html <- read_html(sample_url)

# extract information from html file
html_temp <- haploqa_html %>% html_nodes("a") %>% html_attr("href")
url_list <- paste0(url_domain, html_temp[grepl('/sample', html_temp)])

#### data retrieval functions
## main summary table
sample_summary_scrape <- function(html_file, url_list) {
  # convert into table
  html_table <- html_file %>% html_nodes(".table") %>% html_table()
  sum_table <- html_table[[1]][-1] # remove first column, as it's blank
  sum_table <- as.data.frame(sum_table, row.names = F) # convert into dataframe
  # clean up ID columns
  sum_table$`ID (Secondary IDs)` <- gsub("\\s+", "", sum_table$`ID (Secondary IDs)`) # remove spaces
  # split IDs and secondary IDs
  sum_table <- sum_table %>%
    separate(`ID (Secondary IDs)`, c("ID", "Secondary IDs"), "\\(")
  sum_table$`Secondary IDs` <- gsub("\\)", "", sum_table$`Secondary IDs`) #remove parentheses
  sum_table$`Sample Filepath` <- url_list
   
  return(sum_table)
}

## individual sample datasets
sample_individual_scrape <- function(url) {
  # read html of website with sample data
  html_sample <- read_html(url)
  # get class 'btn'
  sample_list <- html_sample %>% 
    html_nodes(".btn") %>% html_attr("href")
  # only retrieving SNP files for now
  fp <- sample_list[grepl('snp', sample_list)]
  # full file name to get from website
  file <- paste0(url_domain, fp)
  
  return(file)
}


#### implementations
# summary table
if (generate_summary == TRUE) {
  print(paste0('Working on summary file:'))
  sum_df <- sample_summary_scrape(haploqa_html, url_list)
  if (output_summary == TRUE) {
    print('Writing to directory')
    write.csv(sum_df, paste0(data_dir, '/collaborative_cross_summary.csv'), row.names = FALSE)
  }
}

# individual samples
inc = 0
for (url in url_list) {
  inc = inc + 1 # increment
  if (generate_sample_individuals == TRUE) {
    file <- sample_individual_scrape(url)
    print(paste0('Working on file ', inc, '/66: ', file))
    if (output_samples == TRUE) {
      print('Writing to directory')
      file_name <- unlist(strsplit(file, '/'))[5]
      GET(file, write_disk(paste0(data_dir, '/', file_name)))
    }
  }
}

