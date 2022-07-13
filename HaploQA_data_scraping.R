library(dplyr)
library(tidyr)
library(rstudioapi)
library(logr)

#### Set toggles
generate_summary <- TRUE # set to TRUE to generate summary table
output_summary <- FALSE # set to TRUE to output summary table to data directory
generate_sample_individuals <- TRUE # set to TRUE to generate individual sample data
output_samples <- FALSE # set to TRUE to output sample data to directory


#### Environment set-up
# get the directory this file is stored in
root <- dirname(getSourceEditorContext()$path)

# source functions used
source(paste0(root,"/web_scrape_functions.R"))

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
  #url <- url_list[1]
  inc = inc + 1 # increment
  if (generate_sample_individuals == TRUE) {
    file <- sample_individual_scrape(url)
    print(paste0('Working on file ', inc, '/66: ', file))
    sample_df_save <- as.data.frame(content(GET(file)))
    if (output_samples == TRUE) {
      print('Writing to directory')
      file_name <- unlist(strsplit(file, '/'))[5]
      GET(file, write_disk(paste0(data_dir, '/', file_name)))
    }
  }
}

