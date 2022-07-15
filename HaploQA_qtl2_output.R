##### HaploQA data reformatting script
### used for qtl2 input
### annotation config saved in annotations_config.csv
### founder strain dictionary saved in founder_strains_table.csv

library(data.table)
library(rstudioapi)
library(qtl2)
library(dplyr)
library(tidyverse)
library(httr)
library(jsonlite)

# source the function script
root <- dirname(getSourceEditorContext()$path)
source(paste0(root,"/web_scrape_functions.R"))

#### change the below params
dir_name <- 'haploqa_collab_cross' # enter directory name
sample_type <- 'Collaborative Cross' # Collaborative Cross/MiniMUGA/GigaMUGA
output_dir_name <- 'cc_qtl2_test'

#### Set toggles
generate_summary <- TRUE # set to TRUE to generate summary table
output_summary <- TRUE # set to TRUE to output summary table to data directory
generate_sample_individuals <- TRUE # set to TRUE to generate individual sample data
output_samples <- TRUE 
qtl2_file_gen = TRUE


### Part 1 - get results from HaploQA
# data output directory
data_dir <- file.path(root, dir_name)
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
    file <- sample_individual_scrape(url, url_domain)
    print(paste0('Working on file ', inc, '/66: ', file))
    sample_df_save <- as.data.frame(content(GET(file)))
    if (output_samples == TRUE) {
      print('Writing to directory')
      file_name <- unlist(strsplit(file, '/'))[5]
      GET(file, write_disk(paste0(data_dir, '/', file_name)))
    }
  }
}


### Part 2 - convert results to qtl2 input format
## create a data directory for qtl2 input data
qtl2_dir <- file.path(root, output_dir_name) # name of desired output folder
# create if folder not exist
dir.create(qtl2_dir, showWarnings = FALSE)

# implement function
file_output <- get_qtl2_input(dir_name, sample_type, output_dir_name, sum_df)

# unpack
if (qtl2_file_gen == TRUE) {
  df_geno <- file_output[[1]]
  write.csv(df_geno, paste0(qtl2_dir, '/test_geno.csv'), row.names = F)
  df_gmap <- file_output[[2]]
  write.csv(df_gmap, paste0(qtl2_dir, '/test_gmap.csv'), row.names = F)
  df_pmap <- file_output[[3]]
  write.csv(df_pmap, paste0(qtl2_dir, '/test_pmap.csv'), row.names = F)
  df_pheno <- file_output[[4]]
  write.csv(df_pheno, paste0(qtl2_dir, '/test_pheno.csv'), row.names = F)
  df_covar <- file_output[[5]]
  write.csv(df_covar, paste0(qtl2_dir, '/test_covar.csv'), row.names = F)
  df_foundergeno <- file_output[[6]]
  write.csv(df_foundergeno, paste0(qtl2_dir, '/test_foundergeno.csv'), row.names = F)
  df_crossinfo <- file_output[[7]]
  # temporary fix - NA can't be here
  df_crossinfo <- df_crossinfo %>% select(-ID) %>% unique() %>% drop_na()
  write.csv(df_crossinfo, paste0(qtl2_dir, '/test_crossinfo.csv'), row.names = F)
}

### control file
control_fp <- paste0(qtl2_dir, '/test.json')
file.create(control_fp)

control_file <- '{
  "description": "HaploQA data - qtl2 test run",
  "crosstype": "risib8",
  "sep": ",",
  "na.strings": ["-", "NA"],
  "comment.char": "#",
  "geno": "test_geno.csv",
  "founder_geno": "test_foundergeno.csv",
  "gmap": "test_gmap.csv",
  "pmap": "test_pmap.csv",
  "pheno": "test_pheno.csv",
  "covar": "test_covar.csv",
  "x_chr": "X",
  "genotypes": {
    "A": 1,
    "H": 2,
    "B": 3
  },
  "geno_transposed": true,
  "founder_geno_transposed": true,
  "sex": {
    "covar": "Sex",
    "female": "female",
    "male": "male"
  },
  "cross_info": {
    "file": "test_crossinfo.csv"
  }
}'

writeLines(control_file, control_fp)

cc_test <- read_cross2(control_fp)

