### contains all functions necessary to retrieve initial input files from haploqoa
### function directory
## 1. sample_summary_scrape - get metadata table for samples from haploqa website
## 2. sample_individual_scrape - get haplotype calls/input individual files 
## 3. read_sample_txt - function to read txt and perform some pre-processing

library(rstudioapi)
root <- dirname(getSourceEditorContext()$path)


# function to retrieve summary table from haploqa website
# basically takes the table shown on main page as-is
# @param html_file (html_document) - html of main page sample website on HaploQA, output of read_html(url)
# @param url_list (list) - url to main page of individuals in said sample, bound on summary table as metadata
#
# @return sum_table (data.frame) - table as shown on HaploQA sample website
# columns of sum_table: ID (sample IDs), secondary ID, Haplotype Candidate (T/F), Strain Name, Sex, % Het Calls, % Hom Calls, % No Call, % Concordance, Sample filepath
sample_summary_scrape <- function(html_file, url_list, marker_type) {
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
  if (!missing(marker_type)) {
    sum_table <- sum_table[sum_table$Platform == marker_type,]
  }
  return(sum_table)
}

# function to retrieve individuals in a given HaploQA sample
# @param url (string) - url of main page of each individual 
# @param url_domain (string) - almost always 'http://haploqa-dev.jax.org/' 
#
# @return file (string) - url to direct download of the file associated with individual
sample_individual_scrape <- function(url, url_domain) {
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

# function to read in the downloaded individual txt files
# need to check whether the columns are consistent across all files
# @param filepath (string) - file path the downloaded txt files are saved in
#
# @return df (data.frame) - individual files as data frame
# columns of df: sample_id,	original_sample_id, snp_id,	chromosome,	position_bp, allele1_fwd,	allele2_fwd, haplotype1, haplotype2
read_sample_txt <- function(filepath) {
  file_df <- fread(filepath)
  if(ncol(file_df) == 9) {
    df <- file_df
  } else {
    df <- NULL # do not output if file is incorrect
    ### print a log/warning
    print(paste0('skipped file with unaligned columns:', filepath))
  }
  return(df)
}