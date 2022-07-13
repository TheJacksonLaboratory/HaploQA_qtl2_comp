##### HaploQA data reformatting script
### used for qtl2 input
### annotation config saved in annotations_config.csv
### founder strain dictionary saved in founder_strains_table.csv

library(data.table)
library(rstudioapi)
library(qtl2)
library(dplyr)
library(tidyverse)
library(reshape2)
library(stats)
library(rjson)
library(ggplot2)


#### change the below params
dir_name <- 'haploqa_collab_cross' # enter directory name
sample_type <- 'Collaborative Cross' # Collaborative Cross/MiniMUGA/GigaMUGA
output_dir_name <- 'cc_qtl2_test'
file_gen = TRUE

# source functions used
root <- dirname(getSourceEditorContext()$path)
source(paste0(root,"/web_scrape_functions.R"))


## create a data directory for qtl2 input data
qtl2_dir <- file.path(root, output_dir_name) # name of desired output folder
# create if folder not exist
dir.create(qtl2_dir, showWarnings = FALSE)


# function to generate files
get_qtl2_input <- function(dir_name, sample_type, output_dir_name){
  # list to store outputs
  file_output <- list()
  
  ### get all data from directory
  config <- fread(paste0(root, '/annotations_config.csv')) %>% filter(array_type == sample_type)
  annot_file <- config$annot_file
  dict_file <- config$dict_file
  
  # direct to the data directory filepath
  data_dir <- file.path(root, dir_name)
  
  # all txt file from directory
  data_files <- dir(data_dir, pattern = '\\.txt$', full.names = TRUE)
  summary_file <- dir(data_dir, pattern = '\\.csv$', full.names = TRUE)
  
  ### combine all files
  df_all <- rbindlist(lapply(data_files, fread)) # use sample_id or original?
  
  ### genotype data
  geno_sub <- df_all %>% select(sample_id, snp_id, allele1_fwd, allele2_fwd) %>%
    mutate(gene_exp = paste(allele1_fwd, allele2_fwd, sep = '')) %>%
    select(sample_id, snp_id, gene_exp)
  df_geno <- dcast(geno_sub, snp_id~sample_id, value.var="gene_exp") # inbred - most are homozygous
  #write.csv(df_geno, paste0(qtl2_dir, '/test_geno.csv'), row.names = F)
  file_output[[1]] <- df_geno
  
  ### gmap and pmap 
  ### annotations
  dir <- "https://raw.githubusercontent.com/kbroman/MUGAarrays/main/UWisc/" # always the same
  annot_df <- read.csv(paste0(dir, annot_file))
  annot_dict <- read.csv(paste0(dir, dict_file))
  # extract only columns needed from sample file
  map_df <- df_all %>% 
    select(snp_id, chromosome, position_bp) %>%
    rename(marker = snp_id, chr = chromosome) %>% unique()
  # extract columns from annotations
  annot_map <- annot_df %>% select(marker, chr, bp_mm10, cM_cox) %>% unique()
  # merge the two dataframes accordingly
  map_all <- merge(map_df, annot_map, by = c('marker')) # NAs are arrays to detect mutations
  # hence, drop NAs (for now)
  
  ## genetic map data (pos in cM)
  # drop NAs, only keep chromosomes that combine 
  df_gmap <- map_all %>% select(marker, chr.x, cM_cox) %>% rename(chr = chr.x, pos = cM_cox) %>% drop_na()
  file_output[[2]] <- df_gmap
  ## physical map data (pos in Mbp)
  df_pmap <- map_all %>% select(marker, chr.x, bp_mm10) %>% rename(chr = chr.x, pos = bp_mm10) %>%
    mutate(pos = pos/1000000) %>% drop_na()
  ### check if gmap and pmap - intersect to make sure markers are consistent
  file_output[[3]] <- df_pmap
  
  ## phenotype data
  # simulate using rnorm/runiform
  df_pheno <- df_all %>% select(sample_id) %>% unique()
  df_pheno$WBC <- floor(runif(nrow(df_pheno), 1, 9)) # copying Dan's work from 2014
  df_pheno$NEUT <- floor(rnorm(nrow(df_pheno), 800, 500))
  file_output[[4]] <- df_pheno
  
  # covariate file - map to sex
  # sample_id, sex
  sum_df <- fread(summary_file) %>% 
    select(ID, Sex, `% Het. Calls`, `% No Call`) %>% rename(sample_id = ID)
  df_cov_sub <- df_all %>% select(sample_id, chromosome)
  df_cov_all <- merge(sum_df, df_cov_sub, by = 'sample_id')
  df_cov_all$`% Het. Calls` <- as.numeric(gsub("%", "", df_cov_all$`% Het. Calls`))
  df_cov_all$`% No Call` <- as.numeric(gsub("%", "", df_cov_all$`% No Call`))
  # new column
  df_cov_all$calc_sex <- NA
  # set cutoffs
  het_cutoff <- as.numeric(quantile(df_cov_all$`% Het. Calls`, 0.4))
  nocall_cutoff_lower <- as.numeric(quantile(df_cov_all$`% No Call`, 0.4))
  nocall_cutoff_higher <- as.numeric(quantile(df_cov_all$`% No Call`, 0.8))
  # females
  df_cov_all$calc_sex[(df_cov_all$`% Het. Calls` >= het_cutoff) & (df_cov_all$`% No Call` >= nocall_cutoff_higher)] <- 'female'
  # males
  df_cov_all$calc_sex[(df_cov_all$`% Het. Calls` <= het_cutoff) & (df_cov_all$`% No Call` <= nocall_cutoff_lower)] <- 'male'
  #df_cov_all[!is.na(df_cov_all$calc_sex),]
  
  df_covar <- df_cov_all %>% select(sample_id, Sex) %>% rename(id = sample_id) %>% unique()
  file_output[[5]] <- df_covar
  
  ### founder genomes
  founders_dict <- fread(paste0(root, '/founder_strains_table.csv')) %>% rename(`Strain Name` = founder_strain)
  founder_strains <- unique(founders_dict$`Strain Name`)
  # web scraping for founder sample 'UNC_Villena_GIGMUGV01_20141012_FinalReport'
  # set urls
  founder_url <- 'https://haploqa.jax.org/tag/UNC_Villena_GIGMUGV01_20141012_FinalReport.html'
  url_domain <- 'https://haploqa.jax.org' 
  # get html and list of url to sample data
  html_table <- read_html(founder_url) %>% html_nodes("a") %>% html_attr("href")
  url_list <- paste0(url_domain, html_table[grepl('/sample', html_table)])
  # loop through urls to extract individual files
  founders_total = data.frame() # container
  for (url in url_list) {
    file <- sample_individual_scrape(url) # call function from another script
    df_temp <- as.data.frame(content(GET(file)))
    founders_total <- rbind(founders_total,df_temp)
  }
  
  # founder strains all have secondary ID with '_'
  founders_df <- founders_total[grep("_", founders_total$original_sample_id), ] %>% 
    mutate(gene_exp = paste(allele1_fwd, allele2_fwd, sep = '')) %>%
    select('sample_id', 'snp_id', 'gene_exp')
  
  sum_df <- sample_summary_scrape(read_html(founder_url), url_list) %>% 
    select('ID', 'Strain Name') %>% rename(sample_id = ID)
  founders_test <- merge(founders_df, sum_df, by = 'sample_id')
  
  founders_test <- merge(founders_test, founders_dict, by = 'Strain Name')
  test_foundergeno <- founders_test %>% select(snp_id, letter, gene_exp) %>% 
    rename(strain = letter, marker = snp_id)
  test_foundergeno <- dcast(test_foundergeno, marker~strain, value.var="gene_exp")
  file_output[[6]] <- test_foundergeno 
  
  # control file
  # see if read_cross2 function works
  ## skip insert
  # calc_genoprob
  
  return(file_output)
}

### 1. overall no-calls can not be more than 10%
### no-calls on Y chrom for sex - fill in unknowns, and validate the known ones.

# implement function
file_output <- get_qtl2_input(dir_name, sample_type, output_dir_name)

# unpack
if (file_gen == TRUE) {
  df_geno <- file_output[[1]]
  write.csv(df_geno, paste0(qtl2_dir, '/test_geno.csv'), row.names = F)
  df_gmap <- file_output[[2]]
  write.csv(df_gmap, paste0(qtl2_dir, '/test_gmap.csv'), row.names = F)
  df_pmap <- file_output[[3]]
  write.csv(df_pmap, paste0(qtl2_dir, '/test_pmap.csv'), row.names = F)
  df_pmap <- file_output[[4]]
  write.csv(df_pmap, paste0(qtl2_dir, '/test_pheno.csv'), row.names = F)
  df_pmap <- file_output[[5]]
  write.csv(df_pmap, paste0(qtl2_dir, '/test_covar.csv'), row.names = F)
  df_pmap <- file_output[[6]]
  write.csv(df_pmap, paste0(qtl2_dir, '/test_foundergeno.csv'), row.names = F)
}

### view some files
minimuga_dir <- "/Users/linc/Documents/GitHub/HaploQA_qtl2_comp/MiniMUGA_qtl2_test"
cc_dir <- "/Users/linc/Documents/GitHub/HaploQA_qtl2_comp/cc_qtl2_test"
file_type <- 'pmap'
cc_sum <- fread(paste0(data_dir, '/collaborative_cross_summary.csv'))

minimuga_test <- fread(paste0(minimuga_dir, '/test_', file_type,'.csv'))
cc_test <- fread(paste0(cc_dir, '/test_', file_type,'.csv'))
head(minimuga_test)
head(cc_test)
