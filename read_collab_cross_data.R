library(data.table)
library(rstudioapi)
library(qtl2)
library(dplyr)
library(tidyverse)
library(reshape2)

## params
root <- dirname(getSourceEditorContext()$path)
# minimuga - mini_uwisc_dict_v3.csv (annotation), mini_uwisc_dict_v3.csv (dict)
# gigamuga/cc - gm_uwisc_v1.csv (annotation), gm_uwisc_dict_v1.csv (dict)
dir_name <- 'haploqa_MiniMUGA' # enter directory name
sample_type <- 'MiniMUGA' # Collaborative Cross/MiniMUGA/GigaMUGA
output_dir_name <- 'MiniMUGA_qtl2_test'
file_gen = TRUE

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
  dir <- "https://raw.githubusercontent.com/kbroman/MUGAarrays/main/UWisc/"
  annot_df <- read.csv(paste0(dir, annot_file))
  annot_dict <- read.csv(paste0(dir, dict_file))
  # extract only columns needed from sample file
  map_df <- df_all %>% 
    select(snp_id, chromosome, position_bp) %>%
    rename(marker = snp_id, chr = chromosome) %>%
    unique()
  # extract columns from annotations
  annot_map <- annot_df %>% select(marker, chr, bp_mm10, cM_cox) %>% unique()
  # merge the two dataframes accordingly
  map_all <- merge(map_df, annot_map, by = c('marker'))
  
  ## genetic map data (pos in cM)
  df_gmap <- map_all %>% select(marker, chr.x, cM_cox) %>% rename(chr = chr.x, pos = cM_cox)
  file_output[[2]] <- df_gmap
  ## physical map data (pos in Mbp)
  df_pmap <- map_all %>% select(marker, chr.x, bp_mm10) %>% rename(chr = chr.x, pos = bp_mm10) %>%
    mutate(pos = pos/1000000)
  file_output[[3]] <- df_pmap
  
  return(file_output)
}

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
}


