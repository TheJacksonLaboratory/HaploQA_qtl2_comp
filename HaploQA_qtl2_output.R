##### HaploQA data reformatting - Launch script
### used for qtl2 input
### annotation config saved in annotations_config.csv
### founder strain dictionary saved in founder_strains_table.csv

library(rstudioapi)
library(qtl2)
library(dplyr)
library(tidyverse)
library(httr)
library(rvest)
library(data.table)
library(qtl2convert)
library(ggrepel)
library(reshape2)
library(stats)
library(ggplot2)
library(stringr)
library(stringi)
library(lsa)

### BXD: cross_info - all 1s, choose one founder genome to be A. Make sure foundergeno and cross_info have consistent labels for genomes.
# geno for founders - http://haploqa-dev.jax.org/tag/UNC_Villena_GIGMUGV01_20141012_FinalReport.html
# pick a color for DBA
#
# 1. run BXD
# 2. dynamic cross_info instead of hard-coding
# 3. number of crossovers, visualize in some way
# 4. make cos_sim plots on high level (mean of a sample across entire genome)
# 5. have Shiny work for BXD as well
# 6. modularize and write documentations for functions
# 7. http://haploqa-dev.jax.org/tag/f2.html


# source the function script
root <- dirname(getSourceEditorContext()$path)
source(paste0(root,"/input_data_prep_functions.R"))

list_pheno <- c('WBC', 'NEUT')

do_results <- haplotype_reconstruction_pipeline('DO', list_pheno, qtl2_file_gen = F, samples_gen = F)

## count pct of missing values in sample geno

cc_results <- haplotype_reconstruction_pipeline('CC', list_pheno, qtl2_file_gen = F, samples_gen = F)

bxd_results <- haplotype_reconstruction_pipeline('BXD', list_pheno, qtl2_file_gen = F, samples_gen = F)

f2_results <- haplotype_reconstruction_pipeline('F2', list_pheno, qtl2_file_gen = F, samples_gen = F)


minimuga_results <- haplotype_reconstruction_pipeline('MiniMUGA', list_pheno, qtl2_file_gen = T, samples_gen = T)



geno_all_comp <- function(sample, sample_type, sample_results) {
  num_chr <- c((1:19),"X")
  print(sample)
#sample <- '35V'
#sample_type <- 'CC'
#sample_results <- cc_results
founder_lookup_table <- fread(file.path(root, 'founder_lookup_table.csv'))
founder_all_rev_lookup <- setNames(founder_lookup_table$founder_codes, founder_lookup_table$founder_id)

config <- fread(paste0(root, '/annotations_config.csv'))
config_sample <- config[config$array_type == sample_type]
qtl2_dir <- file.path(root, config_sample$qtl2_dir)
df_cov <- fread(file.path(qtl2_dir, 'test_covar.csv')) %>% rename('sample_id' = 'id')
data_dir <- paste0(root, '/', config_sample$data_dir)
summary_df <- fread(paste0(data_dir, '/', sample_type, '_summary.csv'))

haploqa_diplotype <- get_haplotypes(summary_df, data_dir)

## make dictionary
if (sample_type %in% c('CC', 'DO')) { ### change this to sample_type in CC or DO, they are special conditions
  founders_dict <- fread(paste0(root, '/founder_strains_table.csv'))
  founder_haplo_lookup <- setNames(founders_dict$letter, founders_dict$founder_strain)
} else {
  unique_haplotypes <- unique(haploqa_diplotype[,c(haplotype1, haplotype2)])
  founder_haplo_lookup <- setNames(LETTERS[seq(1, length(unique(unique_haplotypes)))], unique(unique_haplotypes))
}

raw_geno_df <- haploqa_diplotype %>% select(sample_id, snp_id, haplotype1, haplotype2, chromosome) %>% unique()
raw_geno_df[,c(3,4)] <- as.data.frame(apply(raw_geno_df[,c(3,4)], 2, function(x) founder_haplo_lookup[x]))
raw_geno_df <- merge(raw_geno_df, df_cov, by = 'sample_id')
raw_geno_df[(raw_geno_df$chromosome == 'X') & (raw_geno_df$Sex == 'male')]$haplotype2 <- 'Y'
raw_geno_df$haplotype <- paste(raw_geno_df$haplotype2, raw_geno_df$haplotype1, sep='')
raw_geno_df <- raw_geno_df %>% select(sample_id, snp_id, chromosome, haplotype)
raw_geno_df[,4] <- lapply(raw_geno_df[,4], function(col) {
  rev_col = stri_reverse(col)
  ifelse(rev_col %in% founder_all_rev_lookup, rev_col, col)
})


founder_all_lookup <- setNames(names(founder_all_rev_lookup), founder_all_rev_lookup) # limit the lookup table
founder_all_lookup <- founder_all_lookup[names(founder_all_lookup) %in% unique(raw_geno_df$haplotype)]
founder_all_lookup <- setNames(LETTERS[seq(1, length(unique(raw_geno_df$haplotype)))], seq(1, length(unique(raw_geno_df$haplotype))))

df_geno_all_chr <- list()

for (chr in num_chr) {
  #chr <- 'X'
  
  print(chr)
    ## raw geno
  raw_geno_acgt <- get_raw_geno(sample_type, chromosome = chr)
  raw_geno_chr_acgt <- raw_geno_acgt %>% filter(sample_id == sample)
  
  haploqa_mom <- as.data.frame(apply(as.data.frame(sample_results[['ph_geno_haploqa']][[chr]][,,1][rownames(sample_results[['ph_geno_haploqa']][[chr]][,,1]) == sample,]), 1, function(x) founder_all_lookup[x]))
  haploqa_mom$marker <- rownames(haploqa_mom)
  haploqa_mom <- haploqa_mom %>% rename(haploqa_mom = 1)
  haploqa_dad <- as.data.frame(apply(as.data.frame(sample_results[['ph_geno_haploqa']][[chr]][,,2][rownames(sample_results[['ph_geno_haploqa']][[chr]][,,2]) == sample,]), 1, function(x) founder_all_lookup[x]))
  haploqa_dad$marker <- rownames(haploqa_dad)
  haploqa_dad <- haploqa_dad %>% rename(haploqa_dad = 1)
  sample_geno <- merge(haploqa_mom, haploqa_dad, by = 'marker')
  sample_geno$haplo_diplotype <- paste(sample_geno$haploqa_mom, sample_geno$haploqa_dad, sep='')
  sample_geno <- sample_geno %>% select(marker, haplo_diplotype)
  
  
  haploqa_mom <- as.data.frame(apply(as.data.frame(sample_results[['ph_geno']][[chr]][,,1][rownames(sample_results[['ph_geno']][[chr]][,,1]) == sample,]), 1, function(x) founder_all_lookup[x]))
  haploqa_mom$marker <- rownames(haploqa_mom)
  haploqa_mom <- haploqa_mom %>% rename(haploqa_mom = 1)
  haploqa_dad <- as.data.frame(apply(as.data.frame(sample_results[['ph_geno']][[chr]][,,2][rownames(sample_results[['ph_geno']][[chr]][,,2]) == sample,]), 1, function(x) founder_all_lookup[x]))
  haploqa_dad$marker <- rownames(haploqa_dad)
  haploqa_dad <- haploqa_dad %>% rename(haploqa_dad = 1)
  sample_geno_qtl2 <- merge(haploqa_mom, haploqa_dad, by = 'marker')
  sample_geno_qtl2$qtl2_calls <- paste(sample_geno_qtl2$haploqa_mom, sample_geno_qtl2$haploqa_dad, sep='')
  sample_geno_qtl2 <- sample_geno_qtl2 %>% select(marker, qtl2_calls)
  
  
  
  df_geno_all_acgt <- merge(merge(raw_geno_chr_acgt, sample_geno, by = 'marker', sort = F), sample_geno_qtl2, by = 'marker', sort = F) %>% select(marker, sample_id, chr, pos, everything())
  df_geno_all_acgt$is_crossover <- NA
  
  ind <- which(df_geno_all_acgt$haplo_diplotype != df_geno_all_acgt$qtl2_calls)
  
  list_xo <- list()
  list <- split(ind, cumsum(c(1, diff(ind) != 1)))
  for (i in seq(1, length(list))) {
    list_xo <- c(list[[i]][1], list_xo)
  }
  list_xo <- sort(unlist(list_xo))
  
  rows <- unlist(lapply(list_xo, function(x) (x-20):(x+20)))
  rows <- rows[(rows > 0) & (rows < nrow(df_geno_all_acgt))]
  if (length(rows) > 0) {
    df_geno_all_acgt[list_xo,]$is_crossover <- 'True'
    df_geno_xo <- df_geno_all_acgt[rows,]
    
    df_geno_xo[is.na(df_geno_xo$is_crossover),]$is_crossover <- 'False'
    
    df_geno_xo$original_row <- rows # always 41 - 20 above, 20 below, and crossover itself
    df_geno_xo$crossover_group <- NA
    df_geno_xo$crossover_group <- as.character(df_geno_xo$crossover_group)
    
    iter <- 1
    for (i in which(diff(df_geno_xo$original_row) > 1)) {
     # i <- 41
      
      df_geno_xo[1:i][is.na(df_geno_xo[1:i]$crossover_group)]$crossover_group <- LETTERS[iter]
      iter <- iter+1
    }
    df_geno_xo$crossover_group[is.na(df_geno_xo$crossover_group)] <- LETTERS[iter]
    df_geno_xo$crossover_groups <- paste(df_geno_xo$crossover_group, df_geno_xo$chr, sep='')
    df_geno_xo <- df_geno_xo %>% select(-any_of(c("crossover_group")))
    
  } else {
    df_geno_xo <- data.frame()
  }
  
    
  df_geno_all_chr[[chr]] <- df_geno_xo
  }

  return(df_geno_all_chr)
}


 
sample_geno_1 <- rbindlist(geno_all_comp('35V', 'CC', cc_results))
write.csv(sample_geno_1, file.path(root, 'genocomp_acgt', '35V_cc_geno_comp.csv'))
sample_geno_2 <- rbindlist(geno_all_comp('35X', 'CC', cc_results))
write.csv(sample_geno_2, file.path(root, 'genocomp_acgt', '35X_cc_geno_comp.csv'))


sample_geno_3 <- rbindlist(geno_all_comp('8RS', 'F2', f2_results))
write.csv(sample_geno_3, file.path(root, 'genocomp_acgt', '8RS_f2_geno_comp.csv'))
sample_geno_4 <- rbindlist(geno_all_comp('8RT', 'F2', f2_results))
write.csv(sample_geno_4, file.path(root, 'genocomp_acgt', '8RT_f2_geno_comp.csv'))

sample_geno_5 <- rbindlist(geno_all_comp('5PK', 'BXD', bxd_results))
write.csv(sample_geno_5, file.path(root, 'genocomp_acgt', '5PK_bxd_geno_comp.csv'))
sample_geno_6 <- rbindlist(geno_all_comp('5PD', 'BXD', bxd_results))
write.csv(sample_geno_6, file.path(root, 'genocomp_acgt', '5PD_bxd_geno_comp.csv'))

sample_geno_7 <- rbindlist(geno_all_comp('6UY', 'DO', do_results))
write.csv(sample_geno_7, file.path(root, 'genocomp_acgt', '6UY_do_geno_comp.csv'))
sample_geno_8 <- rbindlist(geno_all_comp('6VE', 'DO', do_results))
write.csv(sample_geno_8, file.path(root, 'genocomp_acgt', '6VE_bxd_geno_comp.csv'))

### get all crossovers, 20 markers above and below
### run this for one sample, all chromosomes, add chr and pos information
## make into one file



#do_haplo_test <- list()



cos_sim_plot(do_results, rds_dir, results_dir, 'DO')




# lookup tables
founder_all_codes <- colnames(cc_results[['pr']]$`X`) # take the X chromosome - this one has everything
all_map_codes <- seq(1:length(founder_all_codes))
founder_all_lookup <- setNames(all_map_codes, founder_all_codes) # make maxmarg
founder_all_rev_lookup <- setNames(founder_all_codes, all_map_codes)

## location crossover comparison
pos_haploqa_cc <- locate_xo(cc_results[['ginf_haploqa']], cc_results[['map']])
pos_qtl2_cc <- cc_results[['pos']]
x_loc_diff <- location_xo_comp(pos_qtl2_cc, pos_haploqa_cc, num_chr)

## error matrix
err_comp(cc_results[['pr']], cc_results[['map']], num_chr, results_dir, sample_type)

## genocode comparison matrix
comp_matrix <- genocode_comp_matrix(cc_results[['map']], cc_results[['ginf']], cc_results[['ginf_haploqa']], founder_all_rev_lookup, num_chr, file_gen = T, sample_type)

### percent genomic difference
chr_pct <- get_geno_pct_diff(cc_results[['ginf']], cc_results[['ginf_haploqa']], cc_results[['summary']], num_chr, founder_all_rev_lookup)

# lookup tables
founder_all_codes <- colnames(do_results[['pr']]$`X`) # take the X chromosome - this one has everything
all_map_codes <- seq(1:length(founder_all_codes))
founder_all_lookup <- setNames(all_map_codes, founder_all_codes) # make maxmarg
founder_df <- data.frame(founder_id = all_map_codes, founder_codes = founder_all_codes)

## location crossover comparison
### number of crossovers on autosomes
pdf(file.path(root, "num_crossover_plots.pdf"))
xo_number_plot(do_results)
xo_number_plot(cc_results)
xo_number_plot(bxd_results)
xo_number_plot(f2_results)
dev.off()


## crossover distances
pdf(file.path(root, "crossover_distance_mean_plots.pdf"))
loc_xo_distance_plot(do_results)
loc_xo_distance_plot(cc_results)
loc_xo_distance_plot(bxd_results)
loc_xo_distance_plot(f2_results)
dev.off()


## error matrix
err_comp(do_results[['pr']], do_results[['map']], num_chr, results_dir, sample_type)

## genocode comparison matrix
comp_matrix <- genocode_comp_matrix(do_results[['map']], do_results[['ginf']], do_results[['ginf_haploqa']], founder_all_rev_lookup, num_chr, file_gen = T, sample_type)


### percent genomic difference
chr_pct <- get_geno_pct_diff(do_results[['ginf']], do_results[['ginf_haploqa']], do_results[['summary']], num_chr, founder_all_rev_lookup)


## pipeline for one
sample_haplotype_reconstruction <- function(sample_type, sample_name) {
  sample_type <- 'MiniMUGA'
  sample_name <- 'JXN'
  ## results directory
  results_dir <- file.path(root, 'results')
  dir.create(results_dir, showWarnings = FALSE) 
  
  num_chr <- c((1:19),"X")
  
  ### config file
  config <- fread(paste0(root, '/annotations_config.csv'))
  
  ### Environment
  config_sample <- config[config$array_type == sample_type]
  marker_type <- config_sample$marker_type
  founders_list <- unlist(strsplit(config_sample$founders_list, ", "))
  n_founders <- length(founders_list)
  ngen <- config_sample$ngen
  sample_url <- config_sample$url
  
  # data output directory
  data_dir <- file.path(root, config_sample$data_dir)
  dir.create(data_dir, showWarnings = FALSE) 
  
  ## create a data directory for qtl2 input data
  qtl2_dir <- file.path(root, config_sample$qtl2_dir) # name of desired output folder
  # create if folder not exist
  dir.create(qtl2_dir, showWarnings = FALSE)
  
  # filepaths to save any rds file
  rds_dir <- file.path(results_dir, 'RDS')
  dir.create(rds_dir, showWarnings = FALSE)
  
  # annotation file
  annot_file <- config_sample$annot_file
  
  ### control file
  control_fp <- paste0(qtl2_dir, '/test.json')
  if (file.exists(control_fp) == FALSE) {
    file.create(control_fp)
  }
  
  ## summary file
  summary_df_fp <- paste0(data_dir, '/', sample_type, '_summary.csv')
  
  ## list of samples to exclude, if any
  exclude_list <- unlist(strsplit(config_sample$exclude_samples, ", "))
  ## text of control file
  ### alleles need to be consistent with genail
  control_file <- paste0('{
  "description": "HaploQA data - qtl2 test run",
  "crosstype": "genail', n_founders, '",
  "sep": ",",
  "na.strings": ["-", "NA"],
  "comment.char": "#",
  "geno": "test_geno.csv",
  "founder_geno": "test_foundergeno.csv",
  "gmap": "test_gmap.csv",
  "pmap": "test_pmap.csv",
  "pheno": "test_pheno.csv",
  "covar": "test_covar.csv",
  "genotypes": {
    "A": 1,
    "H": 2,
    "B": 3
  },
  "x_chr": "X",
  "alleles": [', substring(gsub("[()]", "", toString(list(LETTERS[seq(1, n_founders)]))),2), '], 
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
}')
  
  ###############################################################################
  ### Part 1 - get results from HaploQA
  # set main URL domain
  url_domain <- 'http://haploqa-dev.jax.org/' 
  
  #### implementations
  # summary table
  if (!file.exists(summary_df_fp)) {
    print(paste0('Working on summary file:'))
    #### html file prep
    # read html file
    haploqa_html <- read_html(sample_url)
    # extract information from html file
    html_temp <- haploqa_html %>% html_nodes("a") %>% html_attr("href")
    url_list <- paste0(url_domain, html_temp[grepl('/sample', html_temp)])
    summary_df <- sample_summary_scrape(haploqa_html, url_list, marker_type)
    print('Writing to directory')
    write.csv(summary_df, paste0(data_dir, '/', sample_type, '_summary.csv'), row.names = FALSE)
    
  } else {
    summary_df <- fread(summary_df_fp)
    summary_df <- summary_df[summary_df$ID == sample_name]
  }
  
  #url_ind <- url_ind[!url_ind %in% exclude_list]
  
  if(samples_gen == T) {
    #### html file prep
    # read html file
    haploqa_html <- read_html(sample_url)
    # extract information from html file
    html_temp <- haploqa_html %>% html_nodes("a") %>% html_attr("href")
    url_list <- paste0(url_domain, html_temp[grepl('/sample', html_temp)])
    
    # list of urls to generate samples from
    url_ind <- unique(summary_df$`Sample Filepath`)
    
    # individual samples
    inc = 0
    for (url in url_ind) {
      inc = inc + 1 # increment
      file <- sample_individual_scrape(url, url_domain)
      print(paste0('Working on file ', inc, '/', length(url_ind), ': ', file))
      sample_df_save <- as.data.frame(content(GET(file)))
      print('Writing to directory')
      file_name <- unlist(strsplit(file, '/'))[6]
      GET(file, write_disk(paste0(data_dir, '/', file_name), overwrite = TRUE), show_col_types = FALSE)
    }
  }
  
  file_output <- get_qtl2_input(data_dir, sample_type, annot_file, qtl2_dir, summary_df, list_pheno, ngen, founders_list, marker_type, exclude_list = NULL)
  # unpack
  df_geno <- file_output[[1]]
  df_gmap <- file_output[[2]]
  df_gmap <- sort_chr(df_gmap, c((1:19),"X")) # put chromosomes in order
  df_pmap <- file_output[[3]]
  df_pmap <- sort_chr(df_pmap, c((1:19),"X"))
  df_pheno <- file_output[[4]]
  df_covar <- file_output[[5]]
  df_foundergeno <- file_output[[6]] # with strain id metadata
  df_crossinfo <- file_output[[7]]
  founder_haplo_lookup <- file_output[[8]]
}



  