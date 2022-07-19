library(rvest)
library(httr)
library(data.table)
library(qtl2convert)
library(ggrepel)
library(dplyr)
library(tidyverse)
library(reshape2)
library(stats)
library(ggplot2)

# function to scrape individual files from single URL
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

# function to retrieve summary table
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

# files with only one allele
read_sample_txt <- function(filename) {
  print(filename)
  file_df <- fread(filename)
  if(ncol(file_df) == 9) {
    df <- file_df
  } else {
    df <- NULL # do not output if file is incorrect
    print(paste0('file with problematic columns:', filename))
  }
  return(df)
}


# function to generate and convert files to 1tl2 input format
get_qtl2_input <- function(dir_name, sample_type, output_dir_name, summary_df){
  # list CONTAINER to store outputs
  file_output <- list()
  
  ### get all data from directory
  # read the config file
  config <- fread(paste0(root, '/annotations_config.csv'))
  config <- config %>% filter(array_type == sample_type)
  annot_file <- config$annot_file
  dict_file <- config$dict_file
  
  ### annotations
  dir <- "https://raw.githubusercontent.com/kbroman/MUGAarrays/main/UWisc/" # always the same
  annot_df <- read.csv(paste0(dir, annot_file))
  annot_dict <- read.csv(paste0(dir, dict_file))
  
  # direct to the data directory filepath
  data_dir <- file.path(root, dir_name)
  
  # all txt file from directory
  data_files <- dir(data_dir, pattern = '\\.txt$', full.names = TRUE)
  #summary_file <- dir(data_dir, pattern = '\\.csv$', full.names = TRUE)
  
  ### combine all files
  df_all <- rbindlist(lapply(data_files, read_sample_txt)) # use sample_id or original?
  
  # summary df
  sum_df <- summary_df %>% 
    select(ID, Sex, `% Het. Calls`, `% No Call`) %>% rename(sample_id = ID)
  sum_df$`% Het. Calls` <- as.numeric(gsub("%", "", sum_df$`% Het. Calls`))
  sum_df$`% No Call` <- as.numeric(gsub("%", "", sum_df$`% No Call`))
  
  # eliminate those that has no call > 10%
  df_all <- merge(df_all, sum_df, by = 'sample_id') %>% filter(`% No Call` < 10)
  
  ### genotype data
  geno_sub <- df_all %>% select(sample_id, snp_id, allele1_fwd, allele2_fwd) %>%
    mutate(gene_exp = paste(allele1_fwd, allele2_fwd, sep = '')) %>%
    select(sample_id, snp_id, gene_exp) %>% rename(marker = snp_id)
  df_geno <- dcast(geno_sub, marker~sample_id, value.var="gene_exp") # inbred - most are homozygous
  df_geno[is.na(df_geno)] <- "--"
  # get marker annotations from wisc
  # encode_genome function on geno and founder from qtl2 convert package - args: matrix of genotypes, matrixs of two alleles
  annot_encode_df <- annot_df %>% select(marker, chr, snp) %>%
    separate(snp, c("allele1", "allele2"), sep=cumsum(c(1)))
  # put alleles in order of marker as shown in geno dataframe
  geno_encode_annot <- merge(df_geno, annot_encode_df, by = 'marker', all.x = T) %>% select(allele1, allele2)
  geno_encode_annot <- find_unique_geno(geno_encode_annot, na.strings = c('NA', '-', ""))
  geno_encode_annot[is.na(geno_encode_annot)] <- '-'
  rownames(df_geno) <- df_geno$marker
  df_geno_encoded <- as.data.frame(encode_geno(df_geno[,-1], geno_encode_annot))
  df_geno_encoded$marker <- rownames(df_geno_encoded)
  # reorder columns
  df_geno_encoded <- df_geno_encoded %>% select(marker, everything())
  file_output[[1]] <- df_geno_encoded 
  print('genotype data done')
  
  ### gmap and pmap 
  # extract only columns needed from sample file
  map_df <- df_all %>% 
    select(snp_id, chromosome, position_bp) %>%
    rename(marker = snp_id, chr = chromosome) %>% unique()
  # extract columns from annotations
  annot_map <- annot_df %>% select(marker, chr, bp_mm10, cM_cox) %>% unique()
  # merge the two dataframes accordingly
  map_df <- merge(map_df, annot_map, by = c('marker')) %>% drop_na() # NAs are arrays to detect mutations
  # hence, drop NAs (for now)
  
  ## genetic map data (pos in cM)
  # drop NAs, only keep chromosomes that combine 
  df_gmap <- map_df %>% select(marker, chr.x, cM_cox) %>% rename(chr = chr.x, pos = cM_cox) %>% drop_na()
  file_output[[2]] <- df_gmap
  print('genetic map data done')
  ## physical map data (pos in Mbp)
  df_pmap <- map_df %>% select(marker, chr.x, bp_mm10) %>% rename(chr = chr.x, pos = bp_mm10) %>%
    mutate(pos = pos/1000000) %>% drop_na()
  ### check if gmap and pmap - intersect to make sure markers are consistent
  file_output[[3]] <- df_pmap
  print('physical map data done')
  
  ## phenotype data
  # simulate using rnorm/runiform
  df_pheno <- df_all %>% select(sample_id) %>% unique()
  df_pheno$WBC <- floor(runif(nrow(df_pheno), 1, 9)) # copying Dan's work from 2014
  df_pheno$NEUT <- floor(rnorm(nrow(df_pheno), 800, 500))
  file_output[[4]] <- df_pheno
  print('phenotype data done')
  
  # covariate file - map to sex
  #summary_df <- fread(summary_file)
  # sample_id, sex
  df_cov_all <- merge(sum_df, df_all %>% select(sample_id, chromosome, allele1_fwd, allele2_fwd), by = 'sample_id')

  # remove those with overall no calls more than 10%
  df_cov_all <- df_cov_all[!df_cov_all$`% No Call` > 10,] %>% select(sample_id, Sex, chromosome, allele1_fwd, allele2_fwd)
  df_cov_all_xy <- df_cov_all[df_cov_all$chromosome %in% c('X', 'Y'),]
  # count het_x and nc_y
  nc_y <- df_cov_all_xy[(df_cov_all_xy$allele1_fwd == '-') & (df_cov_all_xy$allele2_fwd == '-'),] %>% # both are '-', no call
    filter(chromosome == 'Y') %>%
    group_by(sample_id) %>% summarise(NC_Y = n())
  het_x <- df_cov_all_xy[(df_cov_all_xy$allele1_fwd != df_cov_all_xy$allele2_fwd),] %>% # heterozygous
    filter(chromosome == 'X') %>%
    group_by(sample_id) %>% summarise(Het_X = n())
  plot_all <- merge(nc_y, het_x, by = 'sample_id')
  # scatterplot to show results
  ggplot(plot_all) + aes(x = Het_X, y = NC_Y, label = sample_id) + geom_point() + 
    geom_text_repel(size=3, max.overlaps = 20) + ggtitle('SNP count - Heterozygous on X vs. No call on Y in each sample')
  # split genders
  females_list <- plot_all[plot_all$NC_Y > 75,]$sample_id
  males_list <- plot_all[plot_all$NC_Y < 75,]$sample_id
  # create new column
  df_cov_all$calc_sex <- NA
  # calculated gender
  df_cov_all$calc_sex[df_cov_all$sample_id %in% females_list] <- 'female'
  df_cov_all$calc_sex[df_cov_all$sample_id %in% males_list] <- 'male'
  # list of difference in calculated gender and original
  sex_diff <- unique(df_cov_all[(df_cov_all$Sex != df_cov_all$calc_sex),]$sample_id)
  print(paste0('Samples with different calculated gender than original entered: ', paste(sex_diff, collapse=', ' )))
  # select columns needed
  df_covar <- df_cov_all %>% select(sample_id, calc_sex) %>% rename(id = sample_id, Sex = calc_sex) %>% unique()
  file_output[[5]] <- df_covar
  print('phenotype covariate data done')
  
  ### founder genomes
  founders_dict <- fread(paste0(root, '/founder_strains_table.csv')) %>% rename(`Strain Name` = founder_strain)
  founder_strains <- unique(founders_dict$`Strain Name`) # get the names of founder strains
  # web scraping for founder sample 'UNC_Villena_GIGMUGV01_20141012_FinalReport'
  # set urls
  founder_url <- 'https://haploqa.jax.org/tag/UNC_Villena_GIGMUGV01_20141012_FinalReport.html'
  url_domain <- 'https://haploqa.jax.org'
  # get html and list of url to sample data
  html_table <- read_html(founder_url) %>% html_nodes("a") %>% html_attr("href")
  url_list <- paste0(url_domain, html_table[grepl('/sample', html_table)])
  iter_len <- length(url_list)
  # loop through urls to extract individual files
  founders_total = data.frame() # container
  inc = 0
  for (url in url_list) {
    file <- sample_individual_scrape(url, url_domain) # call function from another script
    inc = inc+1 # track progress
    print(paste0('working on file ', inc, '/', iter_len, ' ', file))
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
  rownames(test_foundergeno) <- test_foundergeno$marker
  # put the allele annotations in according order
  founder_encode_annot <- merge(annot_encode_df, test_foundergeno, by = 'marker') %>% select(allele1, allele2)
  df_founders_encoded <- as.data.frame(encode_geno(test_foundergeno[,-1], founder_encode_annot))
  df_founders_encoded$marker <- rownames(df_founders_encoded)
  df_founders_encoded <- df_founders_encoded %>% select(marker, everything())
  file_output[[6]] <- df_founders_encoded 
  print('founder geno data done')
  
  # cross info
  # join with summary table and CC file in karl's github
  # match the strain names with the cc cross info csv
  ci_temp <- read.csv('https://raw.githubusercontent.com/rqtl/qtl2data/main/CC/cc_crossinfo.csv', skip=3) %>% 
    separate(id, c("strain_id", "sec_id"), "/") %>% select(-sec_id)
  ci_sum <- summary_df %>% select(`Strain Name`, `ID`)
  ci_sum$strain_id <- str_extract(ci_sum$`Strain Name`, "CC...")
  ci_sum <- ci_sum %>% select(-`Strain Name`)
  df_crossinfo <- merge(ci_temp, ci_sum, by = 'strain_id', all.y = T) %>% 
    select(ID, everything()) %>%
    rename(id = ID) %>% drop_na() # read_cross does not allow NAs here
  file_output[[7]] <- df_crossinfo 
  print('cross info data done')
  
  return(file_output)
}
