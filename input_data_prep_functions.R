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
### add a small document above each function, argument descriptions (ex. url -> text/str/char, 'url that points to the sample')
### description for returned item too - ex. dataframe/string/vector? and what it is
# do this for all functions
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
  sum_table <- sum_table[sum_table$Platform == 'GigaMUGA',]
  return(sum_table)
}

# files with only one allele
read_sample_txt <- function(filename) {
  file_df <- fread(filename)
  if(ncol(file_df) == 9) {
    df <- file_df
  } else {
    df <- NULL # do not output if file is incorrect
    ### print a log/warning
    print(paste0('skipped file with unaligned columns:', filename))
  }
  return(df)
}

# individual function
# one function for geno
qtl2_geno <- function(df, annot_encode_df) {
  ### genotype data
  geno_sub <- df %>% select(sample_id, snp_id, allele1_fwd, allele2_fwd) %>%
    mutate(gene_exp = paste(allele1_fwd, allele2_fwd, sep = '')) %>%
    select(sample_id, snp_id, gene_exp) %>% rename(marker = snp_id)
  df_geno <- dcast(geno_sub, marker~sample_id, value.var="gene_exp") # inbred - most are homozygous
  # merge with annotations
  df_geno <- df_geno %>% filter(marker %in% unique(annot_encode_df$marker)) ### replace with a stopif assert line here
  # put alleles in order of marker as shown in geno dataframe
  geno_encode_annot <- merge(df_geno, annot_encode_df, by = 'marker', all.x = T) %>% select(allele1, allele2) # let it break if there's NAs
  rownames(df_geno) <- df_geno$marker
  df_geno_encoded <- as.data.frame(encode_geno(df_geno[,-1], geno_encode_annot))
  df_geno_encoded$marker <- rownames(df_geno_encoded)
  # reorder columns
  df_geno_encoded <- df_geno_encoded %>% select(marker, everything())

  return(df_geno_encoded)
}

### function for gmap and pmap
qtl2_maps <- function(df, annot_df) {
  df_maps <- list()
  map_df <- df %>% 
    select(snp_id, chromosome, position_bp) %>%
    rename(marker = snp_id, chr = chromosome) %>% unique()
  # extract columns from annotations
  annot_map <- annot_df %>% select(marker, chr, bp_mm10, cM_cox) %>% unique()
  # merge the two dataframes 
  map_df <- merge(map_df, annot_map, by = c('marker')) %>% drop_na() # NAs are arrays to detect mutations
  # hence, drop NAs (for now)
  
  ## genetic map data (pos in cM)
  df_gmap <- map_df %>% select(marker, chr.x, cM_cox) %>% rename(chr = chr.x, pos = cM_cox) %>% drop_na()
  df_maps[[1]] <- df_gmap
  ## physical map data (pos in Mbp)
  df_pmap <- map_df %>% select(marker, chr.x, bp_mm10) %>% rename(chr = chr.x, pos = bp_mm10) %>%
    mutate(pos = pos/1000000) %>% drop_na()
  df_maps[[2]] <- df_pmap
  
  return(df_maps)
}

### function to simulate phenotype data
qtl2_pheno <- function(list_pheno, df) {
  df <- df %>% select(sample_id) %>% unique()
  for (i in seq(1:length(list_pheno))) {
    print('phenotypes:')
    print(list_pheno[i])
    df[[as.character(list_pheno[i])]] <- floor(runif(nrow(df), 1, 500))
  }
  return(df)
}

### covariate data
qtl2_cov <- function(df_unfiltered, split_level, df) {
  split_level <- as.numeric(split_level) # make sure it's a number
  ### if 'x' and 'y' both here, merge and use scatterplot
  if (any(c('Y', 'y') %in% df_unfiltered$chromosome) & any(c('X', 'x') %in% df_unfiltered$chromosome)) {
    df_cov_all_xy <- df_unfiltered %>% select(sample_id, chromosome, allele1_fwd, allele2_fwd) %>% 
      filter(sample_id %in% unique(df_all$sample_id)) %>%
      filter(chromosome %in% c('X', 'Y'))
    # count het_x and nc_y
    nc_y <- df_cov_all_xy[(df_cov_all_xy$allele1_fwd == '-') & (df_cov_all_xy$allele2_fwd == '-'),] %>% # both are '-', no call
      filter(chromosome == 'Y') %>%
      group_by(sample_id) %>% summarise(NC_Y = n())
    het_x <- df_cov_all_xy[(df_cov_all_xy$allele1_fwd != df_cov_all_xy$allele2_fwd),] %>% # heterozygous
      filter(chromosome == 'X') %>%
      group_by(sample_id) %>% summarise(Het_X = n())
    plot_all <- merge(nc_y, het_x, by = 'sample_id', all = T)
    # scatterplot to show results
    plot <- ggplot(plot_all) + aes(x = Het_X, y = NC_Y, label = sample_id) + geom_point() + 
      geom_text_repel(size=3, max.overlaps = 20) + ggtitle('SNP count - Heterozygous on X vs. No call on Y in each sample')
    print('please see plot for gender splitting and adjust accordingly')
    print(plot)
    # split genders
    females_list <- plot_all[plot_all$NC_Y > split_level,]$sample_id
    males_list <- plot_all[plot_all$NC_Y < split_level,]$sample_id
    # create new column
    df_cov_all <- df %>% select(sample_id, Sex)
    df_cov_all$calc_sex <- NA
    # calculated gender
    df_cov_all$calc_sex[df_cov_all$sample_id %in% females_list] <- 'female'
    df_cov_all$calc_sex[df_cov_all$sample_id %in% males_list] <- 'male'
    # list of difference in calculated gender and original
    sex_diff <- unique(df_cov_all[(df_cov_all$Sex != df_cov_all$calc_sex),]$sample_id)
    print(paste0('Samples with different calculated gender than original entered: ', paste(sex_diff, collapse=', ' )))
    # select columns needed
    df_covar <- df_cov_all %>% select(sample_id, calc_sex) %>% rename(id = sample_id, Sex = calc_sex) %>% unique()
  }
  
  if(!any(c('Y', 'y') %in% df_unfiltered$chromosome) & any(c('X', 'x') %in% df_unfiltered$chromosome)) {
    # if only X is present
    df_cov_all_x <- df_unfiltered %>% select(sample_id, chromosome, allele1_fwd, allele2_fwd) %>% 
      filter(sample_id %in% unique(df_all$sample_id)) %>%
      filter(chromosome %in% c('X'))
    # count het_x only
    het_x <- df_cov_all_x[(df_cov_all_x$allele1_fwd != df_cov_all_x$allele2_fwd),] %>% # heterozygous
      filter(chromosome == 'X') %>%
      group_by(sample_id) %>% summarise(Het_X = n()) %>% as.data.frame()
    # density plot to show one variable
    plot <- ggplot(het_x, aes(x = Het_X)) + geom_density(aes(y = ..count..), fill = "lightgray") + 
      geom_vline(aes(xintercept = mean(Het_X)), linetype = "dashed") + annotate(x=mean(het_x$Het_X), y=+Inf, label = round(mean(het_x$Het_X), 4), vjust = 1, geom = 'label')
    print('please see plot for gender splitting and adjust accordingly')
    print(plot)
    # split genders
    females_list <- het_x[het_x$Het_X > split_level,]$sample_id
    males_list <- het_x[het_x$Het_X < split_level,]$sample_id
    # create new column
    df_cov_all <- df %>% select(sample_id, Sex)
    df_cov_all$calc_sex <- NA
    # calculated gender
    df_cov_all$calc_sex[df_cov_all$sample_id %in% females_list] <- 'female'
    df_cov_all$calc_sex[df_cov_all$sample_id %in% males_list] <- 'male'
    # list of difference in calculated gender and original
    sex_diff <- unique(df_cov_all[(df_cov_all$Sex != df_cov_all$calc_sex),]$sample_id)
    print(paste0('Samples with different calculated gender than original entered: ', paste(sex_diff, collapse=', ' )))
    # select columns needed
    df_covar <- df_cov_all %>% select(sample_id, calc_sex) %>% rename(id = sample_id, Sex = calc_sex) %>% unique()
  } 
  
  if (any(c('Y', 'y') %in% df_unfiltered$chromosome) & !any(c('X', 'x') %in% df_unfiltered$chromosome)) {
    # if only Y is present
    df_cov_all_y <- df_unfiltered %>% select(sample_id, chromosome, allele1_fwd, allele2_fwd) %>% 
      filter(sample_id %in% unique(df_all$sample_id)) %>%
      filter(chromosome %in% c('Y'))
    # count no call in y only
    nc_y <- df_cov_all_xy[(df_cov_all_xy$allele1_fwd == '-') & (df_cov_all_xy$allele2_fwd == '-'),] %>% # both are '-', no call
      filter(chromosome == 'Y') %>%
      group_by(sample_id) %>% summarise(NC_Y = n()) %>% as.data.frame()
    # density plot to show one variable
    plot <- ggplot(nc_y, aes(x = NC_Y)) + geom_density(aes(y = ..count..), fill = "lightgray") + 
      geom_vline(aes(xintercept = mean(NC_Y)), linetype = "dashed") + annotate(x=mean(nc_y$NC_Y), y=+Inf, label = round(mean(nc_y$NC_Y), 4), vjust = 1, geom = 'label')
    print('please see plot for gender splitting and adjust accordingly')
    print(plot)
    # split genders
    females_list <- nc_y[nc_y$NC_Y > split_level,]$sample_id
    males_list <- nc_y[nc_y$NC_Y < split_level,]$sample_id
    # create new column
    df_cov_all <- df %>% select(sample_id, Sex)
    df_cov_all$calc_sex <- NA
    # calculated gender
    df_cov_all$calc_sex[df_cov_all$sample_id %in% females_list] <- 'female'
    df_cov_all$calc_sex[df_cov_all$sample_id %in% males_list] <- 'male'
    # list of difference in calculated gender and original
    sex_diff <- unique(df_cov_all[(df_cov_all$Sex != df_cov_all$calc_sex),]$sample_id)
    print(paste0('Samples with different calculated gender than original entered: ', paste(sex_diff, collapse=', ' )))
    # select columns needed
    df_covar <- df_cov_all %>% select(sample_id, calc_sex) %>% rename(id = sample_id, Sex = calc_sex) %>% unique()
  } 
  
  if (!any(c('Y', 'y') %in% df_unfiltered$chromosome) & !any(c('X', 'x') %in% df_unfiltered$chromosome)) {
    # if x and y both not present, skip gender mapping and map 'unknown's by ratio
    df_cov_all <- df %>% select(sample_id, Sex) %>% unique()
    female_ratio <- nrow(df_cov_all[df_cov_all$Sex == 'female'])/nrow(df_cov_all)
    map_female_num <- floor(nrow(df_cov_all[df_cov_all$Sex == 'unknown']) * female_ratio)
    df_cov_all$Sex[df_cov_all$Sex == 'unknown'][1:map_female_num] <- 'female'
    df_cov_all$Sex[df_cov_all$Sex == 'unknown'] <- 'male' # rest is male
    df_covar <- df_cov_all
  }
  
  return(df_covar)
}


qtl2_foundergeno <- function(df, founder_url, url_list, founders_dict, annot_encode_df){
  founders_df <- df[grep("_", df$original_sample_id), ] %>% 
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
  return(df_founders_encoded)
}

qtl2_ci <- function(summary_df) {
  # join with summary table and CC file in karl's github
  ### if no founder info, filter out.
  # match the strain names with the cc cross info csv
  ci_temp <- read.csv('https://raw.githubusercontent.com/rqtl/qtl2data/main/CC/cc_crossinfo.csv', skip=3) %>% 
    separate(id, c("strain_id", "sec_id"), "/") %>% select(-sec_id)
  ci_sum <- summary_df %>% select(`Strain Name`, `ID`)
  ci_sum$strain_id <- str_extract(ci_sum$`Strain Name`, "CC...")
  ci_sum <- ci_sum %>% select(-`Strain Name`)
  df_crossinfo <- merge(ci_temp, ci_sum, by = 'strain_id', all.y = T) %>% 
    select(ID, everything()) %>%
    rename(id = ID) %>% drop_na() # read_cross does not allow NAs here
  return(df_crossinfo)
}

### get founder strain data from HaploQA if not exist
get_founder_data <- function(founder_url, url_list) {
  # get length to keep track of progress
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
 return(founders_total)
}


# function to generate and convert files to qtl2 input format
get_qtl2_input <- function(data_dir, sample_type, qtl2_output_dir, summary_df, list_pheno) {
  # list container to store outputs
  file_output <- list()
  
  ### get all data from directory
  # read the config file
  config <- fread(paste0(root, '/annotations_config.csv'))
  config <- config %>% filter(array_type == sample_type)
  annot_file <- config$annot_file
  
  ### read annotations
  dir <- "https://raw.githubusercontent.com/kbroman/MUGAarrays/main/UWisc/" # always the same
  annot_df <- read.csv(paste0(dir, annot_file))
  # encode_genome function on geno and founder from qtl2 convert package - args: matrix of genotypes, matrixs of two alleles
  annot_encode_df <- annot_df %>% select(marker, chr, snp) %>% # select allele columns
    separate(snp, c("allele1", "allele2"), sep=cumsum(c(1)))
  
  # all txt file from directory
  data_files <- dir(data_dir, pattern = '\\.txt$', full.names = TRUE)
  ### combine all files
  df_all_init <- rbindlist(lapply(data_files, read_sample_txt)) # save Y chrom here for gender plotting
  
  # clean up summary df
  sum_df <- summary_df %>% 
    select(ID, Sex, `% Het. Calls`, `% No Call`, Platform) %>% rename(sample_id = ID)
  sum_df$`% Het. Calls` <- as.numeric(gsub("%", "", sum_df$`% Het. Calls`))
  sum_df$`% No Call` <- as.numeric(gsub("%", "", sum_df$`% No Call`))
  
  # eliminate those that has no call > 10%
  ### use checkifnot to make sure the markers/SNP are in the annotation file (annotation file is the boss)
  df_all <- merge(df_all_init, sum_df, by = 'sample_id') %>% filter(`% No Call` < 10) %>% filter(!chromosome %in% c('0', 'Y', 'M'))
  
  ### genotype data
  df_geno <- qtl2_geno(df_all, annot_encode_df)
  # save to output
  file_output[[1]] <- df_geno 
  print('genotype data done')
  
  ### gmap and pmap 
  df_maps <- qtl2_maps(df_all, annot_df)
  # save to output
  file_output[[2]] <- df_maps[[1]] # gmap
  file_output[[3]] <- df_maps[[2]] # pmap
  
  ## phenotype data
  # simulate using rnorm/runiform
  df_pheno <- qtl2_pheno(list_pheno, df_all)
  file_output[[4]] <- df_pheno
  print('phenotype data done')
  
  # covariate file - map to sex
  # sample_id, sex
  df_covar <- qtl2_cov(df_all_init, 80, df_all)
  file_output[[5]] <- df_covar
  print('phenotype covariate data done')
  
  ### founder genomes
  # domain of haploqa. Founder data exist on non-dev version
  url_domain <- 'https://haploqa.jax.org'
  # set urls
  founder_url <- 'https://haploqa.jax.org/tag/UNC_Villena_GIGMUGV01_20141012_FinalReport.html'
  # dictionary
  founders_dict <- fread(paste0(root, '/founder_strains_table.csv')) %>% rename(`Strain Name` = founder_strain)
  founder_strains <- unique(founders_dict$`Strain Name`) # get the names of founder strains
  # get html and list of url to sample data
  html_table <- read_html(founder_url) %>% html_nodes("a") %>% html_attr("href")
  url_list <- paste0(url_domain, html_table[grepl('/sample', html_table)])
  
  # get founder data
  fp_founders <- file.path(qtl2_output_dir, 'UNC_Villena_founder_samples.csv')
  if (file.exists(fp_founders)) {
    founders_total <- fread(fp_founders)
  } else { 
    founders_total <- get_founder_data(founder_url)
    write.csv(founders_total, paste0(qtl2_output_dir, '/UNC_Villena_founder_samples.csv'), row.names = F)
  }
  
  df_founders_encoded <- qtl2_foundergeno(founders_total, founder_url, url_list, founders_dict, annot_encode_df)
  # store in output
  file_output[[6]] <- df_founders_encoded 
  print('founder geno data done')
  
  # cross info
  df_crossinfo <- qtl2_ci(summary_df)
  # store in output
  file_output[[7]] <- df_crossinfo 
  print('cross info data done')
  
  return(file_output)
}
