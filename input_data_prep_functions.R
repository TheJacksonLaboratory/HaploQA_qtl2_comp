# contains all functions necessary to retrieve data from HaploQA and convert into qtl2 input format
#
#
library(rstudioapi)
root <- dirname(getSourceEditorContext()$path)
#
# function to retrieve summary table
# @param html_file (html_document) - html of main page sample website on HaploQA, output of read_html(url)
# @param url_list (list) - url to main page of individuals in said sample, bound on summary table as metadata
#
# @return sum_table (data.frame) - table as shown on HaploQA sample website
# columns of sum_table: ID (sample IDs), secondary ID, Haplotype Candidate (T/F), Strain Name, Sex, % Het Calls, % Hom Calls, % No Call, % Concordance, Sample filepath
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

# function to retrieve individuals in a sample
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

# function to generate genotype data for qtl2 input
# @param df (data.frame) - combined dataframe containing all individual files with only necessary columns selected, merged with summary table and chr 0,Y,M and no calls>10% removed
# columns of df: sample_id, original_sample_id, snp_id, chromosome, position_bp, allele1_fwd, allele2_fwd, haplotype1, haplotype2, Sex, % Het. Calls, % No Call, Platform
# @param annot_encode_df (data.frame) - processed UWisc annotation file with alleles
# columns of annot_encode_df: marker, chr, allele1, allele2
#
# @return df_geno_encoded (data.frame) - converted genotype data in qtl2 required format
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

# function to generate genetic and physical map data for qtl2 input
# @param df (data.frame) - combined dataframe containing all individual files with only necessary columns selected, merged with summary table and chr 0,Y,M and no calls>10% removed
# columns of df: sample_id, original_sample_id, snp_id, chromosome, position_bp, allele1_fwd, allele2_fwd, haplotype1, haplotype2, Sex, % Het. Calls, % No Call, Platform
# @param annot_df (data.frame) - raw UWisc annotation file
# columns of annot_df: marker, chr, bp_mm10, cM_cox, cM_g2f1, strand, snp, unique, multi, unmapped, n_blast_hits, n_blast_chr, probe
#
# @return df_map (list) - nested list with two dataframes, one dataframe for gmap and another for pmap
# columns of df_pmap/gmap: marker, chr, pos
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

# function to simulate phenotype data for qtl2 input
# @param list_pheno (list) - list of phenotype names to be simulated
# @param df (data.frame) - combined dataframe containing all individual files with only necessary columns selected, merged with summary table and chr 0,Y,M and no calls>10% removed
# columns of df: sample_id, original_sample_id, snp_id, chromosome, position_bp, allele1_fwd, allele2_fwd, haplotype1, haplotype2, Sex, % Het. Calls, % No Call, Platform
#
#
# @return df_pheno (data.frame) - simulated phenotype data in qtl2 required format
qtl2_pheno <- function(list_pheno, df) {
  df_pheno <- df %>% select(sample_id) %>% unique()
  for (i in seq(1:length(list_pheno))) {
    print('phenotypes:')
    print(list_pheno[i])
    df_pheno[[as.character(list_pheno[i])]] <- floor(runif(nrow(df_pheno), 1, 500))
  }
  
  return(df_pheno)
}

# function to generate phenotype covariate data for qtl2 input
# @param df_raw (data.frame) - raw combined dataframe containing all individual file, nothing removed
# @param split_level (int) - value to split genders
# @param df (data.frame) - combined dataframe containing all individual files with only necessary columns selected, merged with summary table and chr 0,Y,M and no calls>10% removed
# columns of df: sample_id, original_sample_id, snp_id, chromosome, position_bp, allele1_fwd, allele2_fwd, haplotype1, haplotype2, Sex, % Het. Calls, % No Call, Platform
# @param plot_type (string) - type of plot (percent or count)
#
# @return df_covar (list) - phenotype covariate data in qtl2 required format
qtl2_cov <- function(df_raw, split_level, df, plot_type) { # list of unique sample ids instead of df
  split_level <- as.numeric(split_level) # make sure it's a number
  ### if 'x' and 'y' both here, merge and use scatterplot
  # use assert to check if x and y chromosomes both exist
  df_cov_all_xy <- df_raw %>% 
    select(sample_id, snp_id, chromosome, allele1_fwd, allele2_fwd) %>% 
    filter(sample_id %in% unique(df$sample_id)) %>%
    filter(chromosome %in% c('X', 'Y')) %>%
    filter(!snp_id %like% 'PAR')
  if (plot_type %like% 'percent') {
    df_cov_all_xy$alleles <- paste(df_cov_all_xy$allele1_fwd, df_cov_all_xy$allele2_fwd, sep='')
    df_cov_all_xy$alleles_nc <- ifelse(df_cov_all_xy$alleles == '--', 1, 0)
    nc_y <- df_cov_all_xy %>% filter(chromosome == 'Y') %>% group_by(sample_id) %>% summarise(NC_Y = mean(alleles_nc)) %>% as.data.frame()
    df_cov_all_xy$alleles_het <- ifelse(df_cov_all_xy$allele1_fwd != df_cov_all_xy$allele2_fwd, 1, 0)
    het_x <- df_cov_all_xy %>% filter(chromosome == 'X') %>% group_by(sample_id) %>% summarise(Het_X = mean(alleles_het)) %>% as.data.frame()
  }
  if (plot_type %like% 'count') {
    nc_y <- df_cov_all_xy[(df_cov_all_xy$allele1_fwd == '-') & (df_cov_all_xy$allele2_fwd == '-'),] %>% # both are '-', no call
      filter(chromosome == 'Y') %>%
      group_by(sample_id) %>% summarise(NC_Y = n())
    het_x <- df_cov_all_xy[(df_cov_all_xy$allele1_fwd != df_cov_all_xy$allele2_fwd),] %>% # heterozygous
      filter(chromosome == 'X') %>%
      group_by(sample_id) %>% summarise(Het_X = n())
  }
  plot_all <- merge(nc_y, het_x, by = 'sample_id', all = T)
  # scatterplot to show results
  plot <- ggplot(plot_all) + aes(x = Het_X, y = NC_Y, label = sample_id) + geom_point() +
    geom_text_repel(size=3, max.overlaps = 50) + ggtitle(paste0('SNP ', plot_type, ' - Heterozygous on X vs. No call on Y in each sample'))
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
  
  return(df_covar)
}

# function to generate founder genotype data for qtl2 input
# @param df (data.frame) - combined dataframe containing all individual files with only necessary columns selected, merged with summary table and chr 0,Y,M and no calls>10% removed
# columns of df: sample_id, original_sample_id, snp_id, chromosome, position_bp, allele1_fwd, allele2_fwd, haplotype1, haplotype2, Sex, % Het. Calls, % No Call, Platform
# @param founder_url (string) - url to main page of sample that contains founder individuals, for now it's 'https://haploqa.jax.org/tag/UNC_Villena_GIGMUGV01_20141012_FinalReport.html'
# @param url_list (list) - list of urls to main page of founder individuals, retrieved from founder_url
# @param founders_dict (data.frame) - custom built dictionary for founder strain metadata
# columns of founders_dict: founder_strain, letter, color, type, url
# @param annot_encode_df (data.frame) - processed UWisc annotation file with alleles
# columns of annot_encode_df: marker, chr, allele1, allele2
#
#
# @return df_founders_encoded (data.frame) - encoded founder geno data for qtl2 input
# columns of df_founders_encoded: marker, A, B, C, D, E, F, G, H
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

# function to generate cross info data for qtl2 input
# @param summary_df (data.frame) - table as shown on HaploQA sample website, same as output of sample_summary_scrape
# columns of summary_df: ID (sample IDs), secondary ID, Haplotype Candidate (T/F), Strain Name, Sex, % Het Calls, % Hom Calls, % No Call, % Concordance, Sample filepath
#
# @return founders_total (data.frame) - cross info data with sample IDs and mapped strain info for qtl2 input
# columns of founders_total: id, A, B, C, D, E, F, G, H
qtl2_ci <- function(summary_df) {
  # join with summary table and CC file in karl's github
  ### if no founder info, filter out.
  # match the strain names with the cc cross info csv
  ci_sum <- summary_df %>% select(`ID`) %>% unique() %>% rename(id = ID)
  ci_sum$ngen <- 4
  ci_sum[ ,c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')] <- 1
  df_crossinfo <- ci_sum # read_cross does not allow NAs here
  
  return(df_crossinfo)
}

# function to get individual founder strain data from HaploQA, if not exist
# @param founder_url (string) - url to main page of sample that contains founder individuals, for now it's 'https://haploqa.jax.org/tag/UNC_Villena_GIGMUGV01_20141012_FinalReport.html'
# @param url_list (list) - list of urls to main page of founder individuals, retrieved from founder_url
#
# @return founders_total (data.frame) - combined dataframe with all founder individual files
# columns of founders_total: sample_id, original_sample_id, snp_id, chromosome, position_bp, allele1_fwd, allele2_fwd, haplotype1, haplotype2
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

# function to retrieve only the dataframe that all qtl2 input dataframes were based on for testing purposes, and/or the haplotype columns
# @param summary_df (data.frame) - table as shown on HaploQA sample website, same as output of sample_summary_scrape
# @param data_dir (string) - filepath to where individual files (outputs of sample_individual_scrape) are stored
#
# @return df_all (data.frame) - dataframe with haplotype columns
get_haplotypes <- function(summary_df, data_dir) {
  # all txt file from directory
  data_files <- dir(data_dir, pattern = '\\.txt$', full.names = TRUE)
  ### combine all files
  df_raw <- rbindlist(lapply(data_files, read_sample_txt)) # save Y chrom here for gender plotting
  
  # clean up summary df
  sum_df <- summary_df %>% 
    select(ID, Sex, `% Het. Calls`, `% No Call`, Platform) %>% rename(sample_id = ID)
  sum_df$`% Het. Calls` <- as.numeric(gsub("%", "", sum_df$`% Het. Calls`))
  sum_df$`% No Call` <- as.numeric(gsub("%", "", sum_df$`% No Call`))
  
  # eliminate those that has no call > 10%
  ### use checkifnot to make sure the markers/SNP are in the annotation file (annotation file is the boss)
  df_all <- merge(df_raw, sum_df, by = 'sample_id') %>% filter(`% No Call` < 10) %>% filter(!chromosome %in% c('0', 'Y', 'M'))
  
  return(df_all)
}


# function to call the above functions, generate and convert files to qtl2 input format 
# @param data_dir (string) - filepath to where individual files (outputs of sample_individual_scrape) are stored
# @param sample_type (string) - type of sample being analyzed - GigaMUGA, MiniMUGA, Collaborative Cross, etc.
# @param qtl2_output_dir (string) - directory to store qtl2 output files
# @param summary_df (data.frame) - table as shown on HaploQA sample website, same as output of sample_summary_scrape
# columns of summary_df: ID (sample IDs), secondary ID, Haplotype Candidate (T/F), Strain Name, Sex, % Het Calls, % Hom Calls, % No Call, % Concordance, Sample filepath
# @param list_pheno (list) -  - list of phenotype names to be simulated
#
# @return file_output (list) - nested list containing all files necessary for qtl2 input
# specifically: df_geno, df_gmap, df_pmap, df_pheno, df_covar, df_foundergeno, df_crossinfo
# see individual functions above for column names of each dataframe
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
  df_raw <- rbindlist(lapply(data_files, read_sample_txt)) # save Y chrom here for gender plotting
  
  # clean up summary df
  sum_df <- summary_df %>% 
    select(ID, Sex, `% Het. Calls`, `% No Call`, Platform) %>% rename(sample_id = ID)
  sum_df$`% Het. Calls` <- as.numeric(gsub("%", "", sum_df$`% Het. Calls`))
  sum_df$`% No Call` <- as.numeric(gsub("%", "", sum_df$`% No Call`))
  
  # eliminate those that has no call > 10%
  ### use checkifnot to make sure the markers/SNP are in the annotation file (annotation file is the boss)
  df_all <- merge(df_raw, sum_df, by = 'sample_id') %>% filter(`% No Call` < 10) %>% filter(!chromosome %in% c('0', 'Y', 'M'))
  stopifnot("There are markers in qtl2 dataframe not present in annotation file, check validity or exclude such markers" = any(!unique(df_all$snp_id) %in% annot_df$marker) == FALSE) # if there is any marker from df that's not in annotation, raise error
  
  ### genotype data
  df_geno <- qtl2_geno(df_all, annot_encode_df)
  stopifnot("NAs present in genotype data" = any(is.na(df_geno)) == FALSE) # if there is any NAs, raise error
  # save to output
  file_output[[1]] <- df_geno 
  print('genotype data done')
  
  ### gmap and pmap 
  df_maps <- qtl2_maps(df_all, annot_df)
  stopifnot("NAs present in genetic map data" = any(is.na(df_maps[[1]])) == FALSE) # if there is any NAs, raise error
  stopifnot("NAs present in physical map data" = any(is.na(df_maps[[2]])) == FALSE) # if there is any NAs, raise error
  # save to output
  file_output[[2]] <- df_maps[[1]] # gmap
  file_output[[3]] <- df_maps[[2]] # pmap
  
  ## phenotype data
  # simulate using rnorm/runiform
  df_pheno <- qtl2_pheno(list_pheno, df_all)
  stopifnot("NAs present in phenotype data" = any(is.na(df_pheno)) == FALSE) # if there is any NAs, raise error
  file_output[[4]] <- df_pheno
  print('phenotype data done')
  
  # covariate file - map to sex
  # sample_id, sex
  df_covar <- qtl2_cov(df_raw, 80, df_all, 'percent')
  stopifnot("NAs present in phenotype covariate data" = any(is.na(df_covar)) == FALSE)
  stopifnot("Sample IDs in covariate data do not align with input" = length(setdiff(unique(df_covar$id), unique(df_all$sample_id))) == 0)
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
  fp_founders <- file.path(root, 'UNC_Villena_founder_samples.csv')
  if (file.exists(fp_founders)) {
    founders_total <- fread(fp_founders)
  } else { 
    founders_total <- get_founder_data(founder_url, url_list)
    write.csv(founders_total, paste0(root, '/UNC_Villena_founder_samples.csv'), row.names = F)
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

# function to sort chromosomes into custom order
# @param df (data.frame) - selected dataframe to be sorted
# @param sort_order (list) - order to sort the chromosomes
#
# @return df (data.frame) - dataframe with chromosomes sorted in the selected order
sort_chr <- function(df, sort_order) {
  df$chr <- factor(df$chr, levels = sort_order)
  df <- df[order(df$chr),] 
  return(df)
}

## if rds file don't exist for the sample type, read cross
save_rds <- function(rds_dir, qtl2_dir, sample_type, results_dir) {
  # check if control file exist
  control_fp <- paste0(qtl2_dir, '/test.json')
  print(control_fp)
  if (file.exists(control_fp) == FALSE) {
    stop(paste0('Cross file not found in ', qtl2_dir))
  }

  cross_list <- list()
  # read cross file
  cross_file <- read_cross2(control_fp)
  cross_list[['cross']] <- cross_file
  
  ### save map individually as well
  map <- cross_file$gmap
  map_fp <- paste0(rds_dir, "/map_", sample_type, ".rds")
  saveRDS(map, file = map_fp)
  cross_list[['map']] <- map
  
  # calculations
  ### save rds files: pr, ginf, ph_geno, pos
  # genoprob
  pr_fp <- paste0(rds_dir, "/pr_", sample_type, ".rds")
  print(pr_fp)
  if (!file.exists(pr_fp)) {
    print('genoprob file does not exist, running calculations')
    pr <- calc_genoprob(cross=cross_file, map=map, error_prob=0.002, cores = 2)
    print('saving genoprob file')
    cross_list[['pr']] <- pr
    saveRDS(pr, file = pr_fp)
  } else {
    print('genoprob file exists, reading in file')
    pr <- readRDS(pr_fp)
    cross_list[['pr']] <- pr
  }
  
  # maxmarg
  ginf_fp <- paste0(rds_dir, "/ginf_", sample_type, ".rds")
  if(!file.exists(ginf_fp)) {
    print('maxmarg file does not exist, running calculations')
    ginf <- maxmarg(pr, minprob = 0.01) # lower minprob to 0.01
    cross_list[['ginf']] <- ginf
    print('saving maxmarg file')
    saveRDS(ginf, file = ginf_fp)
  } else {
    print('maxmarg file exists, reading in file')
    ginf <- readRDS(ginf_fp)
    cross_list[['ginf']] <- ginf
  }
  
  # phased geno
  ph_geno_fp <- paste0(rds_dir, "/ph_geno_", sample_type, ".rds")
  if(!file.exists(ph_geno_fp)) {
    print('phased geno file does not exist, running calculations')
    ph_geno <- guess_phase(cross_file, ginf)
    print('saving phased file')
    cross_list[['ph_geno']] <- ph_geno
    saveRDS(ph_geno, file = ph_geno_fp)
  } else {
    print('phased geno file exists, reading in file')
    ph_geno <- readRDS(ph_geno_fp)
    cross_list[['ph_geno']] <- ph_geno
  }
  
  # crossover locations
  pos_fp <- paste0(rds_dir, "/pos_", sample_type, ".rds")
  if (!file.exists(pos_fp)) {
    print('location crossover file does not exist, running calculations')
    pos <- locate_xo(ginf, map) 
    cross_list[['pos']] <- pos
    print('saving location crossover file')
    saveRDS(pos, file = pos_fp)
  } else {
    print('location crossover file exists, reading in file')
    pos <- readRDS(pos_fp)
    cross_list[['pos']] <- pos
  }
  

  return(cross_list)
}

rds_loop <- function(samples_for_rds, config, results_dir, rds_dir, num_chr) {
  list_result_all <- list()
  for (sample in samples_for_rds) {
    print(sample)
    sample_config <- config[array_type == sample]
    
    # according directories
    data_dir <- file.path(root, sample_config$data_dir)
    qtl2_dir <- file.path(root, sample_config$qtl2_dir)
    cross_list <- save_rds(rds_dir, qtl2_dir, sample, results_dir)
    
    summary_df <- fread(paste0(data_dir, '/', sample, '_summary.csv'))
    cross_list[['summary']] <- summary_df
    cross_list[['data_dir']] <- data_dir
    
    founder_all_codes <- colnames(cross_list[['pr']]$`X`)
    
    ginf_haploqa <- maxmarg_sim(summary_df, data_dir, num_chr, founder_all_codes, cross = cross_list[['cross']])
    cross_list[['ginf_haploqa']] <- ginf_haploqa
    
    ph_geno_haploqa <- guess_phase_sim(ginf_haploqa, cross_list[['cross']], rds_dir, sample)
    cross_list[['ph_geno_haploqa']] <- ph_geno_haploqa
    
    list_result_all[[sample]] <- cross_list
  } 
  
  
  return(list_result_all)
}


make_cross <- function(crosstype, chromosomes, df, x_only) {
  attr_chr <- setNames(num_chr == 'X', chromosomes)
  if(x_only == T) {
    attr(df, "is_x_chr") <- attr_chr
  } else {
    attr(df, "crosstype") <- crosstype
    attr(df, "alleles") <- c("A", "B", "C", "D", "E", "F", "G", "H")
    attr_chr <- setNames(num_chr == 'X', chromosomes)
    attr(df, "is_x_chr") <- attr_chr
    attr(df, "class") <- c('viterbi', 'list')
  }
  return(df)
}


shiny_viz_input <- function(sample_type, rds_dir) {
  phased_geno_comp_list <- list()
  map <- readRDS(paste0(rds_dir, '/map_', sample_type, '.rds')) # should share the same map
  ph_geno_haplo <- readRDS(paste0(rds_dir, '/ph_geno_haploqa_', sample_type, '.rds'))
  ph_geno_qtl2 <- readRDS(paste0(rds_dir, '/ph_geno_', sample_type, '.rds'))
  phased_geno_comp_list[['haplo']] <- ph_geno_haplo
  phased_geno_comp_list[['qtl2']] <- ph_geno_qtl2
  phased_geno_comp_list[['map']] <- map
  
  return(phased_geno_comp_list)
}


maxmarg_sim <- function(summary_df, data_dir, num_chr, founder_all_codes, cross) {
  ## look-up tables
  founders_dict <- fread(paste0(root, '/founder_strains_table.csv'))
  founder_haplo_lookup <- setNames(founders_dict$letter, founders_dict$founder_strain)

  all_map_codes <- seq(1:length(founder_all_codes))
  founder_all_lookup <- setNames(all_map_codes, founder_all_codes) # make maxmarg
  founder_all_rev_lookup <- setNames(founder_all_codes, all_map_codes)
  
  
  df_all <- get_haplotypes(summary_df, data_dir)
  df_test <- df_all %>% select(sample_id, snp_id, haplotype1, haplotype2, chromosome) %>% unique()
  df_test[,c(3,4)] <- as.data.frame(apply(df_test[,c(3,4)], 2, function(x) founder_haplo_lookup[x]))
  df_test$haplotype <- paste(df_test$haplotype2, df_test$haplotype1, sep='')
  
  haploqa_maxmarg <- df_test %>% select(sample_id, snp_id, chromosome, haplotype)
  
  haploqa_maxmarg[,4] <- lapply(haploqa_maxmarg[,4], function(col) {
    rev_col = stri_reverse(col)
    ifelse(rev_col %in% founder_all_codes, rev_col, col)
  })
  
  haploqa_maxmarg[,4] <- as.data.frame(apply(haploqa_maxmarg[,c(4)], 2, function(x) founder_all_lookup[x]))
  
  ## save per chromosome
  ginf_haploqa <- list()
  for (i in num_chr) {
    #i <- 'X'
    print(i)
    df <- haploqa_maxmarg %>% filter(chromosome == i) %>% dcast(sample_id ~ snp_id, value.var = 'haplotype') %>% as.data.frame()
    rownames(df) <- df$sample_id
    
    df <- df %>% select(-c(sample_id)) %>% as.matrix()
    
    col_order <- names(cross$gmap[[i]])
    df <- df[,col_order]
    ginf_haploqa[[i]] <- df
  }
  
  ginf_haploqa <- make_cross('genail8', num_chr, ginf_haploqa, x_only = F)
  
  return(ginf_haploqa)
}

guess_phase_sim <- function(ginf_haploqa, cross, rds_dir, sample_type) {
  ph_geno_haploqa <- guess_phase(cross, ginf_haploqa)
  # save this separately for easy shiny input
  saveRDS(ph_geno_haploqa, paste0(rds_dir, '/ph_geno_haploqa_', sample_type,'.rds'))
  return(ph_geno_haploqa)
}


location_xo_comp <- function(qtl2_pos, haploqa_pos, num_chr) {
  x_loc_diff <- list()
  for (chr in num_chr) {
    print(chr)
    pos_h_chr <- pos_haploqa_cc[[chr]]
    pos_chr <- pos_qtl2_cc[[chr]]
    nmar_h <- sapply(pos_h_chr, length)
    nmar <- sapply(pos_chr, length) # none has the same length...
    #x_loc_diff[[chr]] <- list()
    for (i in seq(1:length(pos_chr))) {
      x <- unlist(pos_h_chr[i])
      y <- unlist(pos_chr[i])
      if (length(x) > length(y)) {
        df <- apply(abs(outer(x, y, '-')), 2, min)
      }
      if (length(x) < length(y)) {
        df <- apply(abs(outer(y, x, '-')), 2, min)
      }
      x_loc_diff[[chr]][[i]]<- df
    }
  }
  return(x_loc_diff)
}


err_comp <- function(pr, map, num_chr) {
  ### have this take sample ID instead
  for (chr in num_chr) {
    #chr <- 1
    print(chr)
    err <- (1 - apply(pr[[chr]], c(1,3), max)) %>% t() %>% as.data.frame()
    err$marker <- rownames(err)
    err_map <- map[[chr]] %>% as.data.frame() %>% rename(pos = '.')
    err_map$marker <- rownames(err_map)
    
    err_chr <- merge(err, err_map, by = 'marker')
    df <- err_chr[,c(2, ncol(err_chr))]
    # facet wrap each chromosome
    plot <- ggplot(df) + aes(x = df[,2], y = df[,1]) + geom_line() + xlab('pos') + ylab(names(df)[1])
    
    print(plot)
  }
  return(df)
}

comp_df_int <- function(ginf, lookup_table, map_chr) {
  df <- as.data.frame(ginf[[chr]])
  df_chr <- apply(df, 2, function(x) founder_lookup_table[x]) # column level
  rownames(df_chr) <- rownames(df)
  df_chr <- as.data.frame(t(df_chr))
  df_chr$marker <- rownames(df_chr)
  df_chr_comp <- merge(df_chr, map_chr, by = 'marker')
  
  df_comp <- df_chr_comp[1:(ncol(df_chr_comp))] %>% select(-c(marker))
  
  return(df_comp)
}

genocode_comp_matrix <- function(map, ginf_qtl2, ginf_haploqa, founder_lookup_table, num_chr, file_gen) {
  comp_matrix <- list()
  for (chr in num_chr) {
    ## qtl2
    chr <- '1'
    print(chr)
    map_chr <- as.data.frame(map[[chr]]) %>% rename('pos' = 'map[[chr]]')
    map_chr$marker <- rownames(map_chr)

    haploqa_comp <- comp_df_int(ginf_haploqa, founder_lookup_table, map_chr)
    qtl2_comp <- comp_df_int(ginf_qtl2, founder_lookup_table, map_chr)
    
    ### loop through each matching column within both dataframes
    result_df <- data.frame()
    for (col in seq(1:(ncol(qtl2_comp)-1))) {
      #print(col)
      df_match <- data.frame('qtl2_ind' = qtl2_comp[,col], 'haploqa_ind' = haploqa_comp[,col])
      t1 <- dcast(df_match, qtl2_ind ~ haploqa_ind, value.var = 'haploqa_ind', fun.aggregate = length)
      result_df <- bind_rows(result_df, t1) %>% group_by(qtl2_ind) %>% summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>% as.data.frame()
    }
    
    result_df <- result_df %>% rename(geno_code = qtl2_ind) %>%
      select(geno_code, AA, BB, CC, DD, EE, FF, GG, HH, sort(names(.)))
    
    comp_matrix[[chr]] <- result_df
    
    if(file_gen == T) {
      comp_fp <- file.path(root, 'geno_comp_results')
      dir.create(comp_fp, showWarnings = FALSE) 
      write.csv(result_df, paste0(comp_fp, '/geno_comp_matrix_', chr, '.csv'), row.names = F)
    }
  }
  return(comp_matrix)
}


get_geno_pct_diff <- function(ginf_haploqa, ginf_qtl2, summary_df, num_chr) {
  haploqa_comp <- comp_df_int(ginf_haploqa, founder_lookup_table, map_chr)
  qtl2_comp <- comp_df_int(ginf_qtl2, founder_lookup_table, map_chr)
  
  for (chr in num_chr) {
    chr_pct <- c()
    print(chr)
    for (ind in seq(1:(ncol(qtl2_comp)-1))) {
      #ind <- 1
      #ind_index <- ind+1
      ind_qtl2 <- qtl2_comp[,ind]
      ind_haplo <- haploqa_comp[,ind]
      diff <- data.frame('qtl2' = ind_qtl2, 'haploqa' = ind_haplo, 'pos' = qtl2_comp$pos)
      diff <- diff %>% arrange(pos) %>% filter(qtl2 != haploqa)
      #write.csv(diff, paste0(comp_dir, '/diff_ind_', ind, '_chr_', chr, '.csv'), row.names = F)
      
      diff_pct <- +(!((haploqa_comp[1:(ncol(haploqa_comp)-1)] == qtl2_comp[1:(ncol(qtl2_comp)-1)]) * 1))
      ind_df <- as.data.frame(diff_pct[,ind])
      pct <- as.numeric(colSums(ind_df) / nrow(ind_df))
      chr_pct <- c(chr_pct, pct)
    }
    
    chr_diff <- data.frame('individuals' = colnames(qtl2_comp[1:(ncol(qtl2_comp)-1)]), 'genome_pct_diff' = chr_pct)
    chr_diff$chr <- chr
    chr_diff$ind <- seq(1:(ncol(qtl2_comp)-1))
    
    het_df <- summary_df %>% select(ID, `% Het. Calls`) %>% rename(individuals = ID, het_pct = '% Het. Calls')
    plot_diff <- merge(chr_diff, het_df, by = 'individuals')
    plot_diff$het_pct <- as.numeric(gsub("%", "", plot_diff$het_pct)) / 100
    # histogram 
    ggplot(plot_diff) + aes(x = genome_pct_diff, y = het_pct) + geom_point()
    ggplot(plot_diff, aes(x=genome_pct_diff)) + geom_histogram(binwidth = .01)
    
    #write.csv(chr_diff, paste0(results_dir, '/pct_diff_chr_', chr, '.csv'), row.names = F)
  }
  return(chr_pct)
}


