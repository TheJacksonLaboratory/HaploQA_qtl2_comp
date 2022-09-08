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
# rows
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
  df_gmap <- sort_chr(df_gmap, c((1:19),"X"))
  df_gmap <- df_gmap %>% group_by(chr) %>% arrange(pos, .by_group = TRUE) %>% as.data.table()
  df_maps[[1]] <- df_gmap
  marker_order <- unique(df_gmap$marker)
  ## physical map data (pos in Mbp)
  df_pmap <- map_df %>% select(marker, chr.x, bp_mm10) %>% rename(chr = chr.x, pos = bp_mm10) %>%
    mutate(pos = pos/1000000) %>% drop_na()
  df_pmap <- sort_chr(df_pmap, c((1:19),"X"))
  df_pmap <- df_pmap %>% group_by(chr) %>% arrange(pos, .by_group = TRUE) %>% as.data.table()
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
  if (split_level > 1) {
    warning('Split level value entered is largrer than 1. Dividing by 100 to get decimal percentage value')
    split_level <- as.numeric(split_level / 100)
    print(paste0('splitting genders at ', split_level))
  } else {
    print(paste0('splitting genders at ', split_level))
  }
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
  females_list <- plot_all[plot_all$NC_Y > split_level,]$sample_id # high no call on Y - females
  males_list <- plot_all[plot_all$NC_Y < split_level,]$sample_id
  # create new column
  df_cov_all <- df %>% select(sample_id, Sex) %>% unique()
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
qtl2_foundergeno <- function(df, founder_url, url_list, founders_dict, annot_encode_df, founders_list = NULL, marker_order, founder_haplo_lookup){
  founder_codes <- unique(founders_dict$`Strain Name`)
  sum_df <- sample_summary_scrape(read_html(founder_url), url_list) %>% 
    select('ID', 'Strain Name', 'Secondary IDs') %>% rename(sample_id = ID, original_sample_id = 'Secondary IDs')
  if((all(founders_list %in% founder_codes)) & (length(founders_list) == 8)) { ## if 8 founder strains, use founder dictionary to map to according letters
    founders_df <- df[grep("_", df$original_sample_id), ] %>% 
      mutate(gene_exp = paste(allele1_fwd, allele2_fwd, sep = '')) %>%
      select('sample_id', 'snp_id', 'gene_exp')
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
    df_founders_encoded <- df_founders_encoded %>% filter(marker %in% marker_order) %>% arrange(factor(marker, levels = marker_order))
    df_founders_encoded <- df_founders_encoded %>% select(marker, everything())
  } else {
    for (founder in founders_list) {  # if not 8 founder strains, custom select letters
      founders_df <- merge(df, sum_df, by = 'original_sample_id') %>% 
        mutate(gene_exp = paste(allele1_fwd, allele2_fwd, sep = '')) %>%
        rename(sample_id = 'sample_id.y') %>%
        select('sample_id', 'snp_id', 'gene_exp', 'Strain Name', 'original_sample_id')
      
      founders_test <- founders_df[founders_df$`Strain Name` %in% founders_list,]
      haplo_lookup <- data.frame(strain=names(founder_haplo_lookup), letter=founder_haplo_lookup) %>% rename('Strain Name'=strain)
      founders_test <- merge(founders_test, haplo_lookup, by = 'Strain Name')
      
      test_foundergeno <- founders_test %>% select(snp_id, letter, gene_exp) %>% 
        rename(strain = letter, marker = snp_id)
      test_foundergeno <- dcast(test_foundergeno, marker~strain, value.var="gene_exp")
      rownames(test_foundergeno) <- test_foundergeno$marker
      
      # put the allele annotations in according order
      founder_encode_annot <- merge(annot_encode_df, test_foundergeno, by = 'marker') %>% select(allele1, allele2)
      df_founders_encoded <- as.data.frame(encode_geno(test_foundergeno[,-1], founder_encode_annot))
      df_founders_encoded$marker <- rownames(df_founders_encoded)
      df_founders_encoded <- df_founders_encoded %>% filter(marker %in% marker_order) %>% arrange(factor(marker, levels = marker_order))
      df_founders_encoded <- df_founders_encoded %>% select(marker, everything())
      
  } 
  }
  
  
  return(df_founders_encoded)
}

# function to generate cross info data for qtl2 input
# @param summary_df (data.frame) - table as shown on HaploQA sample website, same as output of sample_summary_scrape
# columns of summary_df: ID (sample IDs), secondary ID, Haplotype Candidate (T/F), Strain Name, Sex, % Het Calls, % Hom Calls, % No Call, % Concordance, Sample filepath
#
# @return founders_total (data.frame) - cross info data with sample IDs and mapped strain info for qtl2 input
# columns of founders_total: id, A, B, C, D, E, F, G, H
qtl2_ci <- function(summary_df, ngen, founders_list) {
  founder_list <- LETTERS[seq(1,8)]
  # join with summary table and CC file in karl's github
  ### if no founder info, filter out.
  # match the strain names with the cc cross info csv
  ci_sum <- summary_df %>% select(`ID`) %>% unique() %>% rename(id = ID)
  ci_sum$ngen <- ngen
  ci_sum[ ,c(founder_list[1:length(founders_list)])] <- 1
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
    df_temp <- as.data.frame(content(GET(file), encoding="UTF-8"))
    founders_total <- rbind(founders_total,df_temp)
  }
  
  df <- read.tsv(textConnection(content(GET(file), 'text')))
  
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
get_qtl2_input <- function(data_dir, sample_type, annot_file, qtl2_output_dir, summary_df, list_pheno, ngen, founders_list, marker_type, exclude_list) {
  # list container to store outputs
  file_output <- list()
  
  ### read annotations
  dir <- "https://raw.githubusercontent.com/kbroman/MUGAarrays/main/UWisc/" # always the same
  annot_df <- read.csv(paste0(dir, annot_file))
  # encode_genome function on geno and founder from qtl2 convert package - args: matrix of genotypes, matrixs of two alleles
  annot_encode_df <- annot_df %>% select(marker, chr, snp) %>% # select allele columns
    separate(snp, c("allele1", "allele2"), sep=cumsum(c(1)))
  
  # all txt file from directory
  data_files <- dir(data_dir, pattern = '\\.txt$', full.names = TRUE)
  ### combine all files
  df_raw <- rbindlist(lapply(data_files, read_sample_txt)) %>% filter(!sample_id %in% exclude_list) # save Y chrom here for gender plotting
  
  # clean up summary df
  sum_df <- summary_df %>% 
    select(ID, Sex, `% Het. Calls`, `% No Call`, Platform) %>% rename(sample_id = ID)
  sum_df$`% Het. Calls` <- as.numeric(gsub("%", "", sum_df$`% Het. Calls`))
  sum_df$`% No Call` <- as.numeric(gsub("%", "", sum_df$`% No Call`))
  
  # eliminate those that has no call > 10%
  ### use checkifnot to make sure the markers/SNP are in the annotation file (annotation file is the boss)
  df_all <- merge(df_raw, sum_df, by = 'sample_id') %>% filter(`% No Call` < 10) %>% filter(!chromosome %in% c('0', 'Y', 'M')) %>% filter(!sample_id %in% exclude_list) 
  stopifnot("There are markers in qtl2 dataframe not present in annotation file, check validity or exclude such markers" = any(!unique(df_all$snp_id) %in% annot_df$marker) == FALSE) # if there is any marker from df that's not in annotation, raise error
  
  
  ## make dictionary
  if (sample_type %in% c('CC', 'DO')) { ### change this to sample_type in CC or DO, they are special conditions
    founders_dict <- fread(paste0(root, '/founder_strains_table.csv'))
    founder_haplo_lookup <- setNames(founders_dict$letter, founders_dict$founder_strain)
    file_output[[8]] <- founder_haplo_lookup
  } else {
    unique_haplotypes <- unique(df_all[,c(haplotype1, haplotype2)])
    founder_haplo_lookup <- setNames(LETTERS[seq(1, length(unique(unique_haplotypes)))], unique(unique_haplotypes))
    file_output[[8]] <- founder_haplo_lookup
  }
  
  
  ### genotype data
  ### to-do: sort the genotype markers into same order of gmap and pmap, and make sure all markers are consistent
  # ensure founders are coded consistently
  df_geno <- qtl2_geno(df_all, annot_encode_df)
  stopifnot("NAs present in genotype data" = any(is.na(df_geno)) == FALSE) # if there is any NAs, raise error
  print('genotype data done')
  
  ### gmap and pmap 
  df_maps <- qtl2_maps(df_all, annot_df)
  stopifnot("NAs present in genetic map data" = any(is.na(df_maps[[1]])) == FALSE) # if there is any NAs, raise error
  stopifnot("NAs present in physical map data" = any(is.na(df_maps[[2]])) == FALSE) # if there is any NAs, raise error
  # save to output
  file_output[[2]] <- df_maps[[1]] # gmap
  file_output[[3]] <- df_maps[[2]] # pmap
  
  # reorder genotype data
  marker_order <- unique(df_maps[[1]]$marker)
  df_geno <- df_geno %>% filter(marker %in% marker_order) %>% arrange(factor(marker, levels = marker_order))
  # save to output
  file_output[[1]] <- df_geno 
  
  ## phenotype data
  # simulate using rnorm/runiform
  df_pheno <- qtl2_pheno(list_pheno, df_all)
  stopifnot("NAs present in phenotype data" = any(is.na(df_pheno)) == FALSE) # if there is any NAs, raise error
  file_output[[4]] <- df_pheno
  print('phenotype data done')
  
  # covariate file - map to sex
  # sample_id, sex
  df_covar <- qtl2_cov(df_raw, 0.8, df_all, 'percent')
  stopifnot("NAs present in phenotype covariate data" = any(is.na(df_covar)) == FALSE)
  stopifnot("Sample IDs in covariate data do not align with input" = length(setdiff(unique(df_covar$id), unique(df_all$sample_id))) == 0)
  file_output[[5]] <- df_covar
  print('phenotype covariate data done')
  
  ### founder genomes
  # domain of haploqa. Founder data exist on non-dev version
  url_domain <- 'http://haploqa.jax.org/'
  # set urls
  founder_url <- 'http://haploqa.jax.org//tag/UNC_Villena_GIGMUGV01_20141012_FinalReport.html'
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
  
  df_founders_encoded <- qtl2_foundergeno(founders_total, founder_url, url_list, founders_dict, annot_encode_df, founders_list, marker_order, founder_haplo_lookup)
  # store in output
  file_output[[6]] <- df_founders_encoded 
  print('founder geno data done')
  
  # cross info
  df_crossinfo <- qtl2_ci(summary_df, ngen, founders_list)
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

# function to store results into rds files for future use
# @param rds_dir (string) - directory for which the rds files should be saved in
# @param qtl2_dir (string) - directory to read the qtl2 input files from
# @param results_dir (string) - directory to read qtl2 computation results from
# @param data_dir (string) - directory to read raw haploqa files from
# @param n_founders (int) - number of founders that exists
# @param founder_haplo_lookup (named list) lookup table to convert founder strains to according codes
#
# @return cross_list (list containing multiple objects) - a list that has all the objects necessary for further research
# such as genetic map, calculated qtl2 objects, phased genotypes for both qtl2 and haploqa, etc.
save_rds <- function(rds_dir, qtl2_dir, sample_type, results_dir, data_dir, n_founders, founder_haplo_lookup) {
  num_chr <- c((1:19),"X")
  print(sample_type)
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
  
  # get summary file
  summary_df <- fread(paste0(data_dir, '/', sample_type, '_summary.csv'))
  cross_list[['summary']] <- summary_df
  
  cross_list[['data_dir']] <- data_dir
  
  founder_lookup_table <- fread(file.path(root, 'founder_lookup_table.csv'))
  founder_all_rev_lookup <- setNames(founder_lookup_table$founder_codes, founder_lookup_table$founder_id)
  
  ginf_haploqa <- maxmarg_sim(summary_df, data_dir, num_chr, founder_all_rev_lookup, cross = cross_list[['cross']], qtl2_dir, sample_type, n_founders, founder_haplo_lookup)
  cross_list[['ginf_haploqa']] <- ginf_haploqa
  
  ph_geno_haploqa <- guess_phase_sim(ginf_haploqa, cross_list[['cross']], rds_dir, sample_type)
  cross_list[['ph_geno_haploqa']] <- ph_geno_haploqa

  return(cross_list)
}


# function to add cross attributes to an object
# @param crosstype (string) - type of cross model
# @param num_chr (string) - chromosome to be mapped
# @param chromosomes (string) - chromosomes present in the model
# @param df (data.frame/list) - object to add cross attributes to
# @param x_only (TRUE/FALSE) - only change the is_x_chr attribute
#
# @return df (data.frame/list) - object with attributes added
make_cross <- function(crosstype, num_chr, chromosomes, df, x_only) {
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

# function to pull only the inputs needed for shiny app
# @param sample_type (string) - type of haploqa sample
# @param rds_dir (string) - directory to read rds files from
#
# @return phased_geno_comp_list (list of dataframes) - list containing all objects needed for shiny
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

# function to simulate the maxmarg phasing computations
# @param summary_df (data.frame) - summary table scraped from haploqa
# @param data_dir (string) - directory where raw input data is stored
# @param num_chr (vector) - chromosomes present, sorted into order
# @param founder_all_rev_lookup (named list) - lookup table to convert number codes (1,2,3,4) into letters(AA, BB, CC, DD)
# @param cross (cross object) - output of read_cross for the sample model
# @param qtl2_dir (string) - directory where qtl2 input data is stored
# @param sample_type (string) - type of sample, such as CC, DO, F2, etc.
# @param n_founders (int) - number of founders exist in the model
# @param founder_haplo_lookup (named list) - lookup table to convert strain names (A/J, C57BL/6J) into codes (A, B)
#
# @return ginf_haploqa (list of dataframes) - phased genotypes indexed by chromosome, same as output of maxmarg() from qtl2
maxmarg_sim <- function(summary_df, data_dir, num_chr, founder_all_rev_lookup, cross, qtl2_dir, sample_type, n_founders, founder_haplo_lookup) {
  
  df_cov <- fread(file.path(qtl2_dir, 'test_covar.csv')) %>% rename('sample_id' = 'id')
  
  df_all <- get_haplotypes(summary_df, data_dir)
  df_test <- df_all %>% select(sample_id, snp_id, haplotype1, haplotype2, chromosome) %>% unique()
  df_test[,c(3,4)] <- as.data.frame(apply(df_test[,c(3,4)], 2, function(x) founder_haplo_lookup[x]))
     
  # gender mapping
  df_test <- merge(df_test, df_cov, by = 'sample_id')
  df_test[(df_test$chromosome == 'X') & (df_test$Sex == 'male')]$haplotype2 <- 'Y'
  
  df_test$haplotype <- paste(df_test$haplotype2, df_test$haplotype1, sep='')
  
  haploqa_maxmarg <- df_test %>% select(sample_id, snp_id, chromosome, haplotype)
  
  haploqa_maxmarg[,4] <- lapply(haploqa_maxmarg[,4], function(col) {
    rev_col = stri_reverse(col)
    ifelse(rev_col %in% founder_all_rev_lookup, rev_col, col)
  })
  
  ### join sex into this df - for males, take first letter, then paste Y after it
  founder_all_lookup <- setNames(names(founder_all_rev_lookup), founder_all_rev_lookup) # reverse the lookup table
  founder_all_lookup <- founder_all_lookup[names(founder_all_lookup) %in% unique(haploqa_maxmarg$haplotype)]
  founder_all_lookup <- setNames(seq(1, length(founder_all_lookup)), names(founder_all_lookup))
  haploqa_maxmarg[,4] <- as.data.frame(apply(haploqa_maxmarg[,c(4)], 2, function(x) founder_all_lookup[x]))
  ### subset the lookup table names to just what is in haploqa_maxmarg$haplotype
  ### then assign numbers seq(1, length())
  
  
  ## save per chromosome
  ginf_haploqa <- list()
  for (i in num_chr) {
    #i <- num_chr[1]
    print(i)
    df <- haploqa_maxmarg %>% filter(chromosome == i) %>% dcast(sample_id ~ snp_id, value.var = 'haplotype') %>% as.data.frame()
    df[,2:ncol(df)] <- mutate_all(df[,2:ncol(df)], function(x) as.numeric(as.character(x)))
    rownames(df) <- df$sample_id
    
    df <- df %>% select(-c(sample_id)) %>% as.matrix()
    
    col_order <- names(cross$gmap[[i]])
    df <- df[,col_order]
    ginf_haploqa[[i]] <- df
  }
  
  attr(ginf_haploqa, "crosstype") <- paste0('genail', n_founders)
  
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
    #chr <- '18'
    print(chr)
    pos_h_chr <- haploqa_pos[[chr]]
    pos_chr <- qtl2_pos[[chr]]
    nmar_h <- sapply(pos_h_chr, length)
    nmar <- sapply(pos_chr, length) # none has the same length...
    #x_loc_diff[[chr]] <- list()
    for (i in seq(1:length(pos_chr))) {
      #i <- 149
      x <- unlist(pos_h_chr[i])
      y <- unlist(pos_chr[i])
      if (length(x) > length(y)) {
        df <- apply(abs(outer(x, y, '-')), 2, min)
      }
      if (length(x) < length(y)) {
        df <- apply(abs(outer(y, x, '-')), 2, min)
      }
      x_loc_diff[[chr]][[i]] <- ifelse(is.na(mean(df)), 0, mean(df))
    } # make a histogram of some sort, 100 bins
  }
  return(x_loc_diff)
}


err_comp <- function(pr, map, num_chr, results_dir, sample_type) {
  pdf(paste0(results_dir, 'err_comp_plots_', sample_type, '.pdf'))
  facet_all <- data.frame()
  for (chr in num_chr) {
    #chr <- 2
    print(chr)
    err <- (1 - apply(pr[[chr]], c(1,3), max)) %>% t() %>% as.data.frame()
    err$marker <- rownames(err)
    err_map <- map[[chr]] %>% as.data.frame() %>% rename(pos = '.')
    err_map$marker <- rownames(err_map)
    
    err_chr <- merge(err, err_map, by = 'marker')
    
    for (col in seq(2,(ncol(err_chr)-2))) {
      print(col)
      #col <- 2
      df <- err_chr[,c(col, ncol(err_chr))]
      df$sample <- colnames(df)[1]
      colnames(df)[1] <- 'value'
      df$chromosome <- chr
      #plot <- ggplot(df) + aes(x = df[,2], y = df[,1]) + geom_line() + xlab('pos') + ylab(names(df)[1])
      #print(plot)
      facet_all <- rbind(facet_all, df)
    }
  }
    
    for(sample in unique(facet_all$sample)) {
      #sample <- '35V'
      sample_df <- facet_all[facet_all$sample == sample,]
      plot <- ggplot(sample_df) + aes(x = pos, y = value) + geom_line() + facet_wrap(sample~chromosome)
      print(plot)
    }
  dev.off()
  
  #return(df)
}

comp_df_int <- function(chr, ginf, founder_lookup_table, map_chr) {
  df <- as.data.frame(ginf[[chr]])
  df_chr <- apply(df, 2, function(x) founder_lookup_table[x]) # column level
  rownames(df_chr) <- rownames(df)
  df_chr <- as.data.frame(t(df_chr))
  df_chr$marker <- rownames(df_chr)
  df_chr_comp <- merge(df_chr, map_chr, by = 'marker')
  
  df_comp <- df_chr_comp[1:(ncol(df_chr_comp))] %>% select(-c(marker))
  
  return(df_comp)
}

# do this for each strain
genocode_comp_matrix <- function(map, ginf_qtl2, ginf_haploqa, founder_lookup_table, num_chr, file_gen, sample_type) {
  comp_matrix <- list()
  for (chr in num_chr) {
    ## qtl2
    print(chr)
    map_chr <- as.data.frame(map[[chr]]) %>% rename('pos' = 'map[[chr]]')
    map_chr$marker <- rownames(map_chr)

    haploqa_comp <- comp_df_int(chr, ginf_haploqa, founder_lookup_table, map_chr)
    qtl2_comp <- comp_df_int(chr, ginf_qtl2, founder_lookup_table, map_chr)
    
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
      comp_fp <- paste0(root, '/geno_comp_results_', sample_type)
      dir.create(comp_fp, showWarnings = FALSE) 
      write.csv(result_df, paste0(comp_fp, '/geno_comp_matrix_', chr, '.csv'), row.names = F)
    }
  }
  return(comp_matrix)
}


get_geno_pct_diff <- function(ginf_haploqa, ginf_qtl2, summary_df, num_chr, founder_lookup_table) {
  
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

# function to initiate the haplotype reconstruction pipeline
# @param sample_type (string) - type of sample, such as CC, DO, F2, etc.
# @param list_pheno (list) - list of phenotypes to be simulated
# @param qtl2_file_gen (True/False) - toggle to set whether to (re)generate qtl2 input files
# @param samples_gen (True/False) - toggle to set whether to (re)generate individual sample  files
#
# @return results (list of dataframes) - a list containing all objects necessary for computations
haplotype_reconstruction_pipeline <- function(sample_type, list_pheno, qtl2_file_gen, samples_gen) {
  
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
  
  ### Part 2 - convert results to qtl2 input format
  # implement function
  if (qtl2_file_gen == T) {
    file_output <- get_qtl2_input(data_dir, sample_type, annot_file, qtl2_dir, summary_df, list_pheno, ngen, founders_list, marker_type, exclude_list)
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
    
  
    # write out files
    write.csv(df_geno, paste0(qtl2_dir, '/test_geno.csv'), row.names = F)
    write.csv(df_gmap, paste0(qtl2_dir, '/test_gmap.csv'), row.names = F)
    write.csv(df_pmap, paste0(qtl2_dir, '/test_pmap.csv'), row.names = F)
    write.csv(df_pheno, paste0(qtl2_dir, '/test_pheno.csv'), row.names = F)
    write.csv(df_covar, paste0(qtl2_dir, '/test_covar.csv'), row.names = F)
    write.csv(df_foundergeno, paste0(qtl2_dir, '/test_foundergeno.csv'), row.names = F)
    write.csv(df_crossinfo, paste0(qtl2_dir, '/test_crossinfo.csv'), row.names = F)
    writeLines(control_file, control_fp)
  }
  
  ### get results for the sample
  results <- save_rds(rds_dir, qtl2_dir, sample_type, results_dir, data_dir, n_founders, founder_haplo_lookup)
  
  return(results)
}

get_raw_geno <- function(sample_type, chromosome) {
  #sample_type <- 'F2'
  #chromosome <- '1'
  ### config file
  config <- fread(paste0(root, '/annotations_config.csv'))
  
  ### Environment
  config_sample <- config[config$array_type == sample_type]
  data_dir <- config_sample$data_dir
  # data output directory
  data_dir <- file.path(root, config_sample$data_dir)
  
  map_df <- fread(file.path(root, config_sample$qtl2_dir, 'test_gmap.csv'))
  
  # annotation file
  annot_file <- config_sample$annot_file
  
  ## summary file
  summary_df_fp <- paste0(data_dir, '/', sample_type, '_summary.csv')
  summary_df <- fread(summary_df_fp)
  
  dir <- "https://raw.githubusercontent.com/kbroman/MUGAarrays/main/UWisc/" # always the same
  annot_df <- read.csv(paste0(dir, annot_file))
  # encode_genome function on geno and founder from qtl2 convert package - args: matrix of genotypes, matrixs of two alleles
  annot_encode_df <- annot_df %>% select(marker, chr, snp) %>% # select allele columns
    separate(snp, c("allele1", "allele2"), sep=cumsum(c(1)))
  
  # all txt file from directory
  data_files <- dir(data_dir, pattern = '\\.txt$', full.names = TRUE)
  ### combine all files
  df_raw <- rbindlist(lapply(data_files, read_sample_txt)) %>% filter(!sample_id %in% exclude_list) # save Y chrom here for gender plotting
  
  # clean up summary df
  sum_df <- summary_df %>% 
    select(ID, Sex, `% Het. Calls`, `% No Call`, Platform) %>% rename(sample_id = ID)
  sum_df$`% Het. Calls` <- as.numeric(gsub("%", "", sum_df$`% Het. Calls`))
  sum_df$`% No Call` <- as.numeric(gsub("%", "", sum_df$`% No Call`))
  
  # eliminate those that has no call > 10%
  ### use checkifnot to make sure the markers/SNP are in the annotation file (annotation file is the boss)
  df_all <- merge(df_raw, sum_df, by = 'sample_id', sort = F) %>% filter(`% No Call` < 10) %>% filter(!chromosome %in% c('0', 'Y', 'M')) %>% filter(!sample_id %in% exclude_list) 
  
  geno_sub <- df_all %>% select(sample_id, snp_id, allele1_fwd, allele2_fwd) %>%
    mutate(gene_exp = paste(allele1_fwd, allele2_fwd, sep = '')) %>%
    select(sample_id, snp_id, gene_exp) %>% rename(marker = snp_id)
  df_geno <- dcast(geno_sub, marker~sample_id, value.var="gene_exp")
  raw_geno <- merge(geno_sub, map_df, by = 'marker', sort = F) %>% filter(chr == chromosome) %>% arrange(pos)
  
  return(raw_geno)
  
}

plot_onegeno_test <- function(geno, geno1, map, ind=1, chr=NULL,
                              col=NULL, na_col="white",
                              swap_axes=FALSE,
                              border="black", shift=FALSE,
                              chrwidth=0.5, ...) {
  
  # ignore class of geno object
  geno <- unclass(geno)
  geno1 <- unclass(geno1)
  
  # drop all but the target individual
  if(length(ind)>1) {
    ind <- ind[1]
    warning("Only using the first individual")
  }
  
  # separate mom and dad
  for(i in seq_along(geno)) {
    geno[[i]] <- geno[[i]][ind,,,drop=FALSE]
    geno1[[i]] <- geno1[[i]][ind,,,drop=FALSE]
    
  }
  
  # always start at 0
  map <- lapply(map, function(a) a-min(a,na.rm=TRUE))
  
  plot_onegeno_internal <-
    function(geno, map, col=NULL, na_col="white",
             swap_axes=FALSE,
             border="black", bgcolor="gray90",
             chrwidth=0.5,
             xlab=NULL, ylab=NULL,
             xlim=NULL, ylim=NULL, las=1, xaxs=NULL, yaxs=NULL,
             mgp.x=c(2.6,0.5,0), mgp.y=c(2.6,0.5,0), mgp=NULL,
             hlines=NULL, hlines_col="white", hlines_lwd=1, hlines_lty=1,
             vlines=NULL, vlines_col="gray80", vlines_lwd=1, vlines_lty=1,
             ...)
    {
      dots <- list(...)
      
      nchr <- length(map)
      
      # margin parameters
      if(!is.null(mgp)) mgp.x <- mgp.y <- mgp
      
      if(is.null(xlab)) xlab <- "Chromosome"
      if(is.null(ylab)) ylab <- "Position"
      
      if(is.null(xlim)) xlim <- c(0.5, 52)
      if(is.null(ylim)) ylim <- rev(range(unlist(map), na.rm=TRUE))
      
      if(is.null(hlines)) hlines <- pretty(ylim)
      if(is.null(vlines)) vlines <- seq_len(52)
      
      if(is.null(xaxs)) xaxs <- "i"
      if(is.null(yaxs)) yaxs <- "r"
      
      # blank canvas
      plot(0, 0, type="n", xlab="", ylab="",
           xaxs=xaxs, yaxs=yaxs,
           xaxt="n", yaxt="n", xlim=xlim, ylim=ylim)#), ...)
      
      u <- par("usr")
      if(!is.null(bgcolor))
        rect(u[1], u[3], u[2], u[4], col=bgcolor, border=NA)
      
      # include axis labels?
      if(is.null(dots$xaxt)) dots$xaxt <- par("xaxt")
      if(is.null(dots$yaxt)) dots$yaxt <- par("yaxt")
      
      # add x axis unless par(xaxt="n")
      if(dots$xaxt != "n") {
        odd <- seq(1, nchr, by=2)
        axis(side=1, at=(odd*2.5), names(map)[odd],
             mgp=mgp.x, las=las, tick=FALSE)
        if(nchr > 1) {
          even <- seq(2, nchr, by=2)
          axis(side=1, at=(even*2.5), names(map)[even],
               mgp=mgp.x, las=las, tick=FALSE)
        }
      }
      # add y axis unless par(yaxt="n")
      if(dots$yaxt != "n") {
        axis(side=2, at=pretty(ylim), mgp=mgp.y, las=las, tick=FALSE)
      }
      
      
      # grid lines
      if(!(length(vlines)==1 && is.na(vlines))) {
        abline(v=vlines, col=vlines_col, lwd=vlines_lwd, lty=vlines_lty)
      }
      if(!(length(hlines)==1 && is.na(hlines))) {
        abline(h=hlines, col=hlines_col, lwd=hlines_lwd, lty=hlines_lty)
      }
      
      # x and y axis labels
      title(xlab=xlab, mgp=mgp.x)
      title(ylab=ylab, mgp=mgp.y)
      
      max_geno <- max(unlist(geno), na.rm=TRUE)
      if(is.null(col)) {
        if(max_geno <= 8) {
          col <- qtl2::CCcolors
        }
        else {
          warning("With ", max_geno, " genotypes, you need to provide the vector of colors; recycling some")
          col <- rep(qtl2::CCcolors, max_geno)
        }
      }
      else if(max_geno > length(col)) {
        warning("not enough colors; recycling them")
        col <- rep(col, max_geno)
      }
      ### chromosomes 1-19
      inc <- 0
      
      for(i in seq_len(nchr-1)) {
        #i <- 19
        g <- geno[[i]]
        g1 <- geno1[[i]]
        
        # if completely missing the second chr but not the first, treat as if we have just the one
        #   (this is a kludge to deal with males on X chr;
        #    really should use is_x_chr and is_female but we don't have it)
        this_chrwidth <- chrwidth
        if(!is.matrix(g) && !all(is.na(g[,,1])) && all(is.na(g[,,2]))) { # if g needs to be converted, g1 should be too
          g <- rbind(g[,,1]) # make it a row matrix
          g1 <- rbind(g1[,,1])
          this_chrwidth <- this_chrwidth/2
        }
        ### dataframe - g
        # rectangle shape
        x_lower_left <- i+0.75+inc
        x_higher_left <- (i+1.25+inc)-(this_chrwidth/2)
        x_lower_right <- i+1.25+inc
        x_higher_right <- (i+0.75+inc)+(this_chrwidth/2)
        
        rect(x_higher_left, min(map[[i]], na.rm=TRUE),
             x_lower_left, max(map[[i]], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[[i]], na.rm=TRUE),
             x_lower_right, max(map[[i]], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        
        # add geno colors
        addgenorect(g[1,,1], map[[i]], x_higher_left,x_lower_left,
                    col=col, swap_axes=swap_axes)
        addgenorect(g[1,,2], map[[i]], x_higher_right, x_lower_right,
                    col=col, swap_axes=swap_axes)
        # borders
        rect(x_higher_left, min(map[[i]], na.rm=TRUE),
             x_lower_left, max(map[[i]], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[[i]], na.rm=TRUE),
             x_lower_right, max(map[[i]], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        
        ### dataframe - g1
        # rectangle shape
        x_lower_left <- x_lower_left + 1
        x_higher_left <- x_higher_left + 1
        x_lower_right <-x_lower_right + 1
        x_higher_right <- x_higher_right + 1
        
        rect(x_higher_left, min(map[[i]], na.rm=TRUE),
             x_lower_left, max(map[[i]], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[[i]], na.rm=TRUE),
             x_lower_right, max(map[[i]], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        
        # add geno colors
        addgenorect(g1[1,,1], map[[i]], x_higher_left, x_lower_left,
                    col=col, swap_axes=swap_axes)
        addgenorect(g1[1,,2], map[[i]], x_higher_right, x_lower_right,
                    col=col, swap_axes=swap_axes)
        # borders
        rect(x_higher_left, min(map[[i]], na.rm=TRUE),
             x_lower_left, max(map[[i]], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[[i]], na.rm=TRUE),
             x_lower_right, max(map[[i]], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        inc <- inc + 1.5
        
        
        
      }
      
      ## plot x individually
      g <- geno[['X']]
      g1 <- geno1[['X']]
      if(!is.matrix(g) && !all(is.na(g[,,1])) && all(is.na(g[,,2]))) { # if g needs to be converted, g1 should be too
        g <- rbind(g[,,1]) # make it a row matrix
        g1 <- rbind(g1[,,1])
        this_chrwidth <- this_chrwidth/2
      }
      
      inc <- inc + 1
      
      
      if(is.matrix(g)) {
        print('x has only one')
        ### dataframe - g
        x_lower_left <- i+0.75+inc
        x_higher_left <- (i+1.25+inc)-(this_chrwidth/2)
        x_lower_right <- i+1.25+inc
        x_higher_right <- (i+0.75+inc)+(this_chrwidth/2)
        
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        
        addgenorect(g[1,], map[['X']], x_higher_left, x_lower_left,
                    col=col, swap_axes=swap_axes)
        
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        
        ### dataframe - g1
        x_lower_left <- x_lower_left + 1
        x_higher_left <- x_higher_left + 1
        x_lower_right <- x_lower_right + 1
        x_higher_right <- x_higher_right + 1
        
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        
        addgenorect(g1[1,], map[['X']], x_higher_left, x_lower_left,
                    col=col, swap_axes=swap_axes)
        
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        
        
      } else {
        # rectangle shape
        x_lower_left <- i+0.75+inc
        x_higher_left <- (i+1.25+inc)-(this_chrwidth/2)
        x_lower_right <- i+1.25+inc
        x_higher_right <- (i+0.75+inc)+(this_chrwidth/2)
        
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[['X']], na.rm=TRUE),
             x_lower_right, max(map[['X']], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        
        # add geno colors
        addgenorect(g[1,,1], map[['X']], x_higher_left,x_lower_left,
                    col=col, swap_axes=swap_axes)
        addgenorect(g[1,,2], map[['X']], x_higher_right, x_lower_right,
                    col=col, swap_axes=swap_axes)
        # borders
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[['X']], na.rm=TRUE),
             x_lower_right, max(map[['X']], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        
        ### dataframe - g1
        # rectangle shape
        x_lower_left <- x_lower_left + 1
        x_higher_left <- x_higher_left + 1
        x_lower_right <-x_lower_right + 1
        x_higher_right <- x_higher_right + 1
        
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[['X']], na.rm=TRUE),
             x_lower_right, max(map[['X']], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        
        # add geno colors
        addgenorect(g1[1,,1], map[['X']], x_higher_left, x_lower_left,
                    col=col, swap_axes=swap_axes)
        addgenorect(g1[1,,2], map[['X']], x_higher_right, x_lower_right,
                    col=col, swap_axes=swap_axes)
        # borders
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[['X']], na.rm=TRUE),
             x_lower_right, max(map[['X']], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        
      }
      
      
      
      
      box()
    }
  
  plot_onegeno_internal(geno, map, col=col, na_col=na_col,
                        swap_axes=swap_axes, border=border,
                        chrwidth=chrwidth, ...)
  
  
}

# add rectangles for the genotypes
addgenorect <- function(geno, map, x1, x2, col, swap_axes=FALSE) {
  intervals <- geno2intervals(geno, map)
  if(is.null(intervals) || nrow(intervals) < 1) return(NULL)
  
  for(i in seq_len(nrow(intervals))) {
    if(swap_axes) {
      rect(intervals[i,1], x1,
           intervals[i,2], x2,
           col=col[intervals[i,3]],
           border=NA, lend=1, ljoin=1)
    } else{
      rect(x1, intervals[i,1],
           x2, intervals[i,2],
           col=col[intervals[i,3]],
           border=NA, lend=1, ljoin=1)
    }
  }
}


# convert vector of integer genotypes to intervals with common genotypes
# (start, end, genotype)
geno2intervals <- function(geno, map) {
  if(all(is.na(geno))) return(NULL)
  
  stopifnot(length(geno) == length(map))
  
  # drop missing values
  map <- map[!is.na(geno)]
  geno <- geno[!is.na(geno)]
  
  d <- diff(geno)
  xo_int <- which(d != 0)
  
  data.frame(lo=map[c(1,xo_int+1)],
             hi=map[c(xo_int, length(map))],
             geno=geno[c(xo_int, length(map))])
  
}



xo_number_plot <- function(sample_result) {
  #sample_result <- cc_results
  sample_type <- toupper(strsplit(deparse(substitute(sample_result)), '_')[[1]][1])
  pos_haploqa_do <- locate_xo(sample_result[['ginf_haploqa']], sample_result[['map']])
  pos_qtl2_do <- sample_result[['pos']]
  
  
  xo_all <- data.frame()
  for (i in seq(1, length(pos_haploqa_do[[1]]))) { # should be the same for all chromosomes
    print(i)
    #i <- 1
    xo_haplo <- list()
    xo_qtl2 <- list()
    for (chr in seq(1,19)) {
      #chr <- 1
      xo_haplo <- c(xo_haplo, length(pos_haploqa_do[[chr]][[i]]))
      xo_qtl2 <- c(xo_haplo, length(pos_qtl2_do[[chr]][[i]]))
    }
    xo_haplo_total <- sum(unlist(xo_haplo)) # number of crossovers for one sample across entire genome (all chromosomes)
    xo_qtl2_total <- sum(unlist(xo_qtl2)) 
    df <- data.frame(sample = names(pos_haploqa_do$`1`)[i], total_sample_xo_haplo = xo_haplo_total, total_sample_xo_qtl2 = xo_qtl2_total)
    xo_all <- rbind(xo_all, df)
  }
  
  xo_all$total_sample_xo_haplo <- as.numeric(xo_all$total_sample_xo_haplo)
  xo_all$total_sample_xo_qtl2 <- as.numeric(xo_all$total_sample_xo_qtl2)
  #sd_col <- sapply(xo_all[c(2, 3)], sd)
  high_lim <- max(xo_all[,c(2,3)])
  low_lim <- min(xo_all[,c(2,3)])
  
  plot <- ggplot(xo_all) + aes(x = total_sample_xo_haplo, y = total_sample_xo_qtl2) + geom_point() +
    xlab('total number of crossovers on autosomes - HaploQA') + ylab('total number of crossovers on autosomes - qtl2') + 
    xlim(low_lim, high_lim) + ylim(low_lim, high_lim) + ggtitle(paste0('Sample: ', sample_type, ' - each point represents one individual')) +
    geom_abline(intercept = 0, slope = 1, linetype="dashed")
  
  print(plot)
  
}


loc_xo_distance_plot <- function(sample_result) {
  #sample_result <- cc_results
  num_chr <- c((1:19),"X")
  sample_type <- toupper(strsplit(deparse(substitute(sample_result)), '_')[[1]][1])
  pos_haploqa_do <- locate_xo(sample_result[['ginf_haploqa']], sample_result[['map']])
  pos_qtl2_do <- sample_result[['pos']]
  x_loc_diff <- location_xo_comp(pos_qtl2_do, pos_haploqa_do, num_chr)
  
  
  loc_xo_dist_all <- data.frame()
  for (i in seq(1, length(pos_haploqa_do[[1]]))) { # should be the same for all chromosomes
    print(i)
    #i <- 1
    xo_dist_sample <- list()
    for (chr in seq(1,19)) {
      #chr <- 1
      xo_dist_sample <- c(xo_dist_sample, x_loc_diff[[chr]][[i]])
    }
    df <- data.frame(sample = names(pos_haploqa_do$`1`)[i], loc_xo_dist = mean(unlist(xo_dist_sample)))
    loc_xo_dist_all <- rbind(loc_xo_dist_all, df)
  }
  
  st_err <- sd(loc_xo_dist_all$loc_xo_dist)
  
  plot <- ggplot(loc_xo_dist_all) + aes(x = seq(1, length(sample)), y = loc_xo_dist) + geom_point() + 
    xlab('Sample index') + ylab('mean of crossover distance') +
    ggtitle(paste0('Mean of crossover distance between qtl2 and HaploQA - sample ', sample_type)) +
    geom_errorbar(aes(ymin=loc_xo_dist-st_err, ymax=loc_xo_dist+st_err), width=.2, position=position_dodge(.9)) 
  
  print(plot)
}



