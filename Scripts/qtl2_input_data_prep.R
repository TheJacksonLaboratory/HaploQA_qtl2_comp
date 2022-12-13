### Contains all functions to generate qtl2 input data

library(rstudioapi)
root <- dirname(getSourceEditorContext()$path)

source(paste0(root,"/Scripts/input_data_retrieval.R"))
source(paste0(root,"/Scripts/utils.R"))


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
  ### cut markers more than 10%
  # merge with annotations
  df_geno <- df_geno %>% filter(marker %in% unique(annot_encode_df$marker)) ### replace with a stopif assert line here
  # put alleles in order of marker as shown in geno dataframe
  geno_encode_annot <- merge(df_geno, annot_encode_df, by = 'marker', all.x = T) %>% select(allele1, allele2) # let it break if there's NAs
  rownames(df_geno) <- df_geno$marker
  df_geno_encoded <- as.data.frame(encode_geno(df_geno[,-1], geno_encode_annot))
  df_geno_encoded$marker <- rownames(df_geno)
  # reorder columns
  df_geno_encoded <- df_geno_encoded %>% select(marker, everything())
  colnames(df_geno_encoded) <- colnames(df_geno)
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
qtl2_foundergeno <- function(df, founder_url, url_list, founders_dict, annot_encode_df, founders_list = NULL, marker_order, founder_haplo_lookup, sample_type){
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
    founders_list <- names(founder_haplo_lookup)# if not 8 founder strains, custom select letters
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
    df_founders_encoded$marker <- rownames(test_foundergeno)
    df_founders_encoded <- df_founders_encoded %>% filter(marker %in% marker_order) %>% arrange(factor(marker, levels = marker_order))
    df_founders_encoded <- df_founders_encoded %>% select(marker, everything())
    
  }
  
  
  return(df_founders_encoded)
}

# function to generate cross info data for qtl2 input
# @param summary_df (data.frame) - table as shown on HaploQA sample website, same as output of sample_summary_scrape
# columns of summary_df: ID (sample IDs), secondary ID, Haplotype Candidate (T/F), Strain Name, Sex, % Het Calls, % Hom Calls, % No Call, % Concordance, Sample filepath
#
# @return founders_total (data.frame) - cross info data with sample IDs and mapped strain info for qtl2 input
# columns of founders_total: id, A, B, C, D, E, F, G, H
qtl2_ci <- function(summary_df, ngen, founder_haplo_lookup) {
  founder_list <- LETTERS[seq(1,8)]
  # join with summary table and CC file in karl's github
  ### if no founder info, filter out.
  # match the strain names with the cc cross info csv
  ci_sum <- summary_df %>% select(`ID`) %>% unique() %>% rename(id = ID)
  ci_sum$ngen <- ngen
  ci_sum[ ,c(founder_list[1:length(founder_haplo_lookup)])] <- 1
  df_crossinfo <- ci_sum # read_cross does not allow NAs here
  
  return(df_crossinfo)
}

# function to get individual founder strain data from HaploQA, if not exist
# @param founder_url (string) - url to main page of sample that contains founder individuals, for now it's 'https://haploqa.jax.org/tag/UNC_Villena_GIGMUGV01_20141012_FinalReport.html'
# @param url_list (list) - list of urls to main page of founder individuals, retrieved from founder_url
#
# @return founders_total (data.frame) - combined dataframe with all founder individual files
# columns of founders_total: sample_id, original_sample_id, snp_id, chromosome, position_bp, allele1_fwd, allele2_fwd, haplotype1, haplotype2
get_founder_data <- function(founder_url, url_list, sample_type, data_dir, founder_filename) {
  # get length to keep track of progress
  iter_len <- length(url_list)
  # loop through urls to extract individual files
  #founders_total = data.frame() # container
  founders_sample_dir <- paste0(data_dir, '/', sample_type, '_founders/')
  dir.create(founders_sample_dir, showWarnings = FALSE) 
  
  fp_founders <- file.path(root, founder_filename)
  if (file.exists(fp_founders)) {
    print('founder sample data exists')
    founders_total <- fread(fp_founders)
  } else { 
    print('founder sample data does not exist')
    inc = 0
    st_error_urls <- list()
    for (url in url_list) {
      file <- sample_individual_scrape(url, url_domain) # call function from another script
      inc = inc+1 # track progress
      print(paste0('working on file ', inc, '/', iter_len, ' ', file))
      # screen time errors
      skip <- FALSE
      tryCatch(as.data.frame(content(GET(file), encoding="UTF-8")), error = function(e) { skip <<- TRUE})
      if(!skip) { 
        df_temp <- as.data.frame(content(GET(file), encoding="UTF-8"))
      } else {
        print(paste0('skipped ', url, ' due to screentime error'))
        st_error_urls <- c(st_error_urls, url)
        next
      } 
      
      file_name <- unlist(strsplit(file, '/'))[6]
      print('Writing to directory')
      #founders_total <- rbind(founders_total,df_temp)
      
      GET(file, write_disk(file.path(founders_sample_dir, file_name), overwrite = TRUE), show_col_types = FALSE)
      
    }
    
    founder_data_files <- dir(founders_sample_dir, pattern = '\\.txt$', full.names = TRUE)
    ### combine all files
    founders_total <- rbindlist(lapply(founder_data_files, read_sample_txt))
    write.csv(founders_total, fp_founders, row.names = F)
  }
  
  print(paste0('urls skipped: ', st_error_urls))
  
  return(founders_total)
  
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
get_qtl2_input <- function(data_dir, sample_type, annot_file, qtl2_output_dir, summary_df, list_pheno, ngen, founders_list, marker_type, exclude_list, founder_url, founder_filename) {
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
  df_raw <- rbindlist(lapply(data_files, read_sample_txt)) #%>% filter(!sample_id %in% exclude_list) # save Y chrom here for gender plotting
  
  ## exclude the ones with bad string names
  if(sample_type == 'CC') {
    summary_df <- summary_df %>% filter(grepl('CC', `Strain Name`)) 
    
  }
  
  # clean up summary df
  sum_df <- summary_df %>% 
    select(ID, Sex, `% Het. Calls`, `% No Call`, Platform) %>% rename(sample_id = ID)
  sum_df$`% Het. Calls` <- as.numeric(gsub("%", "", sum_df$`% Het. Calls`))
  sum_df$`% No Call` <- as.numeric(gsub("%", "", sum_df$`% No Call`))
  
  # eliminate those that has no call > 10%
  ### use checkifnot to make sure the markers/SNP are in the annotation file (annotation file is the boss)
  df_all <- merge(df_raw, sum_df, by = 'sample_id') %>% 
    filter(`% No Call` < 10) %>% filter(!chromosome %in% c('0', 'Y', 'M')) %>% 
    filter(!sample_id %in% exclude_list) 
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
  url_domain <- 'http://haploqa-dev.jax.org/'
  # dictionary
  founders_dict <- fread(paste0(root, '/founder_strains_table.csv')) %>% rename(`Strain Name` = founder_strain)
  founder_strains <- unique(founders_dict$`Strain Name`) # get the names of founder strains
  # get html and list of url to sample data
  html_table <- read_html(founder_url) %>% html_nodes("a") %>% html_attr("href")
  url_list <- paste0(url_domain, html_table[grepl('/sample', html_table)])
  
  # get founder data for gigamuga
  if (sample_type %in% c('CC', 'DO', 'BXD', 'F2', 'MURGIGV01')) {
    founders_total <- get_founder_data(founder_url, url_list, sample_type, data_dir, paste0(marker_type, '_founders.csv'))
    #write.csv(founders_total, fp_founders, row.names = F)
    
    df_founders_encoded <- qtl2_foundergeno(founders_total, founder_url, url_list, founders_dict, annot_encode_df, founders_list, marker_order, founder_haplo_lookup, sample_type)
    file_output[[6]] <- df_founders_encoded 
  }
  
  
  print('founder geno data done')
  
  # cross info
  df_crossinfo <- qtl2_ci(summary_df, ngen, founder_haplo_lookup)
  
  # store in output
  file_output[[7]] <- df_crossinfo 
  print('cross info data done')
  
  return(file_output)
}

