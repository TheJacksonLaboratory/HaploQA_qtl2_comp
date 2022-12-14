# contains all functions that aligns haploqa results with that of qtl2 and save the results
# this is necessary because the computed haploqa results needs to be in the same format as qtl2 (with genoprob, phased geno, etc.)
### function directory
## 1. maxmarg_sim - simulate maxmarg phasing computations for haploqa
## 2. save_rds - save all results into rds


library(rstudioapi)
root <- dirname(getSourceEditorContext()$path)


source(paste0(root,"/Scripts/utils.R"))

# function to simulate the qtl2 maxmarg phasing computations
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
# @return ginf_haploqa (list of matrices) - phased genotypes indexed by chromosome, same as output of maxmarg() from qtl2
maxmarg_sim <- function(summary_df, data_dir, num_chr, founder_all_rev_lookup, cross, qtl2_dir, sample_type, n_founders, founder_haplo_lookup, pr) {
  
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
  
  #pr <- cross[['pr']]
  
  ### join sex into this df - for males, take first letter, then paste Y after it
  #founder_all_lookup <- setNames(names(founder_all_rev_lookup), founder_all_rev_lookup) # reverse the lookup table
  #founder_all_lookup <- founder_all_lookup[names(founder_all_lookup) %in% unique(haploqa_maxmarg$haplotype)]
  #founder_all_lookup <- setNames(seq(1, length(founder_all_lookup)), names(founder_all_lookup))
  #haploqa_maxmarg[,4] <- as.data.frame(apply(haploqa_maxmarg[,c(4)], 2, function(x) founder_all_lookup[x]))
  ### subset the lookup table names to just what is in haploqa_maxmarg$haplotype
  ### then assign numbers seq(1, length())
  
  
  ## save per chromosome
  ginf_haploqa <- list()
  for (i in num_chr) {
     #i <- '1'
    
    # founder table
    geno_codes <- colnames(pr[[i]])
    #founder_all_lookup <- setNames(names(founder_all_rev_lookup), founder_all_rev_lookup)
    founder_all_lookup <- setNames(seq(1, length(geno_codes)), geno_codes)
    print(founder_all_lookup)
    
    haploqa_maxmarg_chr <- haploqa_maxmarg %>% filter(chromosome == i)
    
    haploqa_maxmarg_chr[,4] <- as.data.frame(apply(haploqa_maxmarg_chr[,c(4)], 2, function(x) founder_all_lookup[x]))
    
    print(i)
    df <- haploqa_maxmarg_chr %>% dcast(sample_id ~ snp_id, value.var = 'haplotype') %>% as.data.frame()
    df[,2:ncol(df)] <- mutate_all(df[,2:ncol(df)], function(x) as.numeric(as.character(x)))
    rownames(df) <- df$sample_id
    
    df <- df %>% select(-c(sample_id)) 
    col_order <- names(cross$gmap[[i]])
    df <- df[,match(col_order, names(df))]
    #df <- df[,col_order]
    ginf_haploqa[[i]] <- df 
  }
  
  attr(ginf_haploqa, "crosstype") <- paste0('genail', length(founder_haplo_lookup))
  
  return(ginf_haploqa)
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
qtl2_metrics_comp <- function(rds_dir, qtl2_dir, sample_type, results_dir, data_dir, n_founders, founder_haplo_lookup, truth_model) {
  num_chr <- c((1:19),"X")
  print(sample_type)
  if (truth_model == T) {
    sample_type = paste0(sample_type, '_truth')
  }
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
  summary_df <- fread(paste0(data_dir, '/', strsplit(sample_type, '_')[[1]][1], '_summary.csv'))
  cross_list[['summary']] <- summary_df
  
  cross_list[['data_dir']] <- data_dir
  
  cross_list[['founder_lookup']] <- founder_haplo_lookup
  
  founder_lookup_table <- fread(file.path(root, 'founder_lookup_table.csv'))
  founder_all_rev_lookup <- setNames(founder_lookup_table$founder_codes, founder_lookup_table$founder_id)
  
  if (truth_model == F) { # only need HaploQA for genail models
    # simulate maxmarg to calculate haploqs ginf 
    ginf_haploqa_fp <- paste0(rds_dir, "/ginf_haploqa_", sample_type, ".rds")
    if (!file.exists(ginf_haploqa_fp)) {
      print('haploqa ginf file does not exist, running calculations from maxmarg simulation')
      ginf_haploqa <- maxmarg_sim(summary_df, data_dir, num_chr, founder_all_rev_lookup, cross = cross_list[['cross']], qtl2_dir, sample_type, n_founders, founder_haplo_lookup, pr)
      ## save to rds file
      cross_list[['ginf_haploqa']] <- ginf_haploqa
      saveRDS(ginf_haploqa, file = ginf_haploqa_fp)
    } else {
      print('haploqa ginf file exists, reading in file')
      ginf_haploqa <- readRDS(ginf_haploqa_fp)
      cross_list[['ginf_haploqa']] <- ginf_haploqa
    }
    
    # simulate guess_phase for haploqa to calculate phased geno
    ph_geno_haploqa_fp <- paste0(rds_dir, "/ph_geno_haploqa_", sample_type, ".rds")
    if (!file.exists(ph_geno_haploqa_fp)) {
      print('haploqa phased genotype file does not exist, running calculations from guess phase simulation')
      ph_geno_haploqa <- guess_phase(cross_list[['cross']], ginf_haploqa)
      cross_list[['ph_geno_haploqa']] <- ph_geno_haploqa
      saveRDS(ph_geno_haploqa, file = ph_geno_haploqa_fp)
    } else {
      print('haploqa phased genotype file exists, reading in file')
      ph_geno_haploqa <- readRDS(ph_geno_haploqa_fp)
      cross_list[['ph_geno_haploqa']] <- ph_geno_haploqa
    }
  }
  
  
  return(cross_list)
}

