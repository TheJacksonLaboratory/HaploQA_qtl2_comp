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


# source the function script
root <- dirname(getSourceEditorContext()$path)
source(paste0(root,"/input_data_prep_functions.R"))

sample_type <- 'DO' # CC/DO/MiniMUGA/GigaMUGA/BXD
list_pheno <- c('WBC', 'NEUT')

config <- fread(paste0(root, '/annotations_config.csv'))
config_sample <- config[config$array_type == sample_type]

sample_url <- config_sample$url
ngen <- 3

## results directory
results_dir <- file.path(root, 'results')
dir.create(results_dir, showWarnings = FALSE) 

# filepaths to save any rds file
rds_dir <- file.path(results_dir, 'RDS')
dir.create(rds_dir, showWarnings = FALSE)

qtl2_file_gen <- FALSE
num_chr <- c((1:19),"X")

## text of control file
control_file <- '{
  "description": "HaploQA data - qtl2 test run",
  "crosstype": "genail2",
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
  "alleles": ["A", "B", "C", "D", "E", "F", "G", "H"],
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
}'

### fix founder contributions
# founder for bxd comes from https://haploqa.jax.org/tag/UNC_Villena_GIGMUGV01_20141012_FinalReport.html
do_results <- get_sample_result(sample_type, sample_url, results_dir, list_pheno, qtl2_file_gen = F, samples_gen = F, control_file, num_chr, ngen)



### 1. genoprob_to_alleleprob - condense into founder prob (save to rds)
founder_all_codes <- colnames(do_results[['pr']]$`X`) # take the X chromosome - this one has everything
all_map_codes <- seq(1:length(founder_all_codes))
#founder_all_lookup <- setNames(all_map_codes, founder_all_codes) # make maxmarg
founder_all_rev_lookup <- setNames(founder_all_codes, all_map_codes)

test_geno <- do_results[['ginf_haploqa']]
do_qtl2_test <- genoprob_to_alleleprob(do_results[['pr']])
map <- do_results[['map']]

geno_align <- function(df) {
  founder_codes <- LETTERS[seq(1,8)]
  col_diff <- setdiff(founder_codes, colnames(df))
  dum_array <- array(0, dim = c(nrow(df), length(col_diff)))
  colnames(dum_array) <- col_diff
  df <- cbind(df, dum_array)
  df <- df[ ,founder_codes]
  
  return(df)
}


temp <- do_qtl2_test$`19`
temp <- temp * temp

temp_save <- temp

temp <- apply(temp, c(1,3), sum)
temp <- sqrt(temp)

#do_haplo_test <- list()
do_haplo_test <- list()
for (i in num_chr) {
  #i <- '19'
  print(i)
  df <- test_geno[[i]]
  
  do_haplo_test[[i]] <- array(dim = c(nrow(df), 8, ncol(df)))
  rownames(do_haplo_test[[i]]) <- rownames(df)
  colnames(do_haplo_test[[i]]) <- LETTERS[seq(1,8)]
  for (col in seq(1, ncol(df))) {
    #col <- 328
    print(col)
    
    df_test <- data.frame(sample_id = rownames(df), col_val = df[,c(col)])
    df_test[,2] <- as.data.frame(apply(df_test[c(2)], 2, function(x) founder_all_rev_lookup[x]))
    df_test <- df_test %>% separate(col_val, c("allele1", "allele2"), sep=cumsum(c(1)))
    
    # haplo
    a <- table(df_test[,c(1,2)])
    b <- table(df_test[,c(1,3)])
    
    if(length(colnames(b)) < 8) {
      b <- geno_align(b)
    } 
    if((length(colnames(a)) < 8)) {
      a <- geno_align(a)
    }
    
    haplo_comp <- (a + b) / 2
    
    do_haplo_test[[i]][,,col] <- haplo_comp
    
  }
}

saveRDS(do_haplo_test, paste0(rds_dir, 'haploqa_founder_prob_', sample_type, '.rds'))


cos_sim_all <- list()
for (i in num_chr) {
  #i <- '1'
  df <- test_geno[[i]]
  map_i <- map[[i]]
  print(i)
  cos_sim_chr <- list()
  for (col in seq(1, ncol(df))) {
    #col <- 3
    map_val <- as.numeric(map_i[col])
    # qtl2
    #print(col)
    qtl2_comp <- do_qtl2_test[[i]][,,col]
    ### plug the pre-calculated square roots
    haplo_comp <- do_haplo_test[[i]][,,col]
    
    a = rowSums(qtl2_comp * haplo_comp)
    b = (sqrt(rowSums(qtl2_comp * qtl2_comp))) * (sqrt(rowSums(haplo_comp * haplo_comp)))
    
    cos_sim_list = a/b

    cos_marker_df <- data.frame(sample_id = rownames(haplo_comp), cos_sim = unlist(cos_sim_list), pos = map_val)
    cos_marker_df$chromosome <- i
    cos_sim_chr[[col]] <- as.data.table(cos_marker_df)
  }
  cos_sim_all[[i]] <- rbindlist(cos_sim_chr)
}  

cos_sim_all <- rbindlist(cos_sim_all)

saveRDS(cos_sim_all, paste0(rds_dir, '/cos_sim_', sample_type, '.rds'))

for (sample in unique(cos_sim_all$sample_id)) {
  sample <- '6UY'
  sample_df <- cos_sim_all[(cos_sim_all$sample_id == sample),]
  sample_df$chromosome <- factor(sample_df$chromosome, num_chr)
  ggplot(sample_df) + aes(x = pos, y = cos_sim) + geom_line() + ggtitle(paste0('Sample ID: ', sample)) + 
    facet_wrap(.~chromosome)
  
}




exclude_list <- c('http://haploqa-dev.jax.org//sample/595e3b183ac36a089e040c82.html',
                  'http://haploqa-dev.jax.org//sample/595e3b1a3ac36a089e040c83.html',
                  'http://haploqa-dev.jax.org//sample/595e3b1c3ac36a089e040c84.html',
                  'http://haploqa-dev.jax.org//sample/595e3b1f3ac36a089e040c86.html',
                  'http://haploqa-dev.jax.org//sample/595e3b243ac36a089e040c88.html',
                  'http://haploqa-dev.jax.org//sample/595e3b253ac36a089e040c89.html',
                  'http://haploqa-dev.jax.org//sample/595e3b273ac36a089e040c8a.html', 
                  'http://haploqa-dev.jax.org//sample/595e3b2a3ac36a089e040c8b.html',
                  'http://haploqa-dev.jax.org//sample/595e3b2c3ac36a089e040c8c.html',
                  'http://haploqa-dev.jax.org//sample/595e3b303ac36a089e040c8e.html',
                  'http://haploqa-dev.jax.org//sample/595e3b323ac36a089e040c8f.html',
                  'http://haploqa-dev.jax.org//sample/595e3b333ac36a089e040c90.html',
                  'http://haploqa-dev.jax.org//sample/595e3b353ac36a089e040c91.html',
                  'http://haploqa-dev.jax.org//sample/595e3b3a3ac36a089e040c94.html',
                  'http://haploqa-dev.jax.org//sample/595e3b3d3ac36a089e040c96.html',
                  'http://haploqa-dev.jax.org//sample/595e3b3f3ac36a089e040c97.html',
                  'http://haploqa-dev.jax.org//sample/595e3b413ac36a089e040c98.html',
                  'http://haploqa-dev.jax.org//sample/595e3b433ac36a089e040c99.html'
                  )



### Part 3 - haplotype reconstructions
### CC
cc_results <- save_rds(rds_dir, qtl2_dir, sample_type, results_dir, data_dir)



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

### DO
do_results <- get_sample_result(sample_type, sample_url, results_dir, list_pheno, qtl2_file_gen = F, samples_gen = F, control_file, num_chr, ngen)


# lookup tables
founder_all_codes <- colnames(do_results[['pr']]$`X`) # take the X chromosome - this one has everything
all_map_codes <- seq(1:length(founder_all_codes))
founder_all_lookup <- setNames(all_map_codes, founder_all_codes) # make maxmarg
founder_all_rev_lookup <- setNames(founder_all_codes, all_map_codes)

## location crossover comparison
pos_haploqa_do <- locate_xo(do_results[['ginf_haploqa']], do_results[['map']])
pos_qtl2_do <- do_results[['pos']]
x_loc_diff <- location_xo_comp(pos_qtl2_do, pos_haploqa_do, num_chr)

## error matrix
err_comp(do_results[['pr']], do_results[['map']], num_chr, results_dir, sample_type)

## genocode comparison matrix
comp_matrix <- genocode_comp_matrix(do_results[['map']], do_results[['ginf']], do_results[['ginf_haploqa']], founder_all_rev_lookup, num_chr, file_gen = T, sample_type)


### percent genomic difference
chr_pct <- get_geno_pct_diff(do_results[['ginf']], do_results[['ginf_haploqa']], do_results[['summary']], num_chr, founder_all_rev_lookup)


  