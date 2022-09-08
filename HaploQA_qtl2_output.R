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

cc_results <- haplotype_reconstruction_pipeline('CC', list_pheno, qtl2_file_gen = F, samples_gen = F)

bxd_results <- haplotype_reconstruction_pipeline('BXD', list_pheno, qtl2_file_gen = F, samples_gen = F)

f2_results <- haplotype_reconstruction_pipeline('F2', list_pheno, qtl2_file_gen = F, samples_gen = F)


minimuga_results <- haplotype_reconstruction_pipeline('MiniMUGA', list_pheno, qtl2_file_gen = T, samples_gen = T)



geno_all_comp <- function(sample, sample_type) {
#sample <- '8RS'
founder_lookup_table <- fread(file.path(root, 'founder_lookup_table.csv'))
founder_all_rev_lookup <- setNames(founder_lookup_table$founder_codes, founder_lookup_table$founder_id)

config_sample <- config[config$array_type == sample_type]
qtl2_dir <- file.path(root, config_sample$qtl2_dir)
df_cov <- fread(file.path(qtl2_dir, 'test_covar.csv')) %>% rename('sample_id' = 'id')


haploqa_diplotype <- get_haplotypes(summary_df, data_dir)
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
  #chr <- '1'
  print(chr)
    ## raw geno
  raw_geno_acgt <- get_raw_geno(sample_type, chromosome = chr)
  raw_geno_chr_acgt <- raw_geno_acgt %>% filter(sample_id == sample)
  
  f2_haploqa_mom <- as.data.frame(apply(as.data.frame(f2_results[['ph_geno_haploqa']][[chr]][,,1][sample,]), 1, function(x) founder_all_lookup[x]))
  f2_haploqa_mom$marker <- rownames(f2_haploqa_mom)
  f2_haploqa_mom <- f2_haploqa_mom %>% rename(f2_haploqa_mom = 1)
  f2_haploqa_dad <- as.data.frame(apply(as.data.frame(f2_results[['ph_geno_haploqa']][[chr]][,,2][sample,]), 1, function(x) founder_all_lookup[x]))
  f2_haploqa_dad$marker <- rownames(f2_haploqa_dad)
  f2_haploqa_dad <- f2_haploqa_dad %>% rename(f2_haploqa_dad = 1)
  f2_geno <- merge(f2_haploqa_mom, f2_haploqa_dad, by = 'marker')
  f2_geno$haplo_diplotype <- paste(f2_geno$f2_haploqa_mom, f2_geno$f2_haploqa_dad, sep='')
  f2_geno <- f2_geno %>% select(marker, haplo_diplotype)
  
  
  f2_haploqa_mom <- as.data.frame(apply(as.data.frame(f2_results[['ph_geno']][[chr]][,,1][sample,]), 1, function(x) founder_all_lookup[x]))
  f2_haploqa_mom$marker <- rownames(f2_haploqa_mom)
  f2_haploqa_mom <- f2_haploqa_mom %>% rename(f2_haploqa_mom = 1)
  f2_haploqa_dad <- as.data.frame(apply(as.data.frame(f2_results[['ph_geno']][[chr]][,,2][sample,]), 1, function(x) founder_all_lookup[x]))
  f2_haploqa_dad$marker <- rownames(f2_haploqa_dad)
  f2_haploqa_dad <- f2_haploqa_dad %>% rename(f2_haploqa_dad = 1)
  f2_geno_qtl2 <- merge(f2_haploqa_mom, f2_haploqa_dad, by = 'marker')
  f2_geno_qtl2$qtl2_calls <- paste(f2_geno_qtl2$f2_haploqa_mom, f2_geno_qtl2$f2_haploqa_dad, sep='')
  f2_geno_qtl2 <- f2_geno_qtl2 %>% select(marker, qtl2_calls)
  
  
  
  df_geno_all_acgt <- merge(merge(raw_geno_chr_acgt, f2_geno, by = 'marker', sort = F), f2_geno_qtl2, by = 'marker', sort = F) %>% select(marker, sample_id, chr, pos, everything())
  df_geno_all_acgt[df_geno_all_acgt$haplo_diplotype != df_geno_all_acgt$qtl2_calls,]
  df_geno_all_acgt$is_crossover <- NA
  
  ind <- which(df_geno_all_acgt$haplo_diplotype != df_geno_all_acgt$qtl2_calls)
  
  list_xo <- list()
  list <- split(ind, cumsum(c(1, diff(ind) != 1)))
  for (i in seq(1, length(list))) {
    list_xo <- c(list[[i]][1], list_xo)
  }
  list_xo <- sort(unlist(list_xo))
  
  rows <- unlist(lapply(list_xo, function(x) (x-20):(x+20)))
  rows <- rows[rows > 0]
  df_geno_all_acgt[list_xo,]$is_crossover <- 'True'
  df_geno_xo <- df_geno_all_acgt[rows,]
  
  df_geno_xo[is.na(df_geno_xo$is_crossover),]$is_crossover <- 'False'
  
  df_geno_xo$original_row <- rows # always 41 - 20 above, 20 below, and crossover itself
    
  df_geno_all_chr[[chr]] <- df_geno_xo
  }

  return(df_geno_all_chr)
}


 
sample_geno <- geno_all_comp('8RS', 'F2')




### get all crossovers, 20 markers above and below
### run this for one sample, all chromosomes, add chr and pos information
## make into one file
## then start minimu

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

do_haplo_test <- readRDS(paste0(rds_dir, '/haploqa_founder_prob_', sample_type, '.rds'))

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
sample_type <- 'DO'
cos_sim_all <- readRDS(paste0(rds_dir, '/cos_sim_', sample_type, '.rds'))
df <- cos_sim_all %>% group_by(sample_id) %>% summarise(mean_chr = mean(cos_sim)) %>% as.data.frame()
#df$chromosome <- factor(df$chromosome, num_chr)
### make into line plot
ggplot(df) + aes(y = mean_chr, x = seq(1, length(mean_chr))) + geom_line() + xlab('sample_index (1-277)') # take mean of entire genome
pdf(paste0(results_dir, '/cos_sim_', sample_type, '.pdf'))
for (sample in unique(cos_sim_all$sample_id)) {
  sample <- '6UY'
  print(sample)
  sample_df <- cos_sim_all[(cos_sim_all$sample_id == sample),]
  sample_df %>% group_by(sample_id) %>% summarise(mean_chr = mean(cos_sim)) %>% as.data.frame()
  sample_df$chromosome <- factor(sample_df$chromosome, num_chr)
  plot <- ggplot(sample_df) + aes(x = pos, y = cos_sim) + geom_line() + ggtitle(paste0('Sample ID: ', sample)) + 
    facet_wrap(.~chromosome)
  print(plot)
  
}

dev.off()



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


  