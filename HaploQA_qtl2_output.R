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
library(readr)

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

do_results <- haplotype_reconstruction_pipeline('DO', list_pheno, qtl2_file_gen = F, samples_gen = F, truth_model = F)

do_truth_results <- haplotype_reconstruction_pipeline('DO', list_pheno, qtl2_file_gen = F, samples_gen = F, truth_model = T)

## count pct of missing values in sample geno

cc_results <- haplotype_reconstruction_pipeline('CC', list_pheno, qtl2_file_gen = F, samples_gen = F, truth_model = F)
cc_truth_results <- haplotype_reconstruction_pipeline('CC', list_pheno, qtl2_file_gen = T, samples_gen = F, truth_model = T)

bxd_results <- haplotype_reconstruction_pipeline('BXD', list_pheno, qtl2_file_gen = F, samples_gen = F, truth_model = F)
bxd_truth_results <- haplotype_reconstruction_pipeline('BXD', list_pheno, qtl2_file_gen = T, samples_gen = F, truth_model = T)

f2_results <- haplotype_reconstruction_pipeline('F2', list_pheno, qtl2_file_gen = F, samples_gen = F, truth_model = F)
f2_truth_results <- haplotype_reconstruction_pipeline('F2', list_pheno, qtl2_file_gen = T, samples_gen = F, truth_model = T)


## write results
get_genoprob_example <- function(sample_result, sample_type) {
  num_chr <- c((1:19),"X")
  ### config file
  config <- fread(paste0(root, '/annotations_config.csv'))
  
  ### Environment
  config_sample <- config[config$array_type == sample_type]
  ## create a data directory for qtl2 input data
  qtl2_dir <- file.path(root, config_sample$qtl2_dir)
  foundergeno <- fread(file.path(qtl2_dir, 'test_foundergeno.csv'))
  map <- fread(file.path(qtl2_dir, 'test_gmap.csv'))
  founder_geno_map <- merge(foundergeno, map, by = 'marker')
  genoprob_all <- list()
  for (chromosome in num_chr) {
    #chromosome <- num_chr[1]
    print(chromosome)
    founder_map_chr <- founder_geno_map %>% filter(chr == chromosome)
    #print(founder_map_chr)
    chr_genoprob <- sample_result[['pr']][[chromosome]]
    
    chr_list <- list()
    for(i in dimnames(chr_genoprob)[[3]]) {
      #i <- dimnames(chr_genoprob)[[3]][1]
      m <- chr_genoprob[,,i] %>% as.data.frame()
      m$marker <- i
      m$sample_name <- rownames(m)
      chr_list[[i]] <- m
    }
    chr_df <- rbindlist(chr_list)
    df <- merge(chr_df, founder_map_chr, by = 'marker')
    genoprob_all[[chromosome]] <- df
  }
  genoprob_all <- rbindlist(genoprob_all, fill = TRUE)
  return(genoprob_all)
}


cc_genoprob <- get_genoprob_example(cc_results, 'CC')
do_genoprob <- get_genoprob_example(do_results, 'DO')
bxd_genoprob <- get_genoprob_example(bxd_results, 'BXD')
f2_genoprob <- get_genoprob_example(f2_results, 'F2')

saveRDS(cc_genoprob, file = file.path(root, 'genoprob_cc_founder.rds'))
saveRDS(do_genoprob, file = file.path(root, 'genoprob_do_founder.rds'))
saveRDS(bxd_genoprob, file = file.path(root, 'genoprob_bxd_founder.rds'))
saveRDS(f2_genoprob, file = file.path(root, 'genoprob_f2_founder.rds'))


minimuga_results <- haplotype_reconstruction_pipeline('MiniMUGA', list_pheno, qtl2_file_gen = T, samples_gen = T)


qtl2_diff <- list()
haploqa_diff <- list()
do_comp_csv <- geno_all_comp('DO', do_results, do_truth_results)
do_comp_csv <- rbindlist(do_comp_csv)
#ind_qtl2_diff <- which(do_comp_csv$qtl2_calls_truth != do_comp_csv$qtl2_calls)
#ind_haploqa_diff <- which(do_comp_csv$qtl2_calls_truth != do_comp_csv$haplo_diplotype)
#do_comp_csv$is_different <- NA
#do_comp_csv[ind,]$is_different <- TRUE
qtl2_diffs <- do_comp_csv %>% group_by(sample_id) %>% summarise(qtl2_pct_diff = sum(qtl2_calls_truth != qtl2_calls)/n()) %>% as.data.frame()
haploqa_diffs <- do_comp_csv %>% group_by(sample_id) %>% summarise(haploqa_pct_diff = sum(qtl2_calls_truth != haplo_diplotype)/n()) %>% as.data.frame()
qtl2_haploqa_diffs <- do_comp_csv %>% group_by(sample_id) %>% summarise(haploqa_pct_diff = sum(qtl2_calls != haplo_diplotype)/n()) %>% as.data.frame()

sample_diffs <- merge(merge(qtl2_diffs, haploqa_diffs, by = 'sample_id'), qtl2_haploqa_diffs, by = 'sample_id')
write.csv(sample_diffs, file.path(root, 'do_truth_comp.csv'), row.names = F)

### use this comparison for shiny rows
cc_comp_csv <- geno_all_comp('CC', cc_results, cc_truth_results)
cc_comp_csv <- rbindlist(cc_comp_csv)
#ind <- which(cc_comp_csv$qtl2_calls_truth != cc_comp_csv$qtl2_calls | cc_comp_csv$qtl2_calls_truth != cc_comp_csv$haplo_diplotype)
#cc_comp_csv$is_different <- NA
#cc_comp_csv[ind,]$is_different <- TRUE
#ind_qtl2_diff <- which(cc_comp_csv$qtl2_calls_truth != cc_comp_csv$qtl2_calls)
#ind_haploqa_diff <- which(cc_comp_csv$qtl2_calls_truth != cc_comp_csv$haplo_diplotype)
#qtl2_diff <- append(qtl2_diff, as.numeric(length(ind_qtl2_diff)/nrow(cc_comp_csv)))
#haploqa_diff <- append(haploqa_diff, as.numeric(length(ind_haploqa_diff)/nrow(cc_comp_csv)))
qtl2_diffs <- cc_comp_csv %>% group_by(sample_id) %>% summarise(qtl2_pct_diff = sum(qtl2_calls_truth != qtl2_calls)/n()) %>% as.data.frame()
haploqa_diffs <- cc_comp_csv %>% group_by(sample_id) %>% summarise(haploqa_pct_diff = sum(qtl2_calls_truth != haplo_diplotype)/n()) %>% as.data.frame()
sample_diffs <- merge(qtl2_diffs, haploqa_diffs, by = 'sample_id')
write.csv(sample_diffs, file.path(root, 'cc_truth_comp.csv'), row.names = F)


bxd_comp_csv <- geno_all_comp('BXD', bxd_results, bxd_truth_results) 
bxd_comp_csv <- rbindlist(bxd_comp_csv)
#haploqa_diff <- append(haploqa_diff, as.numeric(length(ind_haploqa_diff)/nrow(cc_comp_csv)))
qtl2_diffs <- bxd_comp_csv %>% group_by(sample_id) %>% summarise(qtl2_pct_diff = sum(qtl2_calls_truth != qtl2_calls)/n()) %>% as.data.frame()
haploqa_diffs <- bxd_comp_csv %>% group_by(sample_id) %>% summarise(haploqa_pct_diff = sum(qtl2_calls_truth != haplo_diplotype)/n()) %>% as.data.frame()
sample_diffs <- merge(qtl2_diffs, haploqa_diffs, by = 'sample_id')
write.csv(sample_diffs, file.path(root, 'bxd_truth_comp.csv'), row.names = F)
#ind <- which(bxd_comp_csv$qtl2_calls_truth != bxd_comp_csv$qtl2_calls | bxd_comp_csv$qtl2_calls_truth != bxd_comp_csv$haplo_diplotype)
#bxd_comp_csv$is_different <- NA
#bxd_comp_csv[ind,]$is_different <- TRUE
#rows_diff <- append(rows_diff, as.numeric(length(ind)/nrow(bxd_comp_csv)))
#ind_qtl2_diff <- which(bxd_comp_csv$qtl2_calls_truth != bxd_comp_csv$qtl2_calls)
#ind_haploqa_diff <- which(bxd_comp_csv$qtl2_calls_truth != bxd_comp_csv$haplo_diplotype)
#qtl2_diff <- append(qtl2_diff, as.numeric(length(ind_qtl2_diff)/nrow(bxd_comp_csv)))
#haploqa_diff <- append(haploqa_diff, as.numeric(length(ind_haploqa_diff)/nrow(bxd_comp_csv)))


f2_comp_csv <- geno_all_comp('F2', f2_results, f2_truth_results) 
f2_comp_csv <- rbindlist(f2_comp_csv)
#ind <- which(f2_comp_csv$qtl2_calls_truth != f2_comp_csv$qtl2_calls | f2_comp_csv$qtl2_calls_truth != f2_comp_csv$haplo_diplotype)
#f2_comp_csv$is_different <- NA
#f2_comp_csv[ind,]$is_different <- TRUE
#rows_diff <- append(rows_diff, as.numeric(length(ind)/nrow(f2_comp_csv)))
ind_qtl2_diff <- which(f2_comp_csv$qtl2_calls_truth != f2_comp_csv$qtl2_calls)
ind_haploqa_diff <- which(f2_comp_csv$qtl2_calls_truth != f2_comp_csv$haplo_diplotype)
qtl2_diff <- append(qtl2_diff, as.numeric(length(ind_qtl2_diff)/nrow(f2_comp_csv)))
haploqa_diff <- append(haploqa_diff, as.numeric(length(ind_haploqa_diff)/nrow(f2_comp_csv)))



sum_pct_markers <- data.frame(sample_type = c('DO', 'CC', 'BXD', 'F2'), qtl2_pct_marker_diff = unlist(qtl2_diff), haploqa_pct_marker_diff = unlist(haploqa_diff))


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
write.csv(sample_geno_8, file.path(root, 'genocomp_acgt', '6VE_do_geno_comp.csv'))

### get all crossovers, 20 markers above and below
### run this for one sample, all chromosomes, add chr and pos information
## make into one file



#do_haplo_test <- list()



cos_sim_plot(do_results, rds_dir, results_dir, 'DO')
cos_sim_plot(f2_results, rds_dir, results_dir, 'F2')
cos_sim_plot(cc_results, rds_dir, results_dir, 'CC')
cos_sim_plot(bxd_results, rds_dir, results_dir, 'BXD')


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


minimuga_url <- 'http://haploqa-dev.jax.org/tag/MiniMUGA.html'
html_table <- read_html(minimuga_url) %>% html_nodes("a") %>% html_attr("href")
url_list <- paste0(url_domain, html_table[grepl('/sample', html_table)])
summary_df <- sample_summary_scrape(haploqa_html, url_list, marker_type)
ind_mini <- summary_df[summary_df$`Haplotype Candidate` == 'False',]$`ID`
ind_mini <- ind_mini[ind_mini!='JXX'] # screentime error, skip this one
ind_mini <- ind_mini[ind_mini!='JY9'] # only 1 contributing strain, AIL incompatible, skip this one as well
# use J/NJ
## screentime error: JXX
for (ind in ind_mini) {
  #ind <- 'JXX'
  print(ind)
  sample_haplotype_reconstruction('MiniMUGA', ind, samples_gen = T, qtl2_file_gen = T)
}

  