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

# source the function script
root <- dirname(getSourceEditorContext()$path)
source(paste0(root,"/input_data_prep_functions.R"))

# simulated geno
list_pheno <- c('WBC', 'NEUT')

### GigaMUGA pipelines
do_results <- haplotype_reconstruction_pipeline('DO', list_pheno, qtl2_file_gen = F, samples_gen = F, truth_model = F)
cc_results <- haplotype_reconstruction_pipeline('CC', list_pheno, qtl2_file_gen = F, samples_gen = F, truth_model = F)
bxd_results <- haplotype_reconstruction_pipeline('BXD', list_pheno, qtl2_file_gen = F, samples_gen = F, truth_model = F)
f2_results <- haplotype_reconstruction_pipeline('F2', list_pheno, qtl2_file_gen = F, samples_gen = F, truth_model = F)

## GigaMUGA truth models
do_truth_results <- haplotype_reconstruction_pipeline('DO', list_pheno, qtl2_file_gen = F, samples_gen = F, truth_model = T)
cc_truth_results <- haplotype_reconstruction_pipeline('CC', list_pheno, qtl2_file_gen = T, samples_gen = F, truth_model = T)
bxd_truth_results <- haplotype_reconstruction_pipeline('BXD', list_pheno, qtl2_file_gen = T, samples_gen = F, truth_model = T)
f2_truth_results <- haplotype_reconstruction_pipeline('F2', list_pheno, qtl2_file_gen = T, samples_gen = F, truth_model = T)


### MiniMUGA pipeline
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

#### write results
### genoprobs
cc_genoprob <- get_genoprob_example(cc_results, 'CC')
do_genoprob <- get_genoprob_example(do_results, 'DO')
bxd_genoprob <- get_genoprob_example(bxd_results, 'BXD')
f2_genoprob <- get_genoprob_example(f2_results, 'F2')

saveRDS(cc_genoprob, file = file.path(root, 'genoprob_cc_founder.rds'))
saveRDS(do_genoprob, file = file.path(root, 'genoprob_do_founder.rds'))
saveRDS(bxd_genoprob, file = file.path(root, 'genoprob_bxd_founder.rds'))
saveRDS(f2_genoprob, file = file.path(root, 'genoprob_f2_founder.rds'))


### percentage of marker different between qtl2, haploqa and truth
shiny_pct_fp <- file.path(results_dir, 'shiny_pct_csvs') # directory
dir.create(shiny_pct_fp, showWarnings = F)
## DO
do_comp_csv <- geno_all_comp('DO', do_results, do_truth_results)
write.csv(do_comp_csv, file.path(shiny_pct_fp, 'do_truth_comp.csv'), row.names = F)
## CC
cc_comp_csv <- geno_all_comp('CC', cc_results, cc_truth_results)
write.csv(cc_comp_csv, file.path(shiny_pct_fp, 'cc_truth_comp.csv'), row.names = F)
## BXD
bxd_comp_csv <- geno_all_comp('BXD', bxd_results, bxd_truth_results) 
write.csv(bxd_comp_csv, file.path(shiny_pct_fp, 'bxd_truth_comp.csv'), row.names = F)
## F2
f2_comp_csv <- geno_all_comp('F2', f2_results, f2_truth_results) 



### geno allele comparisons
# plot two individuals per cross
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


### cosine similarity plots
cos_sim_plot(do_results, rds_dir, results_dir, 'DO')
cos_sim_plot(f2_results, rds_dir, results_dir, 'F2')
cos_sim_plot(cc_results, rds_dir, results_dir, 'CC')
cos_sim_plot(bxd_results, rds_dir, results_dir, 'BXD')



## location crossover comparison
pos_haploqa_cc <- locate_xo(cc_results[['ginf_haploqa']], cc_results[['map']])
pos_qtl2_cc <- cc_results[['pos']]
x_loc_diff <- location_xo_comp(pos_qtl2_cc, pos_haploqa_cc, num_chr)

## error matrix
err_comp(cc_results[['pr']], cc_results[['map']], num_chr, results_dir, sample_type)
err_comp(do_results[['pr']], do_results[['map']], num_chr, results_dir, sample_type)

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

## genocode comparison matrix
# lookup tables
founder_all_codes <- colnames(cc_results[['pr']]$`X`) # take the X chromosome - this one has everything
all_map_codes <- seq(1:length(founder_all_codes))
founder_all_lookup <- setNames(all_map_codes, founder_all_codes) # make maxmarg
founder_all_rev_lookup <- setNames(founder_all_codes, all_map_codes)

comp_matrix_do <- genocode_comp_matrix(do_results[['map']], do_results[['ginf']], do_results[['ginf_haploqa']], founder_all_rev_lookup, num_chr, file_gen = T, sample_type)
comp_matrix_cc <- genocode_comp_matrix(cc_results[['map']], cc_results[['ginf']], cc_results[['ginf_haploqa']], founder_all_rev_lookup, num_chr, file_gen = T, sample_type)


### percent genomic difference
chr_pct_do <- get_geno_pct_diff(do_results[['ginf']], do_results[['ginf_haploqa']], do_results[['summary']], num_chr, founder_all_rev_lookup)
chr_pct_cc <- get_geno_pct_diff(cc_results[['ginf']], cc_results[['ginf_haploqa']], cc_results[['summary']], num_chr, founder_all_rev_lookup)




  