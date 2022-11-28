##### HaploQA data reformatting - Launch script
### used for qtl2 input
### annotation config saved in annotations_config.csv
### founder strain dictionary saved in founder_strains_table.csv

### output files for haplotype reconstruction


## everything that starts from BJSJL

# take one founder out (cast), use genail 7
# 6UY DO
# run reconstruction, plot onegeno and run error lod
# find in python where haploqa generates the plots

# 1. DO genail 7, do error lod computations
# 2. document
# 3. haploqa codes
# 4. 
sample_res <- sample_haplotype_reconstruction('MURGIGV01', 'PX', samples_gen = F, qtl2_file_gen = T)


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

# 2, 4, 8, 16 founders
# pipeline for one GigaMUGA mutant sample
# system.time()

# source the function script
root <- dirname(getSourceEditorContext()$path)
source(paste0(root,"/input_data_prep_functions.R"))

# results directory
results_dir <- file.path(root, 'results')
dir.create(results_dir, showWarnings = FALSE) 

# simulated geno
list_pheno <- c('WBC', 'NEUT')

### GigaMUGA pipelines
do_results <- haplotype_reconstruction_pipeline('DO', list_pheno, qtl2_file_gen = F, samples_gen = F)
cc_results <- haplotype_reconstruction_pipeline('CC', list_pheno, qtl2_file_gen = F, samples_gen = F)
bxd_results <- haplotype_reconstruction_pipeline('BXD', list_pheno, qtl2_file_gen = F, samples_gen = F)
f2_results <- haplotype_reconstruction_pipeline('F2', list_pheno, qtl2_file_gen = F, samples_gen = F)

## GigaMUGA truth models
do_truth_results <- truth_model_reconstruction('DO', list_pheno, qtl2_file_gen = T, samples_gen = F)
cc_truth_results <- truth_model_reconstruction('CC', list_pheno, qtl2_file_gen = T, samples_gen = F)
bxd_truth_results <- truth_model_reconstruction('BXD', list_pheno, qtl2_file_gen = T, samples_gen = F)
f2_truth_results <- truth_model_reconstruction('F2', list_pheno, qtl2_file_gen = T, samples_gen = F)


png(file = file.path(root, 'do_qtl2_genoplot.png'), width = 1000, height = 800, res = 110)
qtl2::plot_onegeno(do_results[['ph_geno']], do_results[['map']], ind = 1, shift = TRUE, main = paste0('Geno-wide genotypes from qtl2, individual ', 1), sub = 'Sample: DO')
legend(13.5, 65, legend = c(names(qtl2::CCcolors)), title = 'Founder Colors', fill = qtl2::CCcolors)
dev.off()


png(file = file.path(root, 'cc_qtl2_genoplot.png'), width = 1000, height = 800, res = 110)
qtl2::plot_onegeno(cc_results[['ph_geno']], cc_results[['map']], ind = 1, shift = TRUE, main = paste0('Geno-wide genotypes from qtl2, individual ', 1), sub = 'Sample: CC') 
legend(13.5, 65, legend = c(names(qtl2::CCcolors)), title = 'Founder Colors', fill = qtl2::CCcolors)
dev.off()

png(file = file.path(root, 'bxd_qtl2_genoplot.png'), width = 1000, height = 800, res = 110)
qtl2::plot_onegeno(bxd_results[['ph_geno']], bxd_results[['map']], ind = 1, shift = TRUE, main = paste0('Geno-wide genotypes from qtl2, individual ', 1), sub = 'Sample: F3')
dev.off()

png(file = file.path(root, 'f2_qtl2_genoplot.png'), width = 1000, height = 800, res = 110)
qtl2::plot_onegeno(f2_results[['ph_geno']], f2_results[['map']], ind = 1, shift = TRUE, main = paste0('Geno-wide genotypes from qtl2, individual ', 1), sub = 'Sample: F2') 
dev.off()


### MiniMUGA pipeline
minimuga_url <- 'http://haploqa-dev.jax.org/tag/MiniMUGA.html'
url_domain <- 'http://haploqa-dev.jax.org/'
marker_type <- 'MiniMUGA'
html_table <- read_html(minimuga_url) %>% html_nodes("a") %>% html_attr("href")
url_list <- paste0(url_domain, html_table[grepl('/sample', html_table)])
summary_df <- sample_summary_scrape(read_html(minimuga_url), url_list, marker_type)
ind_mini <- summary_df[summary_df$`Haplotype Candidate` == 'False',]$`ID`
ind_mini <- ind_mini[ind_mini!='JXX'] # screentime error, skip this one
ind_mini <- ind_mini[ind_mini!='JY9'] # only 1 contributing strain, AIL incompatible, skip this one as well
# use J/NJ
## screentime error: JXX
ind_mini <- ind_mini[ind_mini!='JY6'] # temporary - JXU, JY3, JY5
minimuga_results <- list()
mini_pct_res <- list()
for (ind in ind_mini) {
  #ind <- 'JXN'
  print(ind)
  sample_mini_res <- sample_haplotype_reconstruction('MiniMUGA', ind, samples_gen = F, qtl2_file_gen = F)
  
  df <- minimuga_geno_comp(sample_mini_res, ind)
  
  minimuga_results[[ind]] <- sample_mini_res
  mini_pct_res[[ind]] <- df
}


### rename misleading names like save_rds
# revisit structures

mini_pct_res <- rbindlist(mini_pct_res)

mini_plot <- mini_pct_res %>% filter(chromosome != 'X') %>% group_by(sample_id) %>% summarise(qtl2_pct_diff = sum(haplo_diplotype != qtl2_calls)/n()) %>% as.data.frame()
plot_df <- melt(mini_plot)
plot_df$variable <- 'MiniMUGA'
library(ggplot2)    
ggplot(plot_df) + geom_boxplot(aes(x = variable, y = value)) + ggtitle('MiniMUGA - percentage of marker difference between genail and haploqa on autosomes') + xlab('model')


df <- minimuga_geno_comp(minimuga_results[['JXU']], 'JXU')
df1 <- minimuga_geno_comp(minimuga_results[['JXN']], 'JXN')
df2 <- minimuga_geno_comp(minimuga_results[['JXQ']], 'JXQ')

minimuga_geno_comp <- function(sample_results, sample_name) {
  num_chr <- c((1:19),"X")
  #sample_name <- 'JXU'
  sample_type <- 'MiniMUGA'
  #sample_results <- minimuga_results[['JXU']]
  founder_lookup_table <- fread(file.path(root, 'founder_lookup_table.csv'))
  founder_codes_dict <- setNames(founder_lookup_table$founder_codes, founder_lookup_table$founder_id)
  
  geno_codes <- colnames(sample_mini_res[['pr']]$X)
  founder_all_rev_lookup <- setNames(geno_codes, seq(1, length(geno_codes)))
  
  config <- fread(paste0(root, '/annotations_config.csv'))
  config_sample <- config[config$array_type == 'MiniMUGA']
  qtl2_dir <- file.path(root, config_sample$qtl2_dir)
  ### sample folder
  sample_dir <- file.path(qtl2_dir, sample_name)
  
  df_cov <- fread(file.path(sample_dir, 'test_covar.csv')) %>% rename('sample_id' = 'id')
  data_dir <- paste0(root, '/', config_sample$data_dir)
  summary_df <- fread(paste0(data_dir, '/', sample_type, '_summary.csv'))
  
  haploqa_diplotype <- get_haplotypes(summary_df, data_dir) %>% filter(sample_id == sample_name)
  
  unique_haplotypes <- unique(haploqa_diplotype[,c(haplotype1, haplotype2)])
  founder_haplo_lookup <- setNames(LETTERS[seq(1, length(unique(unique_haplotypes)))], unique(unique_haplotypes))
  
  raw_geno_df <- haploqa_diplotype %>% select(sample_id, snp_id, haplotype1, haplotype2, chromosome) %>% unique()
  raw_geno_df[,c(3,4)] <- as.data.frame(apply(raw_geno_df[,c(3,4)], 2, function(x) founder_haplo_lookup[x]))
  raw_geno_df <- merge(raw_geno_df, df_cov, by = 'sample_id')
  raw_geno_df[(raw_geno_df$chromosome == 'X') & (raw_geno_df$Sex == 'male')]$haplotype2 <- 'Y'
  raw_geno_df$haplotype <- paste(raw_geno_df$haplotype2, raw_geno_df$haplotype1, sep='')
  raw_geno_df <- raw_geno_df %>% select(sample_id, snp_id, chromosome, haplotype)
  raw_geno_df[,4] <- lapply(raw_geno_df[,4], function(col) {
    rev_col = stri_reverse(col)
    ifelse(rev_col %in% founder_codes_dict, rev_col, col)
  })
  
  founder_all_lookup <- setNames(LETTERS[seq(1, length(unique(raw_geno_df$haplotype)))], seq(1, length(unique(raw_geno_df$haplotype))))
  
  df_geno_all_chr <- list()
  for (chr in num_chr) {
    #chr <- 'X'
    
    print(chr)
    ## raw geno
    
    map_df <- fread(file.path(sample_dir, 'test_gmap.csv')) %>% rename(chromosome = chr)
    annot_file <- config_sample$annot_file
    dir <- "https://raw.githubusercontent.com/kbroman/MUGAarrays/main/UWisc/" # always the same
    annot_df <- read.csv(paste0(dir, annot_file))
    
    annot_encode_df <- annot_df %>% select(marker, chr, snp) %>% # select allele columns
      separate(snp, c("allele1", "allele2"), sep=cumsum(c(1)))
    
    
    ## raw geno codes
    # all txt file from directory
    data_files <- dir(data_dir, pattern = '\\.txt$', full.names = TRUE)
    ### combine all files
    df_raw <- rbindlist(lapply(data_files, read_sample_txt)) %>% filter(sample_id == sample_name)
    sum_df <- summary_df %>% filter(ID == sample_name) %>% rename(sample_id = ID)
    df_all <- merge(df_raw, sum_df, by = 'sample_id', sort = F) %>% filter(!chr %in% c('0', 'Y', 'M')) #%>% filter(!sample_id %in% exclude_list) 
    geno_sub <- df_all %>% select(sample_id, snp_id, allele1_fwd, allele2_fwd) %>%
      mutate(gene_exp = paste(allele1_fwd, allele2_fwd, sep = '')) %>%
      select(sample_id, snp_id, gene_exp) %>% rename(marker = snp_id)
    df_geno <- dcast(geno_sub, marker~sample_id, value.var="gene_exp")
    raw_geno_chr_acgt <- merge(geno_sub, map_df, by = 'marker', sort = F) %>% filter(chromosome == chr) %>% arrange(pos)
  
    ## qtl2 codes
    qtl2_mini_mom <- as.data.frame(apply(as.data.frame(sample_results[['ph_geno']][[chr]][,,1]), 1, function(x) founder_all_lookup[x]))
    qtl2_mini_dad <- as.data.frame(apply(as.data.frame(sample_results[['ph_geno']][[chr]][,,2]), 1, function(x) founder_all_lookup[x]))
    qtl2_mini_dad[is.na(qtl2_mini_dad)] <- 'Y'
    sample_geno <- data.frame(Map(paste, qtl2_mini_mom, qtl2_mini_dad, MoreArgs = list(sep = "")), check.names = F) %>% rename('qtl2_calls' = 1)
    sample_geno$marker <- names(sample_results[['ph_geno']][[chr]][,,1])
    
    
    ## haploqa
    #if(sample_name == 'JXU') {
     # haplo_mini_mom <- as.data.frame(apply(as.data.frame(sample_results[['ph_geno_haploqa']][[chr]][,,1]), 1, function(x) founder_all_lookup[x])) %>% t()
    #}
    haplo_mini_mom <- as.data.frame(apply(as.data.frame(sample_results[['ph_geno_haploqa']][[chr]][,,1]), 1, function(x) founder_all_lookup[x]))
    haplo_mini_dad <- as.data.frame(apply(as.data.frame(sample_results[['ph_geno_haploqa']][[chr]][,,2]), 1, function(x) founder_all_lookup[x]))
    haplo_mini_dad[is.na(haplo_mini_dad)] <- 'Y'
    sample_geno_haplo <- data.frame(Map(paste, haplo_mini_mom, haplo_mini_dad, MoreArgs = list(sep = "")), check.names = F) %>% rename('haplo_diplotype' = 1)
    sample_geno_haplo$marker <- names(sample_results[['ph_geno_haploqa']][[chr]][,,1])
    
    mini_df <- merge(merge(sample_geno_haplo, sample_geno, by = c('marker'), sort = F), raw_geno_chr_acgt, by = c('marker'), sort = F) %>%
      select(sample_id, marker, chromosome, pos, gene_exp, everything())
    
    mini_df[,6:ncol(mini_df)] <- lapply(mini_df[,6:ncol(mini_df)], function(col) {
      rev_col = stri_reverse(col)
      ifelse(rev_col %in% founder_codes_dict, rev_col, col)
    })
    
    df_geno_all_chr[[chr]] <- mini_df
    
  }
  
  df_geno_all_chr <- rbindlist(df_geno_all_chr)
  
  return(df_geno_all_chr)
}

sample_mini_res[['ph_geno']]
qtl2_mini_mom <- as.data.frame(apply(as.data.frame(sample_mini_res[['ph_geno']][[chr]][,,1]), 1, function(x) founder_all_lookup[x]))
qtl2_mini_dad <- as.data.frame(apply(as.data.frame(sample_mini_res[['ph_geno']][[chr]][,,2]), 1, function(x) founder_all_lookup[x]))
sample_geno <- data.frame(Map(paste, qtl2_mini_mom, qtl2_mini_dad, MoreArgs = list(sep = "")), check.names = F) %>% rename('qtl2_calls' = 1)
sample_geno$marker <- names(sample_mini_res[['ph_geno']][[chr]][,,1])


haplo_mini_mom <- as.data.frame(apply(as.data.frame(sample_mini_res[['ph_geno_haploqa']][[chr]][,,1]), 1, function(x) founder_all_lookup[x]))
haplo_mini_dad <- as.data.frame(apply(as.data.frame(sample_mini_res[['ph_geno_haploqa']][[chr]][,,2]), 1, function(x) founder_all_lookup[x]))
sample_geno_haplo <- data.frame(Map(paste, haplo_mini_mom, haplo_mini_dad, MoreArgs = list(sep = "")), check.names = F) %>% rename('haplo_diplotype' = 1)
sample_geno_haplo$marker <- names(sample_mini_res[['ph_geno_haploqa']][[chr]][,,1])

mini_df <- merge(sample_geno_haplo, sample_geno, by = c('marker'), sort = F)
merge(mini_df, raw_geno_df, on - c(''))

founder_lookup_table <- fread(file.path(root, 'founder_lookup_table.csv'))
founder_codes_dict <- setNames(founder_lookup_table$founder_codes, founder_lookup_table$founder_id)

geno_codes <- colnames(sample_mini_res[['pr']]$X)
founder_all_rev_lookup <- setNames(geno_codes, seq(1, length(geno_codes)))

config <- fread(paste0(root, '/annotations_config.csv'))
config_sample <- config[config$array_type == sample_type]
qtl2_dir <- file.path(root, config_sample$qtl2_dir)
df_cov <- fread(file.path(sample_dir, 'test_covar.csv')) %>% rename('sample_id' = 'id')
data_dir <- paste0(root, '/', config_sample$data_dir)
summary_df <- fread(paste0(data_dir, '/', sample_type, '_summary.csv'))

haploqa_diplotype <- get_haplotypes(summary_df, data_dir)
unique_haplotypes <- unique(haploqa_diplotype[,c(haplotype1, haplotype2)])
founder_haplo_lookup <- setNames(LETTERS[seq(1, length(unique(unique_haplotypes)))], unique(unique_haplotypes))

raw_geno_df <- haploqa_diplotype %>% select(sample_id, snp_id, haplotype1, haplotype2, chromosome) %>% unique()
raw_geno_df[,c(3,4)] <- as.data.frame(apply(raw_geno_df[,c(3,4)], 2, function(x) founder_haplo_lookup[x]))
raw_geno_df <- merge(raw_geno_df, df_cov, by = 'sample_id')
raw_geno_df[(raw_geno_df$chromosome == 'X') & (raw_geno_df$Sex == 'male')]$haplotype2 <- 'Y'
raw_geno_df$haplotype <- paste(raw_geno_df$haplotype2, raw_geno_df$haplotype1, sep='')
raw_geno_df <- raw_geno_df %>% select(sample_id, snp_id, chromosome, haplotype)
raw_geno_df[,4] <- lapply(raw_geno_df[,4], function(col) {
  rev_col = stri_reverse(col)
  ifelse(rev_col %in% founder_codes_dict, rev_col, col)
})

founder_all_lookup <- setNames(LETTERS[seq(1, length(unique(raw_geno_df$haplotype)))], seq(1, length(unique(raw_geno_df$haplotype))))

#haploqa_mom$marker <- colnames(truth_results[['ph_geno_haploqa']][[chr]][,,1])

#### write results
### genoprobs
saveRDS(do_results[['pr']], file = file.path(root, 'genoprob_do.rds'))
saveRDS(cc_results[['pr']], file = file.path(root, 'genoprob_cc.rds'))
saveRDS(bxd_results[['pr']], file = file.path(root, 'genoprob_bxd.rds'))
saveRDS(f2_results[['pr']], file = file.path(root, 'genoprob_f2.rds'))



### percentage of marker different between qtl2, haploqa and truth
shiny_pct_fp <- file.path(results_dir, 'shiny_pct_csvs') # directory
dir.create(shiny_pct_fp, showWarnings = F)
## DO
do_comp_csv <- geno_all_comp('DO', do_results, do_truth_results)
write.csv(do_comp_csv, file.path(shiny_pct_fp, 'do_truth_geno.csv'), row.names = F)
do_comp_pct_csv <- marker_comp_diff(do_comp_csv)
write.csv(do_comp_pct_csv, file.path(shiny_pct_fp, 'do_truth_comp.csv'), row.names = F)
#write.csv(do_comp_csv, file.path(shiny_pct_fp, 'do_truth_comp.csv'), row.names = F)
## CC
cc_comp_csv <- geno_all_comp('CC', cc_results, cc_truth_results)
write.csv(cc_comp_csv, file.path(shiny_pct_fp, 'cc_truth_geno.csv'), row.names = F)
cc_comp_pct_csv <- marker_comp_diff(cc_comp_csv)
write.csv(cc_comp_pct_csv, file.path(shiny_pct_fp, 'cc_truth_comp.csv'), row.names = F)

#write.csv(cc_comp_csv, file.path(shiny_pct_fp, 'cc_truth_comp.csv'), row.names = F)
## BXD
bxd_comp_csv <- geno_all_comp('BXD', bxd_results, bxd_truth_results) 
write.csv(bxd_comp_csv, file.path(shiny_pct_fp, 'bxd_truth_geno.csv'), row.names = F)
bxd_comp_pct_csv <- marker_comp_diff(bxd_comp_csv)
write.csv(bxd_comp_pct_csv, file.path(shiny_pct_fp, 'bxd_truth_comp.csv'), row.names = F)

## F2
f2_comp_csv <- geno_all_comp('F2', f2_results, f2_truth_results) 
write.csv(f2_comp_csv, file.path(shiny_pct_fp, 'f2_truth_geno.csv'), row.names = F)
f2_comp_pct_csv <- marker_comp_diff(f2_comp_csv)
write.csv(f2_comp_pct_csv, file.path(shiny_pct_fp, 'f2_truth_comp.csv'), row.names = F)



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


  