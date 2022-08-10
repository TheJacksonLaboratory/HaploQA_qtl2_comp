library(shiny)
library(qtl2)
library(dplyr)
library(data.table)
library(rstudioapi)
library(stringi)
library(ggplot2)

### requirements:
# 1. config file
# 2. gmap 
# 3. summary_df

# source the function script
root <- dirname(getSourceEditorContext()$path)
source(paste0(root,"/input_data_prep_functions.R"))

# config file
config <- fread(paste0(root, '/annotations_config.csv'))
## results directory
results_dir <- file.path(root, 'results')
dir.create(results_dir, showWarnings = FALSE) 
# to save rds
rds_dir <- file.path(results_dir, 'RDS')
dir.create(rds_dir, showWarnings = FALSE) 

all_results_gen = T
num_chr <- c((1:19),"X")

#### RDS functions
sample_types <- c('CC', 'DO') # right now supports CC, DO, GigaMUGA, MiniMUGA
results_all_file <- file.path(rds_dir, 'results_all.rds')
if(all_results_gen == T) {
  list_results_all <- rds_loop(sample_types, config, results_dir, rds_dir, num_chr) 
  print('saving final results to RDS file')
  saveRDS(list_results_all, file = results_all_file)
} else {
  list_results_all <- readRDS(results_all_file)
}

cc_results <- list_results_all[['CC']]
do_results <- list_results_all[['DO']]

## list components: summary, data_dir, cross, map (same for qtl2 and haploqa)
## pr, ginf, ph_geno, pos (qtl2 only)
## ginf_haploqa, ph_geno_haploqa

## attributes from cc_test
total_ind <- 277 # take this automatically from cc_test?

##############

# lookup tables
founder_codes <- c('AA', 'BB', 'CC', 'DD', 'EE', 'FF', 'GG', 'HH')
single_codes <- sapply(founder_codes, function(x) substr(x, 1, 1))
map_codes <- c(1, 2, 3, 4, 5, 6, 7, 8)
founder_lookup <- setNames(founder_codes, map_codes)
founder_reverse <- setNames(map_codes, single_codes)
founders_dict <- fread(paste0(root, '/founder_strains_table.csv'))
founder_haplo_lookup <- setNames(founders_dict$letter, founders_dict$founder_strain)

founder_all_codes <- colnames(pr$`X`) # take the X chromosome - this one has everything
all_map_codes <- seq(1:length(founder_all_codes))
founder_all_lookup <- setNames(all_map_codes, founder_all_codes) # make maxmarg
founder_all_rev_lookup <- setNames(founder_all_codes, all_map_codes)


## location crossover comparison
pos_haploqa_cc <- locate_xo(cc_results[['ginf_haploqa']], cc_results[['map']])
pos_qtl2_cc <- cc_results[['pos']]

x_loc_diff <- location_xo_comp(pos_qtl2_cc, pos_haploqa_cc, num_chr)

## error matrix
## for each sample, plot each chromosome
## one function for each task
## 

sample_namelist <- rownames(cc_results[['pr']]$`1`)
err_comp(cc_results[['pr']], cc_results[['map']], num_chr)

## genocode comparison matrix
comp_matrix <- genocode_comp_matrix(cc_results[['map']], cc_results[['ginf']], cc_results[['ginf_haploqa']], founder_all_rev_lookup, num_chr, file_gen = F)


### percent genomic difference
chr_pct <- get_geno_pct_diff(cc_results[['ginf']], cc_results[['ginf_haploqa']], cc_results[['summary']], num_chr)

for (chr in num_chr) {
  chr_pct <- c()
  ind_list <- colnames(df_chr_qtl2_comp[2:(ncol(df_chr_haploqa))])
  
  for (ind in seq(1:(total_ind))) {
    ind_index <- ind+1
    #ind <- 192
    print(ind)
    ind_1 <- df_chr_qtl2_comp[,ind_index]
    ind_1_haplo <- df_chr_haploqa_comp[,ind_index]
    diff <- data.frame('qtl2' = ind_1, 'haploqa' = ind_1_haplo, 'pos' = df_chr_qtl2_comp$pos)
    diff <- diff %>% arrange(pos) %>% filter(qtl2 != haploqa)
    #write.csv(diff, paste0(comp_dir, '/diff_ind_', ind, '_chr_', chr, '.csv'), row.names = F)
    
    # compare differences
    diff_pct <- +(!((df_chr_haploqa[1:(ncol(df_chr_haploqa)-1)] == df_chr_qtl2[1:(ncol(df_chr_qtl2)-1)]) * 1))
    ind_df <- as.data.frame(diff_pct[,ind])
    
    pct <- as.numeric(colSums(ind_df) / nrow(ind_df))
    chr_pct <- c(chr_pct, pct)
    
    print(paste0('percent genomic difference between qtl2 and haploqa haplotype results of ind ', ind, ' on chr ', chr, ': ', pct))
  }
  
  chr_diff <- data.frame('individuals' = ind_list, 'genome_pct_diff' = chr_pct)
  chr_diff$chr <- chr
  chr_diff$ind <- seq(1:(total_ind))
  
  het_df <- summary_df %>% select(ID, `% Het. Calls`) %>% rename(individuals = ID, het_pct = '% Het. Calls')
  plot_diff <- merge(chr_diff, het_df, by = 'individuals')
  plot_diff$het_pct <- as.numeric(gsub("%", "", plot_diff$het_pct)) / 100
  # histogram 
  ggplot(plot_diff) + aes(x = genome_pct_diff, y = het_pct) + geom_point()
  ggplot(plot_diff, aes(x=genome_pct_diff)) + geom_histogram(binwidth = .01)
  # percent het (from haploqa)
  #write.csv(chr_diff, paste0(results_dir, '/pct_diff_chr_', chr, '.csv'), row.names = F)
  
}



### save the output of pr (calc_genoprob) in csv files?
#for (chr in num_chr) {
#  print(chr)
#  x <- eval(parse(text = deparse(substitute(chr))))
#  write.csv(pr[[chr]], paste0(results_dir, '/prob_chr_', x, '.csv'), row.names = F)
#}

