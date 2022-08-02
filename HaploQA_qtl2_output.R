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


# source the function script
root <- dirname(getSourceEditorContext()$path)
source(paste0(root,"/input_data_prep_functions.R"))

#### change the below params
dir_name <- 'haploqa_do' # enter directory name
sample_type <- 'DO' # Collaborative Cross/MiniMUGA/GigaMUGA/DO
output_dir_name <- 'do_qtl2_genail'
list_pheno <- c('WBC', 'NEUT')

#### Set toggles
generate_summary <- TRUE # set to TRUE to generate summary table
output_summary <- TRUE # set to TRUE to output summary table to data directory
generate_sample_individuals <- TRUE # set to TRUE to generate individual sample data
output_samples <- TRUE 
qtl2_file_gen <- TRUE


### Environment
# data output directory
data_dir <- file.path(root, dir_name)
dir.create(data_dir, showWarnings = FALSE) 

## create a data directory for qtl2 input data
qtl2_dir <- file.path(root, output_dir_name) # name of desired output folder
# create if folder not exist
dir.create(qtl2_dir, showWarnings = FALSE)

## results directory
results_dir <- file.path(root, 'results')
dir.create(results_dir, showWarnings = FALSE) 

### control file
control_fp <- paste0(qtl2_dir, '/test.json')
if (file.exists(control_fp) == FALSE) {
  file.create(control_fp)
}


###############################################################################
### Part 1 - get results from HaploQA
# set main URL domain
url_domain <- 'http://haploqa-dev.jax.org/' 
# target sample URL
sample_url <- 'http://haploqa-dev.jax.org/tag/J:DO.html' # change to the URL of sample

#### html file prep
# read html file
haploqa_html <- read_html(sample_url)

# extract information from html file
html_temp <- haploqa_html %>% html_nodes("a") %>% html_attr("href")
url_list <- paste0(url_domain, html_temp[grepl('/sample', html_temp)])

#### implementations
# summary table
if (generate_summary == TRUE) {
  print(paste0('Working on summary file:'))
  sum_df <- sample_summary_scrape(haploqa_html, url_list)
  if (output_summary == TRUE) {
    print('Writing to directory')
    write.csv(sum_df, paste0(data_dir, '/', sample_type, '_summary.csv'), row.names = FALSE)
  }
}


### TO-DO: URLs with recurring issues

#skip_urls <- c('http://haploqa-dev.jax.org//sample/5f56294d183b6a06dabbe5ee.html', 
#               'http://haploqa-dev.jax.org//sample/5f562977183b6a06dabbe60c.html',
#               'http://haploqa-dev.jax.org//sample/5f562991183b6a06dabbe61f.html',
#               'http://haploqa-dev.jax.org//sample/5696759d85b062edc4efbe7d.html',
#               'http://haploqa-dev.jax.org//sample/569675b385b062edc4efbe8b.html',
#               'http://haploqa-dev.jax.org//sample/569675b585b062edc4efbe8d.html',
#               'http://haploqa-dev.jax.org//sample/569675e485b062edc4efbeab.html',
#               'http://haploqa-dev.jax.org//sample/569675e785b062edc4efbead.html',
#               'http://haploqa-dev.jax.org//sample/569675ee85b062edc4efbeb2.html',
#               'http://haploqa-dev.jax.org//sample/5696761285b062edc4efbecc.html')


# list of urls to generate samples from
url_ind <- unique(sum_df$`Sample Filepath`)
#url_ind <- url_ind[!url_ind %in% skip_urls]

# individual samples
inc = 0
for (url in url_ind) {
  inc = inc + 1 # increment
  if (generate_sample_individuals == TRUE) {
    file <- sample_individual_scrape(url, url_domain)
    print(paste0('Working on file ', inc, '/', length(url_ind), ': ', file))
    sample_df_save <- as.data.frame(content(GET(file)))
    if (output_samples == TRUE) {
      print('Writing to directory')
      file_name <- unlist(strsplit(file, '/'))[6]
      GET(file, write_disk(paste0(data_dir, '/', file_name), overwrite = TRUE), show_col_types = FALSE)
    }
  }
}

# for testing - read the summary file if not regenerating above, as it's needed for below
summary_df <- fread(paste0(data_dir, '/DO_summary.csv'))

### Part 2 - convert results to qtl2 input format

# implement function
if(qtl2_file_gen == TRUE) {
  file_output <- get_qtl2_input(data_dir, sample_type, qtl2_dir, summary_df, list_pheno)
}


# order the chromosomes?
chr_order <- c((0:19),"X","Y")

### TO-DO: filter out chromosome 0/Y/M
### pull all data out of github

df_geno <- read.csv(paste0(qtl2_dir, '/test_geno.csv'))
df_gmap <- read.csv(paste0(qtl2_dir, '/test_gmap.csv'))

# unpack
df_geno <- file_output[[1]]
df_gmap <- file_output[[2]]
df_pmap <- file_output[[3]]
df_pheno <- file_output[[4]]
df_covar <- file_output[[5]]
df_foundergeno <- file_output[[6]] # with strain id metadata
df_crossinfo <- file_output[[7]]
#df_crossinfo <- df_crossinfo %>% select(-strain_id) # read_cross only takes 8 columns

if (qtl2_file_gen == TRUE) {
  write.csv(df_geno, paste0(qtl2_dir, '/test_geno.csv'), row.names = F)
  write.csv(df_gmap, paste0(qtl2_dir, '/test_gmap.csv'), row.names = F)
  write.csv(df_pmap, paste0(qtl2_dir, '/test_pmap.csv'), row.names = F)
  write.csv(df_pheno, paste0(qtl2_dir, '/test_pheno.csv'), row.names = F)
  write.csv(df_covar, paste0(qtl2_dir, '/test_covar.csv'), row.names = F)
  write.csv(df_foundergeno, paste0(qtl2_dir, '/test_foundergeno.csv'), row.names = F)
  write.csv(df_crossinfo, paste0(qtl2_dir, '/test_crossinfo.csv'), row.names = F)
}

## genail8
# id ngen a,b,c,d,e,f,g,h
control_file <- '{
  "description": "HaploQA data - qtl2 test run",
  "crosstype": "genail8",
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



# output to file
if (qtl2_file_gen == TRUE) {
  writeLines(control_file, control_fp)
}

## dont do the read_cross
# wrangle haploqa raw data into such format, then imitate the data structure so that output is exactly the same as guess_phase
# plot_compare_geno

#control_fp <- '/Users/linc/Documents/GitHub/HaploQA_qtl2_comp/cc_qtl2_genail/test.json'
#df_gmap <- fread('/Users/linc/Documents/GitHub/HaploQA_qtl2_comp/cc_qtl2_genail/test_gmap.csv')
#### Part 3 - qtl2
# read data
cc_test <- read_cross2(control_fp)
map <- cc_test$gmap 

# attributes
total_ind <- 277


num_chr <- c((1:19),"X")

# calculations
pr <- calc_genoprob(cross=cc_test, map=map, error_prob=0.002)
ginf <- maxmarg(pr, minprob = 0.01) # lower minprob to 0.01
ph_geno <- guess_phase(cc_test, ginf)
pos <- locate_xo(ginf, map) # locations of crossovers in each individual on each chr

founder_codes <- c('AA', 'BB', 'CC', 'DD', 'EE', 'FF', 'GG', 'HH')
single_codes <- sapply(founder_codes, function(x) substr(x, 1, 1))
map_codes <- c(1, 2, 3, 4, 5, 6, 7, 8)
founder_lookup <- setNames(founder_codes, map_codes)
founder_reverse <- setNames(map_codes, single_codes)
founders_dict <- fread(paste0(root, '/founder_strains_table.csv'))
founder_haplo_lookup <- setNames(founders_dict$letter, founders_dict$founder_strain)

summary_df <- fread('/Users/linc/Documents/GitHub/HaploQA_qtl2_comp/haploqa_cc_gm/collaborative_cross_summary.csv')
## process haploqa data
df_all <- get_haplotypes(summary_df, data_dir)
df_test <- df_all %>% select(sample_id, snp_id, haplotype1, haplotype2, chromosome) %>% unique()
df_test[,c(3,4)] <- as.data.frame(apply(df_test[,c(3,4)], 2, function(x) founder_haplo_lookup[x]))
df_test$haplotype <- paste(df_test$haplotype1, df_test$haplotype2, sep='')

# map ginf(output of maxmarg) to founder strains in the order of AA, BB, CC, etc.
# then compare with haplotypes from HaploQA
# run for each chromosome
comp_matrix <- list()
for (chr in num_chr) {
  ## qtl2
  chr <- 1
  map_chr <- as.data.frame(map[[chr]]) %>% rename('pos' = 'map[[chr]]')
  map_chr$marker <- rownames(map_chr)
  # rowname is sample id
  df <- as.data.frame(ginf[[chr]])
  df_chr_qtl2 <- apply(df, 2, function(x) founder_lookup[x]) # column level
  rownames(df_chr_qtl2) <- rownames(df)
  
  ## haploqa
  df_chr_haploqa <- df_test %>% filter(chromosome == chr) %>% select(sample_id, snp_id, haplotype) %>% dcast(sample_id~snp_id, value.var = 'haplotype')
  rownames(df_chr_haploqa) <- df_chr_haploqa$sample_id
  df_chr_haploqa <- df_chr_haploqa %>% select(-sample_id)
  
  ## this haplotype info was pulled from unprocessed df, so some snps/samples might not be in qtl2 (due to genotype discrepancies, etc)
  ## make sure all snps are consistent with that of genotype data
  df_chr_haploqa <- df_chr_haploqa[rownames(df_chr_qtl2),colnames(df_chr_qtl2)]
  df_chr_haploqa <- as.data.frame(t(df_chr_haploqa))
  df_chr_haploqa$marker <- rownames(df_chr_haploqa)
  df_chr_haploqa_comp <- merge(df_chr_haploqa, map_chr, by = 'marker')
  
  df_chr_qtl2 <- as.data.frame(t(df_chr_qtl2))
  df_chr_qtl2$marker <- rownames(df_chr_qtl2)
  df_chr_qtl2_comp <- merge(df_chr_qtl2, map_chr, by = 'marker')
  
  qtl2_comp <- df_chr_qtl2_comp[1:(ncol(df_chr_haploqa))] %>% select(-c(marker))
  haploqa_comp <- df_chr_haploqa_comp[1:(ncol(df_chr_haploqa_comp))] %>% select(-c(marker, pos))
  
  col_list <- names(haploqa_comp)
  
  ### loop through each matching column within both dataframes
  result_df <- data.frame()
  for (col in seq(1:ncol(qtl2_comp))) {
    #col <- 1
    print(col)
    df_match <- data.frame('qtl2_ind' = qtl2_comp[,col], 'haploqa_ind' = haploqa_comp[,col])
    t1 <- dcast(df_match, qtl2_ind ~ haploqa_ind, value.var = 'haploqa_ind', fun.aggregate = length)
    result_df <- bind_rows(result_df, t1) %>% group_by(qtl2_ind) %>% summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>% as.data.frame()
  }
  result_df <- result_df %>% rename(geno_code = qtl2_ind) 
  
  result_df <- result_df%>%
    select(geno_code, AA, BB, CC, DD, EE, FF, GG, HH, sort(names(.)))
  
  comp_matrix[[chr]] <- result_df
  #write.csv(result_df, paste0(results_dir, '/comp_matrix_chr_', chr, '.csv'), row.names = F)
  
  ### percent genomic difference
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



# keep the samples that were genotyped for better comparison

df_haplo <- df_test %>% select(sample_id, snp_id, chromosome, haplotype1, haplotype2) %>% rename(marker = snp_id, chr = chromosome)
df_haplo[,c(4,5)] <- as.data.frame(apply(df_haplo[,c(4,5)], 2, function(x) founder_reverse[x])) ## heterozygosity?

df_chr_all <- merge(df_haplo, df_gmap, by = 'marker') # to align with the markers and samples
df_chr_all <- df_chr_all %>% select(marker, sample_id, pos, chr.x, haplotype1, haplotype2) %>% rename(chr = chr.x, mom = haplotype1, dad = haplotype2)

### drop the cross_info samples?

#qtl2_geno <- cc_test$geno

df_haplo_geno <- list()
df_haplo_map <- list()
for (i in num_chr) {
  #i<-1
  print(i)
  #marker_list <- colnames(qtl2_geno[[i]])
  #sample_list <- rownames(qtl2_geno[[i]])
  #i <- 'X'
  ## output should be the same as cc_test$gmap
  df_map <- df_chr_all %>% filter(chr == i) %>% select(marker, pos) %>% unique() %>% as.data.frame() %>% arrange(pos)
  rownames(df_map) <- df_map$marker
  df_map_1 <- df_map %>% select(-marker) %>% t() #%>% as.vector()
  #names(df_map_1) <- df_map$marker
  df_haplo_map[[i]] <- df_map_1
  
  col_order <- df_map$marker
  
  ## output should be the same as cc_test$geno
  df_chr <- df_chr_all %>% filter(chr == i) %>% select(-c(chr, pos))
  #df_chr <- df_chr %>% filter(marker %in% marker_list) %>% filter(sample_id %in% sample_list)
  df_mom <- df_chr %>% select(marker, sample_id, mom) %>% dcast(sample_id~marker, value.var = 'mom')
  rownames(df_mom) <- df_mom$sample_id
  df_mom <- df_mom %>% select(-c(sample_id)) %>% as.matrix()
  df_mom <- df_mom[,col_order]
  df_dad <- df_chr %>% select(marker, sample_id, dad) %>% dcast(sample_id~marker, value.var = 'dad')
  rownames(df_dad) <- df_dad$sample_id
  df_dad <- df_dad %>% select(-c(sample_id)) %>% as.matrix()
  df_dad <- df_dad[,col_order]
  
  df_haplo_geno[[i]] <- array(dim = c(nrow(df_mom), ncol(df_mom),2), dimnames=list(rownames(df_mom), colnames(df_mom), c('mom', 'dad'))) # should be the same
  df_haplo_geno[[i]][,,1] <- df_dad
  df_haplo_geno[[i]][,,2] <- df_mom
  
  
  
  
}



### visualizations
### to-do: structures of files? and write into a function
# genowide genotype
pdf(file = paste0(results_dir, '/geno_plots_comp_cc.pdf'), width=20, height=10, bg="white")
for (ind in seq(1:total_ind)) {
  print(ind)
  #ind <- #
  par(mfrow = c(1, 2))
  plot_onegeno(ph_geno, map, ind = ind, shift = TRUE, main = paste0('qtl2 - Geno-wide genotypes of individual ', ind))  # one individual's geno_wide genotype
  plot_onegeno(df_haplo_geno, df_haplo_map, ind = ind, shift = TRUE, main = paste0('haploqa - Geno-wide genotypes of individual ', ind))  # one individual's geno_wide genotype
}
dev.off()

# genoprob
pdf(file = paste0(results_dir, '/test_plots.pdf'))
for (ind in seq(1:total_ind)) {
  print(paste0('individual: ',ind))
  for (chr in num_chr) {
    print(chr)
    plot_genoprob(pr, cc_test$gmap, ind = ind, chr = chr, main = paste0('Genotype Probabilities of individual ', ind, ' at chromosome ', chr))
    
  }
}
dev.off()

### csv files (?)
# output 
for (chr in num_chr) {
  print(chr)
  x <- eval(parse(text = deparse(substitute(chr))))
  write.csv(pr[[chr]], paste0(results_dir, '/prob_chr_', x, '.csv'), row.names = F)
}

### markdown file that goes through the entire pipeline for gigamuga, minimuga, or other chosen samples


