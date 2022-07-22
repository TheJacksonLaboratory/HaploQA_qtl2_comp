##### HaploQA data reformatting - Launch script
### used for qtl2 input
### annotation config saved in annotations_config.csv
### founder strain dictionary saved in founder_strains_table.csv

library(data.table)
library(rstudioapi)
library(qtl2)
library(dplyr)
library(tidyverse)
library(httr)

# source the function script
root <- dirname(getSourceEditorContext()$path)
source(paste0(root,"/input_data_prep_functions.R"))

#### change the below params
dir_name <- 'haploqa_cc_gm' # enter directory name
sample_type <- 'Collaborative Cross' # Collaborative Cross/MiniMUGA/GigaMUGA
output_dir_name <- 'cc_qtl2_gm_test'
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
results_dir <- file.path(root, 'results_test')
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
sample_url <- 'http://haploqa-dev.jax.org/tag/Collaborative%20Cross.html' # change to the URL of sample

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
    write.csv(sum_df, paste0(data_dir, '/collaborative_cross_summary.csv'), row.names = FALSE)
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
summary_df <- fread(paste0(data_dir, '/collaborative_cross_summary.csv'))

### Part 2 - convert results to qtl2 input format

# implement function
file_output <- get_qtl2_input(data_dir, sample_type, qtl2_dir, summary_df, list_pheno)

# order the chromosomes?
chr_order <- c((0:19),"X","Y")

### TO-DO: filter out chromosome 0/Y/M
### pull all data out of github

# unpack
df_geno <- file_output[[1]]
df_gmap <- file_output[[2]]
df_pmap <- file_output[[3]]
df_pheno <- file_output[[4]]
df_covar <- file_output[[5]]
df_foundergeno <- file_output[[6]]
df_crossinfo <- file_output[[7]]

if (qtl2_file_gen == TRUE) {
  write.csv(df_geno, paste0(qtl2_dir, '/test_geno.csv'), row.names = F)
  write.csv(df_gmap, paste0(qtl2_dir, '/test_gmap.csv'), row.names = F)
  write.csv(df_pmap, paste0(qtl2_dir, '/test_pmap.csv'), row.names = F)
  write.csv(df_pheno, paste0(qtl2_dir, '/test_pheno.csv'), row.names = F)
  write.csv(df_covar, paste0(qtl2_dir, '/test_covar.csv'), row.names = F)
  write.csv(df_foundergeno, paste0(qtl2_dir, '/test_foundergeno.csv'), row.names = F)
  write.csv(df_crossinfo, paste0(qtl2_dir, '/test_crossinfo.csv'), row.names = F)
}


control_file <- '{
  "description": "HaploQA data - qtl2 test run",
  "crosstype": "risib8",
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

#### Part 3 - qtl2
cc_test <- read_cross2(control_fp)
map <- cc_test$gmap 

# calculate probs
pr <- calc_genoprob(cross=cc_test, map=map, error_prob=0.002)

# output plots to pdf
total_ind <- 195
num_chr <- c((1:19),"X")


pdf(file = paste0(results_dir, '/test_plots.pdf'))
for (ind in seq(1:total_ind)) {
  print(paste0('individual: ',ind))
  for (chr in num_chr) {
    print(chr)
    plot_genoprob(pr, cc_test$gmap, ind = ind, chr = chr, main = paste0('Genotype Probabilities of individual ', ind, ' at chromosome ', chr))
    
  }
}
dev.off()


ind <- 1
chr = 2
plot_genoprob(pr, cc_test$gmap, ind = ind, chr = chr, main = paste0('Genotype Probabilities of individual ', ind, ' at chromosome ', chr))
# check homozygosity

for (chr in num_chr) {
  print(chr)
  x <- eval(parse(text = deparse(substitute(chr))))
  write.csv(pr[[chr]], paste0(results_dir, '/prob_chr_', x, '.csv'), row.names = F)
}

### combine guess_phase with plot_onegeno

