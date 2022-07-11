library(data.table)
library(rstudioapi)
library(qtl2)
library(dplyr)
library(tidyverse)
library(reshape2)

### get all data from directory
dir_name <- 'haploqa_MiniMUGA' # enter directory name
root <- dirname(getSourceEditorContext()$path)
# direct to the data directory filepath
data_dir <- file.path(root, dir_name)

# summary file
### to-do: split ID and secondary ID into two columns
sum_df <- fread(paste0(data_dir, '/collaborative_cross_summary.csv'))
#sum_df$file_name <- filename of scraped data

# all txt file from directory
data_files <- dir(data_dir, pattern = '\\.txt$', full.names = TRUE)

#### qtl2 input file format
# https://kbroman.org/qtl2/assets/vignettes/input_files.html#Reading_the_data_files
# 1. genotype data - marker/snp x individual samples
# 2. phenotype file
# 3. covariate file
# 4. genetic map file
# 5. cross info file
# 6. control file

### combine all files
cc_all <- rbindlist(lapply(data_files, fread)) # use sample_id or original?

### genotype data
cc_geno_sub <- cc_all %>% select(sample_id, snp_id, allele1_fwd, allele2_fwd) %>%
  mutate(gene_exp = paste(allele1_fwd, allele2_fwd, sep = '')) %>%
  select(sample_id, snp_id, gene_exp)
cc_geno <- dcast(cc_geno_sub, snp_id~sample_id, value.var="gene_exp") # inbred - most are homozygous

### genetic map data (pos in cM)
# 1cM = 1,000,000bp
cc_gmap <- cc_all %>% select(snp_id, chromosome, position_bp) %>%
  mutate(pos = position_bp/1000000) %>%
  select(snp_id, chromosome, pos)

### physical map data (pos in bp? mbp?)
# 1 Mbp =  1cM?? hmm...
#cc_gmap <- cc_all %>% select(snp_id, chromosome, position_bp) %>%
  #mutate(pos = position_bp/1000000) %>%
  #select(snp_id, chromosome, pos)

### annotations?
dir <- "https://raw.githubusercontent.com/kbroman/MUGAarrays/main/UWisc/"
mm <- read.csv(paste0(dir, "mm_uwisc_v1.csv"))
gm <- read.csv(paste0(dir, "gm_uwisc_v1.csv"))
mm_dict <- read.csv(paste0(dir, "mm_uwisc_dict_v1.csv"))
gm_dict <- read.csv(paste0(dir, "gm_uwisc_dict_v1.csv"))
common <- read.csv(paste0(dir, "mm_gm_commonmark_uwisc_v1.csv"))
common_dict <- read.csv(paste0(dir, "mm_gm_commonmark_uwisc_dict_v1.csv"))

