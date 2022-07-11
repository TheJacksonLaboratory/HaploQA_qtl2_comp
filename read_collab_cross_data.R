library(data.table)
library(rstudioapi)
library(qtl2)
library(dplyr)
library(tidyverse)
library(reshape2)

## params
dir_name <- 'haploqa_collab_cross' # enter directory name
annot_file <- 'gm_uwisc_v1.csv'
dict_file <- 'gm_uwisc_dict_v1.csv'

### get all data from directory
root <- dirname(getSourceEditorContext()$path)
# direct to the data directory filepath
data_dir <- file.path(root, dir_name)

# summary file
### to-do: split ID and secondary ID into two columns
#sum_df <- fread(paste0(data_dir, '/collaborative_cross_summary.csv'))
#sum_df$file_name <- filename of scraped data

# all txt file from directory
data_files <- dir(data_dir, pattern = '\\.txt$', full.names = TRUE)

#### qtl2 input file format
# https://kbroman.org/qtl2/assets/vignettes/input_files.html#Reading_the_data_files
# 1. genotypes (done)
# 2. phenotypes
# 3. phenotype covariates (i.e. tissue type, time points)
# 4. genetic map (done)
# 5. physical map (optional) (done?)
# 6. control file (YAML or JSON format, not CSV) (later)

### combine all files
df_all <- rbindlist(lapply(data_files, fread)) # use sample_id or original?

### genotype data
geno_sub <- df_all %>% select(sample_id, snp_id, allele1_fwd, allele2_fwd) %>%
  mutate(gene_exp = paste(allele1_fwd, allele2_fwd, sep = '')) %>%
  select(sample_id, snp_id, gene_exp)
df_geno <- dcast(geno_sub, snp_id~sample_id, value.var="gene_exp") # inbred - most are homozygous

### genetic map data (pos in cM)
### annotations
dir <- "https://raw.githubusercontent.com/kbroman/MUGAarrays/main/UWisc/"

annot_df <- read.csv(paste0(dir, annot_file))
annot_dict <- read.csv(paste0(dir, dict_file))
# 1cM = 1,000,000bp
map_df <- df_all %>% 
  select(snp_id, chromosome, position_bp) %>%
  rename(marker = snp_id, chr = chromosome) %>%
  unique()

annot_map <- annot_df %>% select(marker, chr, bp_mm10, cM_cox) %>% unique()

map_all <- merge(map_df, annot_map, by = c('marker'))

df_gmap <- map_all %>% select(marker, chr.x, cM_cox) %>% rename(chr = chr.x, pos = cM_cox)

df_pmap <- map_all %>% select(marker, chr.x, bp_mm10) %>% rename(chr = chr.x, pos = bp_mm10)


