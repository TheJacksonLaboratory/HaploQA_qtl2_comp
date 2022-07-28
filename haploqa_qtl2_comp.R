library(rstudioapi)
library(qtl2)
library(dplyr)
library(data.table)

# source the function script
root <- dirname(getSourceEditorContext()$path)
source(paste0(root,"/input_data_prep_functions.R"))

## qtl2 directory
qtl2_dir <- file.path(root, output_dir_name) # name of desired output folder
# create if folder not exist
dir.create(qtl2_dir, showWarnings = FALSE)

## results directory
results_dir <- file.path(root, 'results')
dir.create(results_dir, showWarnings = FALSE) 

data_dir <- file.path(root, dir_name)

control_fp <- paste0(qtl2_dir, '/test.json')
summary_df <- fread(paste0(data_dir, '/collaborative_cross_summary.csv'))

### qtl2
cc_test <- read_cross2(control_fp)
map <- cc_test$gmap 

# attributes
total_ind <- 192
num_chr <- c((1:19),"X")


# calculations
pr <- calc_genoprob(cross=cc_test, map=map, error_prob=0.002)
ginf <- maxmarg(pr, minprob = 0.01) # lower minprob to 0.01
ph_geno <- guess_phase(cc_test, ginf)
pos <- locate_xo(ginf, map) # locations of crossovers in each individual on each chr

founder_codes <- c('AA', 'BB', 'CC', 'DD', 'EE', 'FF', 'GG', 'HH')
map_codes <- c(1, 2, 3, 4, 5, 6, 7, 8)
founder_lookup <- setNames(founder_codes, map_codes)
founder_reverse <- setNames(map_codes, founder_codes)
founders_dict <- fread(paste0(root, '/founder_strains_table.csv'))
founder_haplo_lookup <- setNames(founders_dict$letter, founders_dict$founder_strain)

## process haploqa data
df_all <- get_haplotypes(summary_df, data_dir)
df_test <- df_all %>% select(sample_id, snp_id, haplotype1, haplotype2, chromosome) %>% unique()
df_test[,c(3,4)] <- as.data.frame(apply(df_test[,c(3,4)], 2, function(x) founder_haplo_lookup[x]))
df_test$haplotype <- paste(df_test$haplotype1, df_test$haplotype2, sep='')


# map ginf(output of maxmarg) to founder strains in the order of AA, BB, CC, etc.
# then compare with haplotypes from HaploQA
# run for each chromosome
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
  
  qtl2_comp <- df_chr_qtl2_comp[1:(ncol(df_chr_haploqa))]
  haploqa_comp <- df_chr_haploqa_comp[1:(ncol(df_chr_haploqa_comp))]
  
  temp <- rbind(qtl2_comp, qtl2_comp)
  temp %>% arrange(marker)
  temp %>% group_by(marker) %>% summarise(Freq=n_distinct(`35V`))
  
  dcast(setDT(qtl2_comp), marker~everything(), length)
  dcast(setDT(haploqa_comp), marker~`35V`, length)
  
  table(qtl2_comp)
  
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
sample_list <- rownames(cc_test$geno$`1`)

df_haplo <- df_test %>% select(sample_id, snp_id, chromosome, haplotype) %>% rename(marker = snp_id, chr = chromosome)

df_haplo[,4] <- as.data.frame(apply(df_haplo[,4], 2, function(x) founder_reverse[x])) ## heterozygosity?
df_haplo <- df_haplo[!is.na(df_haplo$haplotype)] ## drop the heterozygous ones for now to get it up and running 
df_haplo <- df_haplo %>% filter(sample_id %in% sample_list)

df_chr_map <- df_gmap %>% filter(marker %in% unique(df_haplo$marker)) %>%
  select(marker, chr, pos)

df_chr_all <- merge(df_haplo, df_gmap, by = 'marker') # to align with the markers and samples
df_chr_all <- df_chr_all %>% select(marker, sample_id, pos, chr.x, haplotype) %>% rename(chr = chr.x)

### drop the cross_info samples?
df_haplo_geno <- list()
df_haplo_map <- list()
for (i in num_chr) {
  print(i)
  #i <- 1
  ## output should be the same as cc_test$geno
  df_chr <- df_chr_all %>% filter(chr == i) %>% select(-c(chr, pos))
  df_chr <- dcast(df_chr, sample_id~marker, value.var = 'haplotype')
  
  rownames(df_chr) <- df_chr$sample_id
  df_chr <- df_chr %>% select(-c(sample_id)) %>% as.matrix()
  
  df_haplo_geno[[i]] <- df_chr
  
  ## output should be the same as cc_test$gmap
  df_map <- df_chr_all %>% filter(chr == i) %>% select(marker, pos) %>% unique() %>% as.data.frame()
  rownames(df_map) <- df_map$marker
  df_map <- df_map %>% select(-marker) %>% t() %>% as.matrix()
  df_haplo_map[[i]] <- df_map
}

df_ch <- df_chr_all %>% filter(chr == 2) %>% select(-c(chr, pos))
temp <- t(df_haplo_geno[[1]]) 
temp$marker <- rownames(temp)
temp_m <- t(df_haplo_map[[1]])
temp_m$marker <- rownames(temp_m)
temp_1 <- merge(temp, temp_m, by = 'marker')




temp$
  t(tem)
t(temp_m)


temp1 <- df_haplo_geno[[2]]

temp1_m <- df_haplo_map[[2]]
ind <- 2
plot_onegeno(df_haplo_geno, df_haplo_map, ind = ind, shift = TRUE, main = paste0('Geno-wide genotypes of individual ', ind))  # one individual's geno_wide genotype

### visualizations
### to-do: structures of files? and write into a function
# genowide genotype
pdf(file = paste0(results_dir, '/test_geno_plots_comp.pdf'), width=20, height=10, bg="white")
for (ind in seq(1:total_ind)) {
  print(ind)
  #ind <- 1
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
library(qtl)
data(map10)
chrL <- round(summary(map10)$length[1:20])
map <- sim.map(len=chrL, n.mar=chrL+1, include.x=TRUE, eq.spacing=TRUE)
n <- 500
bc <- sim.cross(map, n.ind=n, type="bc", missing.prob=0.05)
bc$pheno$sex <- rep(0:1, each=n/2)







