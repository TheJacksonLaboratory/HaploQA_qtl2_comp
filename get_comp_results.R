## contains all functions to calculate accuracy metrics and/or to compare results 
### function directory
# 1. location_xo_comp - calculates the distances between crossover locations of haploqa and qtl2 results
# 2. err_comp
# 3. comp_df_int

library(rstudioapi)
root <- dirname(getSourceEditorContext()$path)

# calculates the distance between crossover locations for the two models
# 1. compares the location crossovers values (which are usually lists) between haploqa and qtl2
# 2. determine which one has more crossovers than the other (length of lists)
# 3. maps the values in the shorter crossover list each to the closest value in the longer list
# 4. finds the absolute value between the pair
#
# @param qtl2_pos (listed matrices) - output of locate_xo() for qtl2
# @param haploqa_pos (listed matrices) - output of locate_xo() for haploqa
# @param num_chr (list) - list of chromosomes
#
# @output x_loc_diff (listed dataframe) - results of the calculations, with the results for each chromosome being an element of the list
location_xo_comp <- function(qtl2_pos, haploqa_pos, num_chr) {
  x_loc_diff <- list()
  for (chr in num_chr) {
    print(chr)
    pos_h_chr <- haploqa_pos[[chr]]
    pos_chr <- qtl2_pos[[chr]]
    nmar_h <- sapply(pos_h_chr, length)
    nmar <- sapply(pos_chr, length) # none has the same length...
    #x_loc_diff[[chr]] <- list()
    for (i in seq(1:length(pos_chr))) {
      #i <- 149
      x <- unlist(pos_h_chr[i])
      y <- unlist(pos_chr[i])
      if (length(x) > length(y)) {
        df <- apply(abs(outer(x, y, '-')), 2, min)
      }
      if (length(x) < length(y)) {
        df <- apply(abs(outer(y, x, '-')), 2, min)
      }
      x_loc_diff[[chr]][[i]] <- ifelse(is.na(mean(df)), 0, mean(df))
    } # make a histogram of some sort, 100 bins
  }
  return(x_loc_diff)
}


err_comp <- function(pr, map, num_chr, results_dir, sample_type) {
  pdf(paste0(results_dir, 'err_comp_plots_', sample_type, '.pdf'))
  facet_all <- data.frame()
  for (chr in num_chr) {
    #chr <- 2
    print(chr)
    err <- (1 - apply(pr[[chr]], c(1,3), max)) %>% t() %>% as.data.frame()
    err$marker <- rownames(err)
    err_map <- map[[chr]] %>% as.data.frame() %>% rename(pos = '.')
    err_map$marker <- rownames(err_map)
    
    err_chr <- merge(err, err_map, by = 'marker')
    
    for (col in seq(2,(ncol(err_chr)-2))) {
      print(col)
      #col <- 2
      df <- err_chr[,c(col, ncol(err_chr))]
      df$sample <- colnames(df)[1]
      colnames(df)[1] <- 'value'
      df$chromosome <- chr
      #plot <- ggplot(df) + aes(x = df[,2], y = df[,1]) + geom_line() + xlab('pos') + ylab(names(df)[1])
      #print(plot)
      facet_all <- rbind(facet_all, df)
    }
  }
  
  for(sample in unique(facet_all$sample)) {
    #sample <- '35V'
    sample_df <- facet_all[facet_all$sample == sample,]
    plot <- ggplot(sample_df) + aes(x = pos, y = value) + geom_line() + facet_wrap(sample~chromosome)
    print(plot)
  }
  dev.off()
  
  #return(df)
}


### simulating maxmarg for some reason?? there's already a *better* function built for it
comp_df_int <- function(chr, ginf, founder_lookup_table, map_chr) {
  df <- as.data.frame(ginf[[chr]])
  df_chr <- apply(df, 2, function(x) founder_lookup_table[x]) # column level
  rownames(df_chr) <- rownames(df)
  df_chr <- as.data.frame(t(df_chr))
  df_chr$marker <- rownames(df_chr)
  df_chr_comp <- merge(df_chr, map_chr, by = 'marker')
  
  df_comp <- df_chr_comp[1:(ncol(df_chr_comp))] %>% select(-c(marker))
  
  return(df_comp)
}

# do this for each strain
genocode_comp_matrix <- function(map, ginf_qtl2, ginf_haploqa, founder_lookup_table, num_chr, file_gen, sample_type) {
  comp_matrix <- list()
  for (chr in num_chr) {
    ## qtl2
    print(chr)
    map_chr <- as.data.frame(map[[chr]]) %>% rename('pos' = 'map[[chr]]')
    map_chr$marker <- rownames(map_chr)
    
    haploqa_comp <- comp_df_int(chr, ginf_haploqa, founder_lookup_table, map_chr)
    qtl2_comp <- comp_df_int(chr, ginf_qtl2, founder_lookup_table, map_chr)
    
    ### loop through each matching column within both dataframes
    result_df <- data.frame()
    for (col in seq(1:(ncol(qtl2_comp)-1))) {
      #print(col)
      df_match <- data.frame('qtl2_ind' = qtl2_comp[,col], 'haploqa_ind' = haploqa_comp[,col])
      t1 <- dcast(df_match, qtl2_ind ~ haploqa_ind, value.var = 'haploqa_ind', fun.aggregate = length)
      result_df <- bind_rows(result_df, t1) %>% group_by(qtl2_ind) %>% summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>% as.data.frame()
    }
    
    result_df <- result_df %>% rename(geno_code = qtl2_ind) %>%
      select(geno_code, AA, BB, CC, DD, EE, FF, GG, HH, sort(names(.)))
    
    comp_matrix[[chr]] <- result_df
    
    if(file_gen == T) {
      comp_fp <- paste0(root, '/geno_comp_results_', sample_type)
      dir.create(comp_fp, showWarnings = FALSE) 
      write.csv(result_df, paste0(comp_fp, '/geno_comp_matrix_', chr, '.csv'), row.names = F)
    }
  }
  return(comp_matrix)
}


get_geno_pct_diff <- function(ginf_haploqa, ginf_qtl2, summary_df, num_chr, founder_lookup_table) {
  
  haploqa_comp <- comp_df_int(ginf_haploqa, founder_lookup_table, map_chr)
  qtl2_comp <- comp_df_int(ginf_qtl2, founder_lookup_table, map_chr)
  
  for (chr in num_chr) {
    chr_pct <- c()
    print(chr)
    for (ind in seq(1:(ncol(qtl2_comp)-1))) {
      #ind <- 1
      #ind_index <- ind+1
      ind_qtl2 <- qtl2_comp[,ind]
      ind_haplo <- haploqa_comp[,ind]
      diff <- data.frame('qtl2' = ind_qtl2, 'haploqa' = ind_haplo, 'pos' = qtl2_comp$pos)
      diff <- diff %>% arrange(pos) %>% filter(qtl2 != haploqa)
      #write.csv(diff, paste0(comp_dir, '/diff_ind_', ind, '_chr_', chr, '.csv'), row.names = F)
      
      diff_pct <- +(!((haploqa_comp[1:(ncol(haploqa_comp)-1)] == qtl2_comp[1:(ncol(qtl2_comp)-1)]) * 1))
      ind_df <- as.data.frame(diff_pct[,ind])
      pct <- as.numeric(colSums(ind_df) / nrow(ind_df))
      chr_pct <- c(chr_pct, pct)
    }
    
    chr_diff <- data.frame('individuals' = colnames(qtl2_comp[1:(ncol(qtl2_comp)-1)]), 'genome_pct_diff' = chr_pct)
    chr_diff$chr <- chr
    chr_diff$ind <- seq(1:(ncol(qtl2_comp)-1))
    
    het_df <- summary_df %>% select(ID, `% Het. Calls`) %>% rename(individuals = ID, het_pct = '% Het. Calls')
    plot_diff <- merge(chr_diff, het_df, by = 'individuals')
    plot_diff$het_pct <- as.numeric(gsub("%", "", plot_diff$het_pct)) / 100
    # histogram 
    ggplot(plot_diff) + aes(x = genome_pct_diff, y = het_pct) + geom_point()
    ggplot(plot_diff, aes(x=genome_pct_diff)) + geom_histogram(binwidth = .01)
    
    #write.csv(chr_diff, paste0(results_dir, '/pct_diff_chr_', chr, '.csv'), row.names = F)
  }
  return(chr_pct)
}



xo_number_plot <- function(sample_result) {
  sample_result <- cc_results
  sample_type <- toupper(strsplit(deparse(substitute(sample_result)), '_')[[1]][1])
  pos_haploqa_do <- locate_xo(sample_result[['ginf_haploqa']], sample_result[['map']])
  pos_qtl2_do <- sample_result[['pos']]
  
  
  xo_all <- data.frame()
  for (i in seq(1, length(pos_haploqa_do[[1]]))) { # should be the same for all chromosomes
    print(i)
    #i <- 1
    xo_haplo <- list()
    xo_qtl2 <- list()
    for (chr in seq(1,19)) {
      #chr <- 1
      xo_haplo <- c(xo_haplo, length(pos_haploqa_do[[chr]][[i]]))
      xo_qtl2 <- c(xo_haplo, length(pos_qtl2_do[[chr]][[i]]))
    }
    xo_haplo_total <- sum(unlist(xo_haplo)) # number of crossovers for one sample across entire genome (all chromosomes)
    xo_qtl2_total <- sum(unlist(xo_qtl2)) 
    df <- data.frame(sample = names(pos_haploqa_do$`1`)[i], total_sample_xo_haplo = xo_haplo_total, total_sample_xo_qtl2 = xo_qtl2_total)
    xo_all <- rbind(xo_all, df)
  }
  
  xo_all$total_sample_xo_haplo <- as.numeric(xo_all$total_sample_xo_haplo)
  xo_all$total_sample_xo_qtl2 <- as.numeric(xo_all$total_sample_xo_qtl2)
  #sd_col <- sapply(xo_all[c(2, 3)], sd)
  high_lim <- max(xo_all[,c(2,3)])
  low_lim <- min(xo_all[,c(2,3)])
  
  plot <- ggplot(xo_all) + aes(x = total_sample_xo_haplo, y = total_sample_xo_qtl2) + geom_point() +
    xlab('total number of crossovers on autosomes - HaploQA') + ylab('total number of crossovers on autosomes - qtl2') + 
    xlim(low_lim, high_lim) + ylim(low_lim, high_lim) + ggtitle(paste0('Sample: ', sample_type, ' - each point represents one individual')) +
    geom_abline(intercept = 0, slope = 1, linetype="dashed")
  
  print(plot)
  
}


loc_xo_distance_plot <- function(sample_result) {
  #sample_result <- cc_results
  num_chr <- c((1:19),"X")
  sample_type <- toupper(strsplit(deparse(substitute(sample_result)), '_')[[1]][1])
  pos_haploqa_do <- locate_xo(sample_result[['ginf_haploqa']], sample_result[['map']])
  pos_qtl2_do <- sample_result[['pos']]
  x_loc_diff <- location_xo_comp(pos_qtl2_do, pos_haploqa_do, num_chr)
  
  
  loc_xo_dist_all <- data.frame()
  for (i in seq(1, length(pos_haploqa_do[[1]]))) { # should be the same for all chromosomes
    print(i)
    #i <- 1
    xo_dist_sample <- list()
    for (chr in seq(1,19)) {
      #chr <- 1
      xo_dist_sample <- c(xo_dist_sample, x_loc_diff[[chr]][[i]])
    }
    df <- data.frame(sample = names(pos_haploqa_do$`1`)[i], loc_xo_dist = mean(unlist(xo_dist_sample)))
    loc_xo_dist_all <- rbind(loc_xo_dist_all, df)
  }
  
  st_err <- sd(loc_xo_dist_all$loc_xo_dist)
  
  plot <- ggplot(loc_xo_dist_all) + aes(x = seq(1, length(sample)), y = loc_xo_dist) + geom_point() + 
    xlab('Sample index') + ylab('mean of crossover distance') +
    ggtitle(paste0('Mean of crossover distance between qtl2 and HaploQA - sample ', sample_type)) +
    geom_errorbar(aes(ymin=loc_xo_dist-st_err, ymax=loc_xo_dist+st_err), width=.2, position=position_dodge(.9)) 
  
  print(plot)
}



# computes and visualizes the cosine similarity between haploqa/qtl2  
#
# @param sample_result (object) - result object, output of pipelines
# @param rds_dir (string) - path to the saved rds files
# @param results_dir (string) - path to the saved results
# @param sample_type (string) - name of the sample being computed for (DO, CC, etc.)
#
# @output (ggplot object) - plots to show the cosine similarity for the selected sample
cos_sim_plot <- function(sample_result, rds_dir, results_dir, sample_type) {
  #sample_result <- bxd_results
  #sample_type <- 'BXD'
  print(sample_type)
  config <- fread(paste0(root, '/annotations_config.csv'))
  config_sample <- config[config$array_type == sample_type]
  founders_list <- unlist(strsplit(config_sample$founders_list, ", "))
  
  num_chr <- seq(1,19) # do only autosomes unless otherwise specified
  ### 1. genoprob_to_alleleprob - condense into founder prob (save to rds)
  founder_all_codes <- colnames(sample_result[['pr']]$`X`) # take the X chromosome - this one has everything
  all_map_codes <- seq(1:length(founder_all_codes))
  #founder_all_lookup <- setNames(all_map_codes, founder_all_codes) # make maxmarg
  founder_all_rev_lookup <- setNames(founder_all_codes, all_map_codes)
  
  test_geno <- sample_result[['ginf_haploqa']]
  qtl2_test <- genoprob_to_alleleprob(sample_result[['pr']])
  map <- sample_result[['map']]
  haploqa_founderprob_fp <- paste0(rds_dir, '/haploqa_founder_prob_', sample_type, '.rds')
  
  if(!file.exists(haploqa_founderprob_fp)) {
    print('founder prob does not exist, running calculations')
    haplo_test <- list()
    for (i in num_chr) {
      #i <- 1
      print(i)
      df <- test_geno[[i]]
      
      haplo_test[[i]] <- array(dim = c(nrow(df), length(founders_list), ncol(df)))
      rownames(haplo_test[[i]]) <- rownames(df)
      colnames(haplo_test[[i]]) <- LETTERS[seq(1,length(founders_list))]
      
      
      for (col in seq(1, ncol(df))) {
        #col <- 1
        print(col)
        
        df_test <- data.frame(sample_id = rownames(df), col_val = df[,c(col)])
        df_test[,2] <- as.data.frame(apply(df_test[c(2)], 2, function(x) founder_all_rev_lookup[x]))
        df_test <- df_test %>% separate(col_val, c("allele1", "allele2"), sep=cumsum(c(1)))
        
        # haplo
        a <- table(df_test[,c(1,2)])
        b <- table(df_test[,c(1,3)])
        
        if(length(colnames(b)) < 8) {
          b <- geno_align(b, length(founders_list))
        } 
        if((length(colnames(a)) < 8)) {
          a <- geno_align(a, length(founders_list))
        }
        
        haplo_comp <- (a + b) / 2
        
        haplo_comp <- haplo_comp[,1:length(unique(c(df_test$allele1, df_test$allele2)))]
        
        haplo_test[[i]][,,col] <- haplo_comp
        
      }
    }
    
    saveRDS(haplo_test, haploqa_founderprob_fp)
  } else {
    print('founder prob exists, reading in file')
    haplo_test <- readRDS(haploqa_founderprob_fp)
  }
  
  cos_sim_diff_fp <- paste0(rds_dir, '/cos_sim_', sample_type, '.rds')
  if(!file.exists(cos_sim_diff_fp)) {
    print('cos similarity does not exist, running calculations')
    cos_sim_all <- list()
    for (i in seq(1, 19)) {
      #i <- 1
      df <- test_geno[[i]]
      map_i <- map[[i]]
      print(i)
      cos_sim_chr <- list()
      for (col in seq(1, ncol(df))) {
        #col <- 3
        map_val <- as.numeric(map_i[col])
        # qtl2
        #print(col)
        qtl2_comp <- qtl2_test[[i]][,,col]
        ### plug the pre-calculated square roots
        haplo_comp <- haplo_test[[i]][,,col]
        
        print(qtl2_comp)
        print(haplo_comp)
        
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
    saveRDS(cos_sim_all, cos_sim_diff_fp)
  } else {
    print('cosine similarity file exists, reading in file')
    cos_sim_all <- readRDS(cos_sim_diff_fp)
  }
  df <- cos_sim_all %>% group_by(sample_id) %>% summarise(mean_chr = mean(cos_sim)) %>% as.data.frame()
  #df$chromosome <- factor(df$chromosome, num_chr)
  ### make into line plot
  ggplot(df) + aes(y = mean_chr, x = seq(1, length(mean_chr))) + geom_line() + xlab('sample_index (1-277)') # take mean of entire genome
  pdf(paste0(results_dir, '/cos_sim_', sample_type, '.pdf'))
  for (sample in unique(cos_sim_all$sample_id)) {
    #sample <- '6UY'
    print(sample)
    sample_df <- cos_sim_all[(cos_sim_all$sample_id == sample),]
    sample_df %>% group_by(sample_id) %>% summarise(mean_chr = mean(cos_sim)) %>% as.data.frame()
    sample_df$chromosome <- factor(sample_df$chromosome, num_chr)
    plot <- ggplot(sample_df) + aes(x = pos, y = cos_sim) + geom_line() + ggtitle(paste0('Sample ID: ', sample)) + 
      facet_wrap(.~chromosome)
    print(plot)
  }
  
  dev.off()
  
}



# dataframe to compare the geno codes between haploqa, qtl2 and truth model
# used to generate the summary table and boxplots
#
# @param sample_type (string) - type of sample being calculated (DO, CC, etc.)
# @param sample_results (object) - result of the sample (output of the pipelines)
# @param truth_results (object) - result of the optimal sample (output of truth_model_reconstruction)
#
# @output df_geno_all_chr (dataframe) - dataframe with geno codes
# columns of output: marker, sample_id, chr, position, gene_exp(raw ACGT geno), haplo_diplotype(haploqa output), qtl2_calls(qtl2 output), qtl2_calls_truth(optimal model output)
geno_all_comp <- function(sample_type, sample_results, truth_results) {
  num_chr <- c((1:19),"X")
  #sample <- '6UY'
  #sample_type <- 'F2'
  #sample_results <- f2_results
  #truth_results <- f2_truth_results
  founder_lookup_table <- fread(file.path(root, 'founder_lookup_table.csv'))
  founder_codes_dict <- setNames(founder_lookup_table$founder_codes, founder_lookup_table$founder_id)
  
  geno_codes <- colnames(truth_results[['pr']]$X)
  founder_all_rev_lookup <- setNames(geno_codes, seq(1, length(geno_codes)))
  
  config <- fread(paste0(root, '/annotations_config.csv'))
  config_sample <- config[config$array_type == sample_type]
  qtl2_dir <- file.path(root, config_sample$qtl2_dir)
  df_cov <- fread(file.path(qtl2_dir, 'test_covar.csv')) %>% rename('sample_id' = 'id')
  data_dir <- paste0(root, '/', config_sample$data_dir)
  summary_df <- fread(paste0(data_dir, '/', sample_type, '_summary.csv'))
  
  foundergeno <- fread(file.path(qtl2_dir, 'test_foundergeno.csv'))
  map <- fread(file.path(qtl2_dir, 'test_gmap.csv'))
  cov <- fread(file.path(qtl2_dir, 'test_covar.csv'))
  founder_geno_map <- merge(foundergeno, map, by = 'marker') 
  
  haploqa_diplotype <- get_haplotypes(summary_df, data_dir)
  
  ## make dictionary
  if (sample_type %in% c('CC', 'DO')) { ### change this to sample_type in CC or DO, they are special conditions
    founders_dict <- fread(paste0(root, '/founder_strains_table.csv'))
    founder_haplo_lookup <- setNames(founders_dict$letter, founders_dict$founder_strain)
  } else {
    unique_haplotypes <- unique(haploqa_diplotype[,c(haplotype1, haplotype2)])
    founder_haplo_lookup <- setNames(LETTERS[seq(1, length(unique(unique_haplotypes)))], unique(unique_haplotypes))
  }
  
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
  
  
  #founder_all_lookup <- setNames(names(founder_all_rev_lookup), founder_all_rev_lookup) # limit the lookup table
  #founder_all_lookup <- founder_all_lookup[names(founder_all_lookup) %in% unique(raw_geno_df$haplotype)]
  founder_all_lookup <- setNames(LETTERS[seq(1, length(unique(raw_geno_df$haplotype)))], seq(1, length(unique(raw_geno_df$haplotype))))
  
  df_geno_all_chr <- list()
  
  for (chr in num_chr) {
    #chr <- 1
    
    print(chr)
    ## raw geno
    raw_geno_acgt <- get_raw_geno(sample_type, chromosome = chr)
    raw_geno_chr_acgt <- raw_geno_acgt #%>% filter(sample_id == sample)
    
    haploqa_mom <- as.data.frame(apply(as.data.frame(sample_results[['ph_geno_haploqa']][[chr]][,,1]), 1, function(x) founder_all_lookup[x]))
    #haploqa_mom$marker <- colnames(sample_results[['ph_geno_haploqa']][[chr]][,,1])
    #haploqa_mom <- haploqa_mom %>% rename(haploqa_mom = 1)
    haploqa_dad <- as.data.frame(apply(as.data.frame(sample_results[['ph_geno_haploqa']][[chr]][,,2]), 1, function(x) founder_all_lookup[x]))
    haploqa_dad[is.na(haploqa_dad)] <- 'Y'
    #haploqa_dad$marker <- colnames(sample_results[['ph_geno_haploqa']][[chr]][,,2])
    sample_geno <- data.frame(Map(paste, haploqa_mom, haploqa_dad, MoreArgs = list(sep = "")), check.names = F)
    sample_geno$marker <- colnames(sample_results[['ph_geno_haploqa']][[chr]][,,1])
    sample_geno <- melt(sample_geno, id.vars = c("marker")) %>% rename(sample_id = variable, haplo_diplotype = value)
    
    haploqa_mom <- as.data.frame(apply(as.data.frame(sample_results[['ph_geno']][[chr]][,,1]), 1, function(x) founder_all_lookup[x]))
    #haploqa_mom$marker <- colnames(sample_results[['ph_geno']][[chr]][,,1])
    #haploqa_mom <- haploqa_mom %>% rename(haploqa_mom = 1)
    haploqa_dad <- as.data.frame(apply(as.data.frame(sample_results[['ph_geno']][[chr]][,,2]), 1, function(x) founder_all_lookup[x]))
    haploqa_dad[is.na(haploqa_dad)] <- 'Y'
    #haploqa_dad$marker <- colnames(sample_results[['ph_geno']][[chr]][,,2])
    #haploqa_dad <- haploqa_dad %>% rename(haploqa_dad = 1)
    sample_geno_qtl2 <- data.frame(Map(paste, haploqa_mom, haploqa_dad, MoreArgs = list(sep = "")), check.names = F)
    sample_geno_qtl2$marker <- colnames(sample_results[['ph_geno']][[chr]][,,1])
    sample_geno_qtl2 <- melt(sample_geno_qtl2, id.vars = c("marker")) %>% rename(sample_id = variable, qtl2_calls = value)
    
    # only need qtl2 for truth models
    if(is.na(dim(truth_results[['ph_geno']][[chr]])[3])) {
      sample_geno_truth <- as.data.frame(apply(as.data.frame(truth_results[['ph_geno']][[chr]]), 1, function(x) founder_all_rev_lookup[x]))
      sample_geno_truth$marker <- colnames(truth_results[['ph_geno']][[chr]])
      sample_geno_truth <- melt(sample_geno_truth, id.vars = c("marker")) %>% rename(sample_id = variable, qtl2_calls_truth = value)
      if(sample_type == 'CC' & chr == 'X') {
        sample_geno_truth <- merge(sample_geno_truth, cov %>% rename(sample_id = id), on = 'sample_id') %>% 
          separate(qtl2_calls_truth,  c("allele1", "allele2"), sep=cumsum(c(1,1))) %>%
          mutate(allele2 = ifelse(Sex == 'male', 'Y', allele2)) %>%
          unite(qtl2_calls_truth, allele1, allele2, sep = '')
        # if male, change to Y
        
      }
      
      #haploqa_df
    } else {
      haploqa_mom <- as.data.frame(apply(as.data.frame(truth_results[['ph_geno']][[chr]][,,1]), 1, function(x) founder_all_lookup[x]))
      #haploqa_mom$marker <- colnames(truth_results[['ph_geno_haploqa']][[chr]][,,1])
      #haploqa_mom <- haploqa_mom %>% rename(haploqa_mom = 1)
      haploqa_dad <- as.data.frame(apply(as.data.frame(truth_results[['ph_geno']][[chr]][,,2]), 1, function(x) founder_all_lookup[x]))
      haploqa_dad[is.na(haploqa_dad)] <- 'Y'
      #haploqa_dad$marker <- colnames(truth_results[['ph_geno_haploqa']][[chr]][,,2])
      #haploqa_dad <- haploqa_dad %>% rename(haploqa_dad = 1)
      sample_geno_truth <- data.frame(Map(paste, haploqa_mom, haploqa_dad, MoreArgs = list(sep = "")), check.names = F)
      sample_geno_truth$marker <- colnames(truth_results[['ph_geno']][[chr]][,,1])
      sample_geno_truth <- melt(sample_geno_truth, id.vars = c("marker")) %>% rename(sample_id = variable, qtl2_calls_truth = value)
      
    }
    
    df_geno_all_acgt <- merge(merge(merge(raw_geno_chr_acgt, sample_geno, by = c('marker', 'sample_id'), sort = F), sample_geno_qtl2, by = c('marker', 'sample_id'), sort = F), sample_geno_truth, by = c('marker', 'sample_id'), sort = F) %>% select(marker, sample_id, chr, pos, everything()) %>% select(-any_of(c('Sex')))
    #df_geno_all_acgt$is_different <- NA
    
    df_geno_all_acgt[,6:ncol(df_geno_all_acgt)] <- lapply(df_geno_all_acgt[,6:ncol(df_geno_all_acgt)], function(col) {
      rev_col = stri_reverse(col)
      ifelse(rev_col %in% founder_codes_dict, rev_col, col)
    })
    
    df_geno_all_chr[[chr]] <- df_geno_all_acgt
  }
  
  df_geno_all_chr <- rbindlist(df_geno_all_chr)
  
  return(df_geno_all_chr)
}


# function to compute the percentage of markers with disagreeing geno codes based on the geno_comp_all function
# two functions need to be separate because need both the geno code object and the quantitative metrics
# @param geno_comp_df (data.frame) - output of geno_comp_all
#
# @output marker_comp_diff (data.frame) - dataframe with calculated percentages of disagreeing geno codes each between qtl2, haploqa and optimal model
marker_comp_diff <- function(geno_comp_df) {
  qtl2_diffs <- geno_comp_df %>% group_by(sample_id) %>% summarise(qtl2_pct_diff = sum(qtl2_calls_truth != qtl2_calls)/n()) %>% as.data.frame()
  haploqa_diffs <- geno_comp_df %>% group_by(sample_id) %>% summarise(haploqa_pct_diff = sum(qtl2_calls_truth != haplo_diplotype)/n()) %>% as.data.frame()
  qtl2_haploqa_diffs <- geno_comp_df %>% group_by(sample_id) %>% summarise(haploqa_qtl2_pct_diff = sum(qtl2_calls != haplo_diplotype)/n()) %>% as.data.frame()
  sample_diffs <- merge(merge(qtl2_diffs, haploqa_diffs, by = 'sample_id'), qtl2_haploqa_diffs, by = 'sample_id')
  return(marker_comp_diff)
}

