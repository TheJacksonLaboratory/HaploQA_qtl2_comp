### Contains three major pipelines for haplotype reconstructions

library(rstudioapi)
root <- dirname(getSourceEditorContext()$path)

# environment set-up
source(paste0(root,"/input_data_retrieval.R"))
source(paste0(root,"/qtl2_input_data_prep.R"))
source(paste0(root,"/qtl2_sim_comp.R"))
source(paste0(root,"/get_comp_results.R"))
source(paste0(root,"/utils.R"))


# function to initiate the haplotype reconstruction pipeline
# @param sample_type (string) - type of sample, such as CC, DO, F2, etc.
# @param list_pheno (list) - list of phenotypes to be simulated
# @param qtl2_file_gen (True/False) - toggle to set whether to (re)generate qtl2 input files
# @param samples_gen (True/False) - toggle to set whether to (re)generate individual sample  files
#
# @return results (list of dataframes) - a list containing all objects necessary for computations
haplotype_reconstruction_pipeline <- function(sample_type, list_pheno, qtl2_file_gen, samples_gen) {
  
  #sample_type <- 'CC'
  #truth_model = T
  ## results directory
  results_dir <- file.path(root, 'results')
  dir.create(results_dir, showWarnings = FALSE) 
  
  num_chr <- c((1:19),"X")
  
  ### config file
  config <- fread(paste0(root, '/annotations_config.csv'))
  
  ### Environment
  config_sample <- config[config$array_type == sample_type]
  marker_type <- config_sample$marker_type
  founders_list <- unlist(strsplit(config_sample$founders_list, ", "))
  n_founders <- length(founders_list)
  ngen <- config_sample$ngen
  sample_url <- config_sample$url
  founder_url <- config_sample$founder_url
  founder_filename <- config_sample$founder_filename
  
  # data output directory
  data_dir <- file.path(root, config_sample$data_dir)
  dir.create(data_dir, showWarnings = FALSE) 
  
  ## create a data directory for qtl2 input data
  qtl2_dir <- file.path(root, config_sample$qtl2_dir) # name of desired output folder
  # create if folder not exist
  dir.create(qtl2_dir, showWarnings = FALSE)
  
  # filepaths to save any rds file
  rds_dir <- file.path(results_dir, 'RDS')
  dir.create(rds_dir, showWarnings = FALSE)
  
  # annotation file
  annot_file <- config_sample$annot_file
  
  ### control file
  control_fp <- paste0(qtl2_dir, '/test.json')
  if (file.exists(control_fp) == FALSE) {
    file.create(control_fp)
  }
  
  ## summary file
  summary_df_fp <- paste0(data_dir, '/', sample_type, '_summary.csv')
  
  ## list of samples to exclude, if any
  exclude_list <- unlist(strsplit(config_sample$exclude_samples, ", "))
  
  
  control_file <- paste0('{
  "description": "HaploQA data - qtl2 test run",
  "crosstype": "genail', n_founders, '",
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
  "alleles": [', substring(gsub("[()]", "", toString(list(LETTERS[seq(1, n_founders)]))),2), '], 
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
}')
  
  ###############################################################################
  ### Part 1 - get results from HaploQA
  # set main URL domain
  url_domain <- 'http://haploqa-dev.jax.org/' 
  
  #### implementations
  # summary table
  if (!file.exists(summary_df_fp)) {
    print(paste0('Working on summary file:'))
    #### html file prep
    # read html file
    haploqa_html <- read_html(sample_url)
    # extract information from html file
    html_temp <- haploqa_html %>% html_nodes("a") %>% html_attr("href")
    url_list <- paste0(url_domain, html_temp[grepl('/sample', html_temp)])
    summary_df <- sample_summary_scrape(haploqa_html, url_list, marker_type)
    print('Writing to directory')
    write.csv(summary_df, paste0(data_dir, '/', sample_type, '_summary.csv'), row.names = FALSE)
    
  } else {
    summary_df <- fread(summary_df_fp)
  }
  
  #url_ind <- url_ind[!url_ind %in% exclude_list]
  
  if(samples_gen == T) {
    #### html file prep
    # read html file
    haploqa_html <- read_html(sample_url)
    # extract information from html file
    html_temp <- haploqa_html %>% html_nodes("a") %>% html_attr("href")
    url_list <- paste0(url_domain, html_temp[grepl('/sample', html_temp)])
    
    # list of urls to generate samples from
    url_ind <- unique(summary_df$`Sample Filepath`)
    
    # individual samples
    inc = 0
    for (url in url_ind) {
      inc = inc + 1 # increment
      file <- sample_individual_scrape(url, url_domain)
      print(paste0('Working on file ', inc, '/', length(url_ind), ': ', file))
      sample_df_save <- as.data.frame(content(GET(file)))
      print('Writing to directory')
      file_name <- unlist(strsplit(file, '/'))[6]
      GET(file, write_disk(paste0(data_dir, '/', file_name), overwrite = TRUE), show_col_types = FALSE)
      
      
    }
  }
  
  ### Part 2 - convert results to qtl2 input format
  # implement function
  qtl2_objects_fp <- paste0(rds_dir, "/qtl2_objects_", sample_type, ".rds")
  if (!file.exists(qtl2_objects_fp)) {
    print('qtl2 input files do not exist, running calculations')
    file_output <- get_qtl2_input(data_dir, sample_type, annot_file, qtl2_dir, summary_df, list_pheno, ngen, founders_list, marker_type, exclude_list, founder_url)
    saveRDS(file_output, file = qtl2_objects_fp)
  } else {
    print('qtl2 input files exist, reading in file')  
    file_output <- readRDS(qtl2_objects_fp)
  }
  
  # unpack
  df_geno <- file_output[[1]]
  df_gmap <- file_output[[2]]
  df_gmap <- sort_chr(df_gmap, c((1:19),"X")) # put chromosomes in order
  df_pmap <- file_output[[3]]
  df_pmap <- sort_chr(df_pmap, c((1:19),"X"))
  df_pheno <- file_output[[4]]
  df_covar <- file_output[[5]]
  df_crossinfo <- file_output[[7]] 
  df_foundergeno <- file_output[[6]] # with strain id metadata
  founder_haplo_lookup <- file_output[[8]]
  
  
  if (qtl2_file_gen == T) {
    # write out files
    write.csv(df_geno, paste0(qtl2_dir, '/test_geno.csv'), row.names = F)
    write.csv(df_gmap, paste0(qtl2_dir, '/test_gmap.csv'), row.names = F)
    write.csv(df_pmap, paste0(qtl2_dir, '/test_pmap.csv'), row.names = F)
    write.csv(df_pheno, paste0(qtl2_dir, '/test_pheno.csv'), row.names = F)
    write.csv(df_covar, paste0(qtl2_dir, '/test_covar.csv'), row.names = F)
    write.csv(df_foundergeno, paste0(qtl2_dir, '/test_foundergeno.csv'), row.names = F)
    write.csv(df_crossinfo, paste0(qtl2_dir, '/test_crossinfo.csv'), row.names = F)
    writeLines(control_file, control_fp)
  }
  
  ### get results for the sample
  results <- qtl2_metrics_comp(rds_dir, qtl2_dir, sample_type, results_dir, data_dir, n_founders, founder_haplo_lookup, truth_model = F)
  
  return(results)
}


# function to initiate the haplotype reconstruction pipeline for one individual sample/mouse
# @param sample_type (string) - type of sample, such as CC, DO, F2, etc.
# @param list_pheno (list) - list of phenotypes to be simulated
# @param qtl2_file_gen (True/False) - toggle to set whether to (re)generate qtl2 input files
# @param samples_gen (True/False) - toggle to set whether to (re)generate individual sample  files
#
# @return results (list of dataframes) - a list containing all objects necessary for computations

## pipeline for one
sample_haplotype_reconstruction <- function(sample_type, sample_name, samples_gen, qtl2_file_gen, num_founder_test = NULL) {
  ## results directory
  results_dir <- file.path(root, 'results')
  dir.create(results_dir, showWarnings = FALSE) 
  
  num_chr <- c((1:19),"X")
  
  ### config file
  config <- fread(paste0(root, '/annotations_config.csv'))
  
  ### Environment
  config_sample <- config[config$array_type == sample_type]
  marker_type <- config_sample$marker_type
  founders_list <- unlist(strsplit(config_sample$founders_list, ", "))
  n_founders <- length(founders_list)
  ngen <- config_sample$ngen
  sample_url <- config_sample$url
  founder_url <- config_sample$founder_url
  founder_filename <- config_sample$founder_filename
  
  # data output directory
  data_dir <- file.path(root, config_sample$data_dir)
  dir.create(data_dir, showWarnings = FALSE) 
  
  ## create a data directory for qtl2 input data
  qtl2_dir <- file.path(root, config_sample$qtl2_dir) # name of desired output folder
  # create if folder not exist
  dir.create(qtl2_dir, showWarnings = FALSE)
  
  # filepaths to save any rds file
  rds_dir <- file.path(results_dir, 'RDS')
  dir.create(rds_dir, showWarnings = FALSE)
  
  # annotation file
  annot_file <- config_sample$annot_file
  
  ### sample folder
  sample_dir <- file.path(qtl2_dir, sample_name)
  dir.create(sample_dir, showWarnings = FALSE)
  
  ### control file
  control_fp <- paste0(sample_dir, '/test.json')
  if (file.exists(control_fp) == FALSE) {
    file.create(control_fp)
  }
  
  
  ## summary file
  summary_df_fp <- paste0(sample_dir, '/', sample_type, '_', sample_name, '_summary.csv')
  
  ## list of samples to exclude, if any
  exclude_list <- unlist(strsplit(config_sample$exclude_samples, ", "))
  
  
  ###############################################################################
  ### Part 1 - get results from HaploQA
  # set main URL domain
  url_domain <- 'http://haploqa-dev.jax.org/' 
  
  #### implementations
  # summary table
  if (!file.exists(summary_df_fp)) {
    print(paste0('Working on summary file:'))
    #### html file prep
    # read html file
    haploqa_html <- read_html(sample_url)
    # extract information from html file
    html_temp <- haploqa_html %>% html_nodes("a") %>% html_attr("href")
    url_list <- paste0(url_domain, html_temp[grepl('/sample', html_temp)])
    summary_df <- sample_summary_scrape(haploqa_html, url_list, marker_type)
    summary_df <- summary_df[summary_df$ID == sample_name,]
    print('Writing to directory')
    write.csv(summary_df, summary_df_fp, row.names = FALSE)
    
  } else {
    summary_df <- fread(summary_df_fp)
    summary_df <- summary_df[summary_df$ID == sample_name,]
  }
  
  #url_ind <- url_ind[!url_ind %in% exclude_list]
  
  if(samples_gen == T) {
    #### html file prep
    # read html file
    haploqa_html <- read_html(sample_url)
    # extract information from html file
    html_temp <- haploqa_html %>% html_nodes("a") %>% html_attr("href")
    url_list <- paste0(url_domain, html_temp[grepl('/sample', html_temp)])
    
    # list of urls to generate samples from
    url_ind <- unique(summary_df$`Sample Filepath`)
    
    # individual samples
    inc = 0
    for (url in url_ind) {
      inc = inc + 1 # increment
      file <- sample_individual_scrape(url, url_domain)
      print(paste0('Working on file ', inc, '/', length(url_ind), ': ', file))
      sample_df_save <- as.data.frame(content(GET(file)))
      print('Writing to directory')
      file_name <- unlist(strsplit(file, '/'))[6]
      GET(file, write_disk(paste0(data_dir, '/', file_name), overwrite = TRUE), show_col_types = FALSE)
    }
  }
  
  file_output <- get_qtl2_input(data_dir, sample_type, annot_file, qtl2_dir, summary_df, list_pheno, ngen, founders_list, marker_type, exclude_list = NULL, founder_url, founder_filename)
  # unpack
  df_geno <- file_output[[1]]
  df_gmap <- file_output[[2]]
  df_gmap <- sort_chr(df_gmap, c((1:19),"X")) # put chromosomes in order
  df_pmap <- file_output[[3]]
  df_pmap <- sort_chr(df_pmap, c((1:19),"X"))
  df_pheno <- file_output[[4]]
  df_covar <- file_output[[5]]
  df_foundergeno <- file_output[[6]] # with strain id metadata
  df_crossinfo <- file_output[[7]]
  founder_haplo_lookup <- file_output[[8]]
  
  
  ## unique foundergeno situation
  if (!is.null(num_founder_test)) {
    ## should only do this for gigamuga?
    ## gigamuga founders
    url_domain <- 'http://haploqa-dev.jax.org/' 
    haploqa_html <- read_html('http://haploqa-dev.jax.org/tag/UNC_Villena_GIGMUGV01_20141012_FinalReport.html')
    # extract information from html file
    html_temp <- haploqa_html %>% html_nodes("a") %>% html_attr("href")
    url_list <- paste0(url_domain, html_temp[grepl('/sample', html_temp)])
    founder_sum <- sample_summary_scrape(haploqa_html, url_list, 'GigaMUGA')
    founder_sum_gm <- founder_sum[founder_sum$`Haplotype Candidate` == 'True',] %>% rename(sample_id = ID, original_sample_id = 'Secondary IDs')
    founder_gm <- fread(file.path(root, 'UNC_Villena_founder_samples.csv'))
    gm_founders <- merge(founder_sum_gm, founder_gm, by = 'original_sample_id') %>% 
      mutate(gene_exp = paste(allele1_fwd, allele2_fwd, sep = '')) %>%
      rename(sample_id = 'sample_id.y') %>%
      select('snp_id', 'gene_exp', 'Strain Name') 
    
    # take cast out for test
    #founders_select <- names(founder_haplo_lookup)[names(founder_haplo_lookup) != 'CAST/EiJ']
    #founders_test_lookup <- data.frame(strain=founders_select, letter=LETTERS[1:length(founders_select)]) %>% rename('Strain Name'=strain)
    
    founders_select <- sample(unique(gm_founders$`Strain Name`), num_founder_test) # randomly selects a selected number of samples
    founders_test_lookup <- data.frame(strain=founders_select, letter=LETTERS[1:length(founders_select)]) %>% rename('Strain Name'=strain)
    
    selected_founders <- merge(gm_founders, founders_test_lookup, by = 'Strain Name')
    
    selected_founders <- selected_founders %>% select(snp_id, letter, gene_exp) %>%
      rename(strain = letter, marker = snp_id)
    selected_founders <- dcast(selected_founders, marker~strain, value.var="gene_exp")
    rownames(selected_founders) <- selected_founders$marker
    
    # encode
    dir <- "https://raw.githubusercontent.com/kbroman/MUGAarrays/main/UWisc/" # always the same
    annot_df <- read.csv(paste0(dir, annot_file))
    
    annot_encode_df <- annot_df %>% select(marker, chr, snp) %>% # select allele columns
      separate(snp, c("allele1", "allele2"), sep=cumsum(c(1)))
    founder_encode_annot <- merge(annot_encode_df, selected_founders, by = 'marker') %>% select(allele1, allele2)
    df_founders_encoded <- as.data.frame(encode_geno(selected_founders[,-1], founder_encode_annot))
    df_founders_encoded$marker <- rownames(selected_founders)
    marker_order <- unique(df_gmap$marker)
    df_founders_encoded <- df_founders_encoded %>% filter(marker %in% marker_order) %>% arrange(factor(marker, levels = marker_order))
    df_foundergeno <- df_founders_encoded %>% select(marker, everything())
  }
  
  if(sample_type == 'MiniMUGA') {
    url_domain <- 'http://haploqa-dev.jax.org/' 
    haploqa_html <- read_html(sample_url)
    # extract information from html file
    html_temp <- haploqa_html %>% html_nodes("a") %>% html_attr("href")
    url_list <- paste0(url_domain, html_temp[grepl('/sample', html_temp)])
    founder_sum <- sample_summary_scrape(haploqa_html, url_list, marker_type)
    founder_sum <- founder_sum[founder_sum$`Haplotype Candidate` == 'True',] %>% rename(sample_id = ID, original_sample_id = 'Secondary IDs')
    fp_founders <- file.path(root, 'MiniMUGA_founder_samples.csv')
    
    founders_total <- get_founder_data(sample_url, founder_sum$`Sample Filepath`, sample_type, data_dir, founder_filename)
    
    founders_list <- names(founder_haplo_lookup)
    
    founders_df <- merge(founders_total, founder_sum, by = 'original_sample_id') %>% 
      mutate(gene_exp = paste(allele1_fwd, allele2_fwd, sep = '')) %>%
      rename(sample_id = 'sample_id.y') %>%
      select('sample_id', 'snp_id', 'gene_exp', 'Strain Name', 'original_sample_id')
    
    founders_test <- founders_df[founders_df$`Strain Name` %in% founders_list,]
    haplo_lookup <- data.frame(strain=names(founder_haplo_lookup), letter=founder_haplo_lookup) %>% rename('Strain Name'=strain)
    founders_test <- merge(founders_test, haplo_lookup, by = 'Strain Name')
    
    test_foundergeno <- founders_test %>% select(snp_id, letter, gene_exp) %>% 
      rename(strain = letter, marker = snp_id) %>% unique()
    test_foundergeno <- dcast(test_foundergeno, marker~strain, value.var="gene_exp")
    rownames(test_foundergeno) <- test_foundergeno$marker
    
    # put the allele annotations in according order
    ### read annotations
    dir <- "https://raw.githubusercontent.com/kbroman/MUGAarrays/main/UWisc/" # always the same
    annot_df <- read.csv(paste0(dir, annot_file))
    # encode_genome function on geno and founder from qtl2 convert package - args: matrix of genotypes, matrixs of two alleles
    annot_encode_df <- annot_df %>% select(marker, chr, snp) %>% # select allele columns
      separate(snp, c("allele1", "allele2"), sep=cumsum(c(1)))
    founder_encode_annot <- merge(annot_encode_df, test_foundergeno, by = 'marker') %>% select(allele1, allele2)
    df_founders_encoded <- as.data.frame(encode_geno(test_foundergeno[,-1], founder_encode_annot))
    df_founders_encoded$marker <- rownames(test_foundergeno)
    marker_order <- unique(df_gmap$marker)
    df_founders_encoded <- df_founders_encoded %>% filter(marker %in% marker_order) %>% arrange(factor(marker, levels = marker_order))
    df_foundergeno <- df_founders_encoded %>% select(marker, everything())
  }
  
  ## text of control file
  ### alleles need to be consistent with genail
  control_file <- paste0('{
  "description": "HaploQA data - qtl2 test run",
  "crosstype": "genail', length(founder_haplo_lookup), '",
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
  "alleles": [', substring(gsub("[()]", "", toString(list(LETTERS[seq(1, ncol(df_foundergeno) - 1)]))),2), '], 
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
}')
  
  if (!is.null(num_founder_test)) {
    ## text of control file
    ### alleles need to be consistent with genail
    control_file <- paste0('{
      "description": "HaploQA data - qtl2 test run",
      "crosstype": "genail', num_founder_test, '",
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
      "alleles": [', substring(gsub("[()]", "", toString(list(LETTERS[seq(1, ncol(df_foundergeno) - 1)]))),2), '], 
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
    }')
    
    cols_list <- LETTERS[seq(1, num_founder_test)]
    df_crossinfo <- df_crossinfo %>% select(id, ngen, all_of(cols_list))
    df_crossinfo[ ,c(cols_list)] <- 1
  }
  
  # output
  ## genail, json alleles, founder_geno, cross_info
  if (qtl2_file_gen == T) {
    # write out files
    write.csv(df_geno, paste0(sample_dir, '/test_geno.csv'), row.names = F)
    write.csv(df_gmap, paste0(sample_dir, '/test_gmap.csv'), row.names = F)
    write.csv(df_pmap, paste0(sample_dir, '/test_pmap.csv'), row.names = F)
    write.csv(df_pheno, paste0(sample_dir, '/test_pheno.csv'), row.names = F)
    write.csv(df_covar, paste0(sample_dir, '/test_covar.csv'), row.names = F)
    write.csv(df_foundergeno, paste0(sample_dir, '/test_foundergeno.csv'), row.names = F)
    write.csv(df_crossinfo, paste0(sample_dir, '/test_crossinfo.csv'), row.names = F)
    writeLines(control_file, control_fp)
  }
  
  founder_lookup_table <- fread(file.path(root, 'founder_lookup_table.csv'))
  founder_all_rev_lookup <- setNames(founder_lookup_table$founder_codes, founder_lookup_table$founder_id)
  
  # calculations
  cross_list <- list()
  cross_file <- read_cross2(control_fp)
  cross_list[['cross']] <- cross_file
  pr <- calc_genoprob(cross=cross_file, map=cross_file$gmap, error_prob=0.002, cores = 2)
  #test_error <- calc_errorlod(cross_file, pr)
  #pdf(paste0(results_dir, '/error_lod_PX_genail3.pdf'))
  #for (chr in num_chr) {
  #  print(chr)
  #  chr_error <- test_error[[chr]] %>% t() %>% as.data.frame()
  #  chr_error$marker <- rownames(chr_error)
  #  pmap <- cross_file$pmap[[chr]] %>% as.data.frame()
  #  pmap$marker <- rownames(pmap)
  #  plot_error <- merge(pmap, chr_error, by = 'marker', sort = F) %>% rename('pmap' = '.')
  #  plot <- ggplot(plot_error) + aes(y = `PX`, x = pmap) + geom_point() + ggtitle(paste0('Mahaffey PX genail3 - Error LOD plot on chromosome ', chr)) + xlab('pmap') + ylab('Error LOD')
  #  print(plot)
  #}
  #dev.off()
  
  cross_list[['pr']] <- pr
  ginf <- maxmarg(pr, minprob = 0.01)
  cross_list[['ginf']] <- ginf
  ph_geno <- guess_phase(cross_file, ginf)
  cross_list[['ph_geno']] <- ph_geno
  pos <- locate_xo(ginf, cross_file$gmap)
  cross_list[['pos']] <- pos
  
  ginf_haploqa <- maxmarg_sim(summary_df, data_dir, num_chr, founder_all_rev_lookup, cross = cross_file, sample_dir, sample_type, n_founders, founder_haplo_lookup, pr)
  cross_list[['ginf_haploqa']] <- ginf_haploqa
  ph_geno_haploqa <- guess_phase(cross_file, ginf_haploqa)
  cross_list[['ph_geno_haploqa']] <- ph_geno_haploqa
  
  #plot_onegeno_test(ph_geno, ph_geno_haploqa, map)
  
  # save for shiny
  # map
  map <- cross_file$gmap
  map_fp <- file.path(sample_dir, 'map.rds')
  saveRDS(map, file = map_fp)
  # qtl2 phased geno
  ph_geno_fp <- file.path(sample_dir, 'ph_geno.rds')
  saveRDS(ph_geno, file = ph_geno_fp)
  # haploqa phased geno
  ph_geno_haploqa_fp <- file.path(sample_dir, 'ph_geno_haploqa.rds')
  saveRDS(ph_geno_haploqa, file = ph_geno_haploqa_fp)
  
  #pdf(paste0(sample_dir, '/', sample_name, '_genoplot.pdf'))
  #plot_onegeno_test(ph_geno, ph_geno_haploqa, cross_file$gmap, ind = 1, 
  #shift = TRUE, main = paste0('Geno-wide genotypes of individual ', '1'),
  #sub = 'Left - qtl2, Right - HaploQA')
  #dev.off()
  #print(cross_list)
  return(cross_list)
}



# function to initiate the haplotype reconstruction pipeline for optimal models
# please refer to annotations_config.csv to determine which cross was used as optimal for each crosstype
# @param sample_type (string) - type of sample, such as CC, DO, F2, etc.
# @param list_pheno (list) - list of phenotypes to be simulated
# @param qtl2_file_gen (True/False) - toggle to set whether to (re)generate qtl2 input files
# @param samples_gen (True/False) - toggle to set whether to (re)generate individual sample  files
#
# @return results (list of dataframes) - a list containing all objects necessary for computations
truth_model_reconstruction <- function(sample_type, list_pheno, qtl2_file_gen, samples_gen) {
  ## results directory
  results_dir <- file.path(root, 'results')
  dir.create(results_dir, showWarnings = FALSE) 
  
  num_chr <- c((1:19),"X")
  
  ### config file
  config <- fread(paste0(root, '/annotations_config.csv'))
  
  ### Environment
  config_sample <- config[config$array_type == sample_type]
  marker_type <- config_sample$marker_type
  founders_list <- unlist(strsplit(config_sample$founders_list, ", "))
  n_founders <- length(founders_list)
  ngen <- config_sample$ngen
  sample_url <- config_sample$url
  founder_url <- config_sample$founder_url
  
  # data output directory
  data_dir <- file.path(root, config_sample$data_dir)
  dir.create(data_dir, showWarnings = FALSE) 
  
  ## create a data directory for qtl2 input data
  qtl2_dir <- paste0(root, '/', config_sample$qtl2_dir, '_truth')
  # create if folder not exist
  dir.create(qtl2_dir, showWarnings = FALSE)
  
  # filepaths to save any rds file
  rds_dir <- file.path(results_dir, 'RDS')
  dir.create(rds_dir, showWarnings = FALSE)
  
  # annotation file
  annot_file <- config_sample$annot_file
  
  ### control file
  control_fp <- paste0(qtl2_dir, '/test.json')
  if (file.exists(control_fp) == FALSE) {
    file.create(control_fp)
  }
  
  ## summary file
  summary_df_fp <- paste0(data_dir, '/', sample_type, '_summary.csv')
  summary_df <- fread(summary_df_fp)
  
  ## list of samples to exclude, if any
  exclude_list <- unlist(strsplit(config_sample$exclude_samples, ", "))
  
  if (sample_type %in% c('CC', 'DO')) {
    control_file <- paste0('{
    "description": "HaploQA data - qtl2 test run",
    "crosstype": "', config_sample$truth_model, '",
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
    "alleles": [', substring(gsub("[()]", "", toString(list(LETTERS[seq(1, n_founders)]))),2), '], 
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
    }')}
  
  if(sample_type %in% c('F2')) {
    control_file <- paste0('{
    "description": "HaploQA data - qtl2 test run",
    "crosstype": "', config_sample$truth_model, '",
    "sep": ",",
    "na.strings": ["-", "NA"],
    "comment.char": "#",
    "geno": "test_geno.csv",
    "gmap": "test_gmap.csv",
    "pmap": "test_pmap.csv",
    "pheno": "test_pheno.csv",
    "covar": "test_covar.csv",
    "genotypes": {
      "AA": 1,
      "AB": 2,
      "BB": 3
    },
    "x_chr": "X",
    "alleles": ["A", "B"],
    "geno_transposed": true,
    "sex": {
      "covar": "Sex",
      "female": "female",
      "male": "male"
    },
    "cross_info": {
      "covar":"cross_direction", "(AxB)x(AxB)":0, "(BxA)x(BxA)":1
    }
    }')}
  
  if(sample_type %in% c('BXD')) {
    control_file <- paste0('{
    "description": "HaploQA data - qtl2 test run",
    "crosstype": "', config_sample$truth_model, '",
    "sep": ",",
    "na.strings": ["-", "NA"],
    "comment.char": "#",
    "geno": "test_geno.csv",
    "gmap": "test_gmap.csv",
    "pmap": "test_pmap.csv",
    "pheno": "test_pheno.csv",
    "covar": "test_covar.csv",
    "genotypes": {
      "AA": 1,
      "AB": 2,
      "BB": 3
    },
    "x_chr": "X",
    "alleles": ["A", "B"],
    "geno_transposed": true,
    "sex": {
      "covar": "Sex",
      "female": "female",
      "male": "male"
    },
    "cross_info": {
      "file": "test_crossinfo.csv"
    }
    }')}
  
  
  qtl2_objects_fp <- paste0(rds_dir, "/qtl2_objects_", sample_type, ".rds")
  stopifnot("qtl2 input file does not exist for the cross this truth model is for, please check input file" = !file.exists(qtl2_objects_fp) == FALSE) 
  file_output <- readRDS(qtl2_objects_fp)
  
  # unpack
  df_geno <- file_output[[1]]
  df_gmap <- file_output[[2]]
  df_gmap <- sort_chr(df_gmap, c((1:19),"X")) # put chromosomes in order
  df_pmap <- file_output[[3]]
  df_pmap <- sort_chr(df_pmap, c((1:19),"X"))
  df_pheno <- file_output[[4]]
  df_covar <- file_output[[5]]
  df_foundergeno <- file_output[[6]]
  df_crossinfo <- file_output[[7]]
  founder_haplo_lookup <- file_output[[8]]
  
  if(sample_type %in% c('F2', 'BXD')) {
    ## encode geno
    new_geno = df_geno
    fgeno = df_foundergeno
    fgeno[fgeno[,2] == 'H',2] = '-'
    fgeno[fgeno[,3] == 'H',3] = '-'
    rownames(new_geno) = new_geno$marker
    new_geno = as.matrix(new_geno[,-1])
    
    wh = which(fgeno[,2] == fgeno[,3] | fgeno[,2] == '-' | fgeno[,3] == '-')
    fgeno[wh, 2] = '00'
    fgeno[wh, 3] = '00'
    new_geno[new_geno == fgeno[,'A']] = 'AA'
    new_geno[new_geno == fgeno[,'B']] = 'BB'
    new_geno[new_geno == 'H']         = 'AB'
    new_geno[new_geno == '-']         = '--'
    
    fgeno[fgeno == '00'] = '-'
    
    wh = which(fgeno[,2] == '-' & fgeno[,3] == '-')
    new_geno[wh,] = matrix('--', nrow = length(wh), ncol = ncol(new_geno))
    
    new_geno = data.frame(marker = rownames(new_geno), new_geno)
    colnames(new_geno) = sub('^X', '', colnames(new_geno))
    
    df_geno <- new_geno
    
    ## covar file
    if(sample_type == 'F2'){
      df_covar$cross_direction <- '(BxA)x(BxA)'
    }
    if(sample_type == 'BXD'){
      df_crossinfo$ngen <- 3
      df_crossinfo$cross_direction <- '0'
      df_crossinfo <- df_crossinfo %>% select(id, ngen, cross_direction)
    }
    
  }
  
  ### cross info
  # with funnel metadata
  if (sample_type == 'CC') {
    sample_ci <- fread('https://raw.githubusercontent.com/kbroman/qtl2data/main/CC/cc_crossinfo.csv')
    sample_ci$strain = substr(sample_ci$id, 1, 5)
    summary_df$strain = substr(substring(summary_df$`Strain Name`, regexpr("CC", summary_df$`Strain Name`)),1,5)
    df_crossinfo <- merge(summary_df, sample_ci, on = 'strain', sort = F, all.x = T) %>% select(ID, LETTERS[1:8]) %>% rename(id = ID)## missing some from funnel
    df_crossinfo <- df_crossinfo %>% #mutate_at(.vars = LETTERS[1:8], funs(ifelse(is.na(.), '1', .)))
      mutate(A = ifelse(is.na(A), '1', A)) %>% 
      mutate(B = ifelse(is.na(B), '2', B)) %>%
      mutate(C = ifelse(is.na(C), '3', C)) %>%
      mutate(D = ifelse(is.na(D), '4', D)) %>%
      mutate(E = ifelse(is.na(E), '5', E)) %>%
      mutate(F = ifelse(is.na(F), '6', F)) %>%
      mutate(G = ifelse(is.na(G), '7', G)) %>%
      mutate(H = ifelse(is.na(H), '8', H))
  }
  
  if (qtl2_file_gen == T) {
    # write out files
    write.csv(df_geno, paste0(qtl2_dir, '/test_geno.csv'), row.names = F)
    write.csv(df_gmap, paste0(qtl2_dir, '/test_gmap.csv'), row.names = F)
    write.csv(df_pmap, paste0(qtl2_dir, '/test_pmap.csv'), row.names = F)
    write.csv(df_pheno, paste0(qtl2_dir, '/test_pheno.csv'), row.names = F)
    write.csv(df_covar, paste0(qtl2_dir, '/test_covar.csv'), row.names = F)
    write.csv(df_foundergeno, paste0(qtl2_dir, '/test_foundergeno.csv'), row.names = F)
    ### f2 and bxd don't need cross_info files
    if (!sample_type %in% c('F2')) {
      write.csv(df_crossinfo, paste0(qtl2_dir, '/test_crossinfo.csv'), row.names = F) # only output for CC and DO
    }
    
    writeLines(control_file, control_fp)
  }
  
  ### get results for the sample
  results <- save_rds(rds_dir, qtl2_dir, sample_type, results_dir, data_dir, n_founders, founder_haplo_lookup, truth_model = T)
  
  return(results)
}
