### Contains all general utility functions used across the board
### function directory
## 1. get_haplotypes
## 2. sort_chr - sort chromosomes from 1-19 and X as needed
## 3. make_cross - adds crosstype attributes to an object to make it a cross (as needed)
## 4. shiny_viz_input
## 5. get_raw_geno - 
## 6. plot_onegeno_test - modified version of plot_onegeno in qtl2, plots two genos instead of one
## 7. addgenorect - util function that adds rectangles based on geno positions, used in plot_onegeno_test


# function to retrieve only the dataframe that all qtl2 input dataframes were based on for testing purposes, and/or the haplotype columns
# @param summary_df (data.frame) - table as shown on HaploQA sample website, same as output of sample_summary_scrape
# @param data_dir (string) - filepath to where individual files (outputs of sample_individual_scrape) are stored
#
# @return df_all (data.frame) - dataframe with haplotype columns
get_haplotypes <- function(summary_df, data_dir) {
  # all txt file from directory
  data_files <- dir(data_dir, pattern = '\\.txt$', full.names = TRUE)
  ### combine all files
  df_raw <- rbindlist(lapply(data_files, read_sample_txt)) # save Y chrom here for gender plotting
  
  # clean up summary df
  sum_df <- summary_df %>% 
    select(ID, Sex, `% Het. Calls`, `% No Call`, Platform) %>% rename(sample_id = ID)
  sum_df$`% Het. Calls` <- as.numeric(gsub("%", "", sum_df$`% Het. Calls`))
  sum_df$`% No Call` <- as.numeric(gsub("%", "", sum_df$`% No Call`))
  
  # eliminate those that has no call > 10%
  df_all <- merge(df_raw, sum_df, by = 'sample_id') %>% filter(`% No Call` < 10) %>% filter(!chromosome %in% c('0', 'Y', 'M'))
  
  return(df_all)
}

# function to sort chromosomes into custom order
# @param df (data.frame) - selected dataframe to be sorted
# @param sort_order (list) - order to sort the chromosomes
#
# @return df (data.frame) - dataframe with chromosomes sorted in the selected order
sort_chr <- function(df, sort_order) {
  df$chr <- factor(df$chr, levels = sort_order)
  df <- df[order(df$chr),] 
  return(df)
}


# function to add cross attributes to an object
# @param crosstype (string) - type of cross model
# @param num_chr (string) - chromosome to be mapped
# @param chromosomes (string) - chromosomes present in the model
# @param df (data.frame/list) - object to add cross attributes to
# @param x_only (TRUE/FALSE) - only change the is_x_chr attribute
#
# @return df (data.frame/list) - object with attributes added
make_cross <- function(crosstype, num_chr, chromosomes, df, x_only) {
  attr_chr <- setNames(num_chr == 'X', chromosomes)
  if(x_only == T) {
    attr(df, "is_x_chr") <- attr_chr
  } else {
    attr(df, "crosstype") <- crosstype
    attr(df, "alleles") <- c("A", "B", "C", "D", "E", "F", "G", "H")
    attr_chr <- setNames(num_chr == 'X', chromosomes)
    attr(df, "is_x_chr") <- attr_chr
    attr(df, "class") <- c('viterbi', 'list')
  }
  return(df)
}

# function to pull only the inputs needed for shiny app
# @param sample_type (string) - type of haploqa sample
# @param rds_dir (string) - directory to read rds files from
#
# @return phased_geno_comp_list (list of dataframes) - list containing all objects needed for shiny
shiny_viz_input <- function(sample_type, rds_dir) {
  phased_geno_comp_list <- list()
  map <- readRDS(paste0(rds_dir, '/map_', sample_type, '.rds')) # should share the same map
  ph_geno_haplo <- readRDS(paste0(rds_dir, '/ph_geno_haploqa_', sample_type, '.rds'))
  ph_geno_qtl2 <- readRDS(paste0(rds_dir, '/ph_geno_', sample_type, '.rds'))
  phased_geno_comp_list[['haplo']] <- ph_geno_haplo
  phased_geno_comp_list[['qtl2']] <- ph_geno_qtl2
  phased_geno_comp_list[['map']] <- map
  
  return(phased_geno_comp_list)
}



# function to get only the raw genotype data for the genotype comparison process on one chromosome
# @param sample_type (string) - type of sample, such as CC, DO, F2, etc.
# @param chromosome (string) - chromosome to reteieve raw genotype data for
#
# @return raw_geno (data.frame) - dataframe of raw genotype data for one chromosome
get_raw_geno <- function(sample_type, chromosome) {
  ### config file
  config <- fread(paste0(root, '/annotations_config.csv'))
  
  ### Environment
  config_sample <- config[config$array_type == sample_type]
  data_dir <- config_sample$data_dir
  # data output directory
  data_dir <- file.path(root, config_sample$data_dir)
  
  map_df <- fread(file.path(root, config_sample$qtl2_dir, 'test_gmap.csv'))
  
  # annotation file
  annot_file <- config_sample$annot_file
  
  ## summary file
  summary_df_fp <- paste0(data_dir, '/', sample_type, '_summary.csv')
  summary_df <- fread(summary_df_fp)
  
  dir <- "https://raw.githubusercontent.com/kbroman/MUGAarrays/main/UWisc/" # always the same
  annot_df <- read.csv(paste0(dir, annot_file))
  # encode_genome function on geno and founder from qtl2 convert package - args: matrix of genotypes, matrixs of two alleles
  annot_encode_df <- annot_df %>% select(marker, chr, snp) %>% # select allele columns
    separate(snp, c("allele1", "allele2"), sep=cumsum(c(1)))
  
  # all txt file from directory
  data_files <- dir(data_dir, pattern = '\\.txt$', full.names = TRUE)
  ### combine all files
  df_raw <- rbindlist(lapply(data_files, read_sample_txt))  # save Y chrom here for gender plotting
  
  # clean up summary df
  sum_df <- summary_df %>% 
    select(ID, Sex, `% Het. Calls`, `% No Call`, Platform) %>% rename(sample_id = ID)
  sum_df$`% Het. Calls` <- as.numeric(gsub("%", "", sum_df$`% Het. Calls`))
  sum_df$`% No Call` <- as.numeric(gsub("%", "", sum_df$`% No Call`))
  
  # eliminate those that has no call > 10%
  ### use checkifnot to make sure the markers/SNP are in the annotation file (annotation file is the boss)
  df_all <- merge(df_raw, sum_df, by = 'sample_id', sort = F) %>% filter(`% No Call` < 10) %>% filter(!chromosome %in% c('0', 'Y', 'M')) #%>% filter(!sample_id %in% exclude_list) 
  
  geno_sub <- df_all %>% select(sample_id, snp_id, allele1_fwd, allele2_fwd) %>%
    mutate(gene_exp = paste(allele1_fwd, allele2_fwd, sep = '')) %>%
    select(sample_id, snp_id, gene_exp) %>% rename(marker = snp_id)
  df_geno <- dcast(geno_sub, marker~sample_id, value.var="gene_exp")
  raw_geno <- merge(geno_sub, map_df, by = 'marker', sort = F) %>% filter(chr == chromosome) %>% arrange(pos)
  
  return(raw_geno)
  
}


# revised version of plot_onegeno function from qtl2 package
# instead of plotting one geno per chromosome, this function plots haploqa and qtl2 at the same time and put the genotypes side by side for each chr
# @param geno (data.frame) - genotype data of first model to plot (qtl2 or haploqa)
# @param geno1 (data.frame) - genotype data of the other model to plot (qtl2 or haploqa, different from geno above)
# @param map (named list) - genetic map, should be the same for both models
# @param ind (int) - the individual to plot for
#
# @output (ggplot object) - geno plot for the selected models
plot_onegeno_test <- function(geno, geno1, map, ind=1, chr=NULL,
                              col=NULL, na_col="white",
                              swap_axes=FALSE,
                              border="black", shift=FALSE,
                              chrwidth=0.5, ...) {
  
  # ignore class of geno object
  geno <- unclass(geno)
  geno1 <- unclass(geno1)
  
  # drop all but the target individual
  if(length(ind)>1) {
    ind <- ind[1]
    warning("Only using the first individual")
  }
  
  # separate mom and dad
  for(i in seq_along(geno)) {
    geno[[i]] <- geno[[i]][ind,,,drop=FALSE]
    geno1[[i]] <- geno1[[i]][ind,,,drop=FALSE]
    
  }
  
  # always start at 0
  map <- lapply(map, function(a) a-min(a,na.rm=TRUE))
  
  plot_onegeno_internal <-
    function(geno, map, col=NULL, na_col="white",
             swap_axes=FALSE,
             border="black", bgcolor="gray90",
             chrwidth=0.5,
             xlab=NULL, ylab=NULL,
             xlim=NULL, ylim=NULL, las=1, xaxs=NULL, yaxs=NULL,
             mgp.x=c(2.6,0.5,0), mgp.y=c(2.6,0.5,0), mgp=NULL,
             hlines=NULL, hlines_col="white", hlines_lwd=1, hlines_lty=1,
             vlines=NULL, vlines_col="gray80", vlines_lwd=1, vlines_lty=1,
             ...)
    {
      dots <- list(...)
      
      nchr <- length(map)
      
      # margin parameters
      if(!is.null(mgp)) mgp.x <- mgp.y <- mgp
      
      if(is.null(xlab)) xlab <- "Chromosome"
      if(is.null(ylab)) ylab <- "Position"
      
      if(is.null(xlim)) xlim <- c(0.5, 52)
      if(is.null(ylim)) ylim <- rev(range(unlist(map), na.rm=TRUE))
      
      if(is.null(hlines)) hlines <- pretty(ylim)
      if(is.null(vlines)) vlines <- seq_len(52)
      
      if(is.null(xaxs)) xaxs <- "i"
      if(is.null(yaxs)) yaxs <- "r"
      
      # blank canvas
      plot(0, 0, type="n", xlab="", ylab="",
           xaxs=xaxs, yaxs=yaxs,
           xaxt="n", yaxt="n", xlim=xlim, ylim=ylim)#), ...)
      
      u <- par("usr")
      if(!is.null(bgcolor))
        rect(u[1], u[3], u[2], u[4], col=bgcolor, border=NA)
      
      # include axis labels?
      if(is.null(dots$xaxt)) dots$xaxt <- par("xaxt")
      if(is.null(dots$yaxt)) dots$yaxt <- par("yaxt")
      
      # add x axis unless par(xaxt="n")
      if(dots$xaxt != "n") {
        odd <- seq(1, nchr, by=2)
        axis(side=1, at=(odd*2.5), names(map)[odd],
             mgp=mgp.x, las=las, tick=FALSE)
        if(nchr > 1) {
          even <- seq(2, nchr, by=2)
          axis(side=1, at=(even*2.5), names(map)[even],
               mgp=mgp.x, las=las, tick=FALSE)
        }
      }
      # add y axis unless par(yaxt="n")
      if(dots$yaxt != "n") {
        axis(side=2, at=pretty(ylim), mgp=mgp.y, las=las, tick=FALSE)
      }
      
      # grid lines
      if(!(length(vlines)==1 && is.na(vlines))) {
        abline(v=vlines, col=vlines_col, lwd=vlines_lwd, lty=vlines_lty)
      }
      if(!(length(hlines)==1 && is.na(hlines))) {
        abline(h=hlines, col=hlines_col, lwd=hlines_lwd, lty=hlines_lty)
      }
      
      # x and y axis labels
      title(xlab=xlab, mgp=mgp.x)
      title(ylab=ylab, mgp=mgp.y)
      
      max_geno <- max(unlist(geno), na.rm=TRUE)
      if(is.null(col)) {
        if(max_geno <= 8) {
          col <- qtl2::CCcolors
        }
        else {
          warning("With ", max_geno, " genotypes, you need to provide the vector of colors; recycling some")
          col <- rep(qtl2::CCcolors, max_geno)
        }
      }
      else if(max_geno > length(col)) {
        warning("not enough colors; recycling them")
        col <- rep(col, max_geno)
      }
      ### chromosomes 1-19
      inc <- 0
      
      for(i in seq_len(nchr-1)) {
        #i <- 19
        g <- geno[[i]]
        g1 <- geno1[[i]]
        
        # if completely missing the second chr but not the first, treat as if we have just the one
        #   (this is a kludge to deal with males on X chr;
        #    really should use is_x_chr and is_female but we don't have it)
        this_chrwidth <- chrwidth
        if(!is.matrix(g) && !all(is.na(g[,,1])) && all(is.na(g[,,2]))) { # if g needs to be converted, g1 should be too
          g <- rbind(g[,,1]) # make it a row matrix
          g1 <- rbind(g1[,,1])
          this_chrwidth <- this_chrwidth/2
        }
        ### dataframe - g
        # rectangle shape
        x_lower_left <- i+0.75+inc
        x_higher_left <- (i+1.25+inc)-(this_chrwidth/2)
        x_lower_right <- i+1.25+inc
        x_higher_right <- (i+0.75+inc)+(this_chrwidth/2)
        
        rect(x_higher_left, min(map[[i]], na.rm=TRUE),
             x_lower_left, max(map[[i]], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[[i]], na.rm=TRUE),
             x_lower_right, max(map[[i]], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        
        # add geno colors
        addgenorect(g[1,,1], map[[i]], x_higher_left,x_lower_left,
                    col=col, swap_axes=swap_axes)
        addgenorect(g[1,,2], map[[i]], x_higher_right, x_lower_right,
                    col=col, swap_axes=swap_axes)
        # borders
        rect(x_higher_left, min(map[[i]], na.rm=TRUE),
             x_lower_left, max(map[[i]], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[[i]], na.rm=TRUE),
             x_lower_right, max(map[[i]], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        
        ### dataframe - g1
        # rectangle shape
        x_lower_left <- x_lower_left + 1
        x_higher_left <- x_higher_left + 1
        x_lower_right <-x_lower_right + 1
        x_higher_right <- x_higher_right + 1
        
        rect(x_higher_left, min(map[[i]], na.rm=TRUE),
             x_lower_left, max(map[[i]], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[[i]], na.rm=TRUE),
             x_lower_right, max(map[[i]], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        
        # add geno colors
        addgenorect(g1[1,,1], map[[i]], x_higher_left, x_lower_left,
                    col=col, swap_axes=swap_axes)
        addgenorect(g1[1,,2], map[[i]], x_higher_right, x_lower_right,
                    col=col, swap_axes=swap_axes)
        # borders
        rect(x_higher_left, min(map[[i]], na.rm=TRUE),
             x_lower_left, max(map[[i]], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[[i]], na.rm=TRUE),
             x_lower_right, max(map[[i]], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        inc <- inc + 1.5
      }
      
      ## plot x individually
      g <- geno[['X']]
      g1 <- geno1[['X']]
      if(!is.matrix(g) && !all(is.na(g[,,1])) && all(is.na(g[,,2]))) { # if g needs to be converted, g1 should be too
        g <- rbind(g[,,1]) # make it a row matrix
        g1 <- rbind(g1[,,1])
        this_chrwidth <- this_chrwidth/2
      }
      
      inc <- inc + 1
      
      if(is.matrix(g)) {
        ### dataframe - g
        x_lower_left <- i+0.75+inc
        x_higher_left <- (i+1.25+inc)-(this_chrwidth/2)
        x_lower_right <- i+1.25+inc
        x_higher_right <- (i+0.75+inc)+(this_chrwidth/2)
        
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        
        addgenorect(g[1,], map[['X']], x_higher_left, x_lower_left,
                    col=col, swap_axes=swap_axes)
        
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        
        ### dataframe - g1
        x_lower_left <- x_lower_left + 1
        x_higher_left <- x_higher_left + 1
        x_lower_right <- x_lower_right + 1
        x_higher_right <- x_higher_right + 1
        
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        
        addgenorect(g1[1,], map[['X']], x_higher_left, x_lower_left,
                    col=col, swap_axes=swap_axes)
        
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        
        
      } else {
        # rectangle shape
        x_lower_left <- i+0.75+inc
        x_higher_left <- (i+1.25+inc)-(this_chrwidth/2)
        x_lower_right <- i+1.25+inc
        x_higher_right <- (i+0.75+inc)+(this_chrwidth/2)
        
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[['X']], na.rm=TRUE),
             x_lower_right, max(map[['X']], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        
        # add geno colors
        addgenorect(g[1,,1], map[['X']], x_higher_left,x_lower_left,
                    col=col, swap_axes=swap_axes)
        addgenorect(g[1,,2], map[['X']], x_higher_right, x_lower_right,
                    col=col, swap_axes=swap_axes)
        # borders
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[['X']], na.rm=TRUE),
             x_lower_right, max(map[['X']], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        
        ### dataframe - g1
        # rectangle shape
        x_lower_left <- x_lower_left + 1
        x_higher_left <- x_higher_left + 1
        x_lower_right <-x_lower_right + 1
        x_higher_right <- x_higher_right + 1
        
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[['X']], na.rm=TRUE),
             x_lower_right, max(map[['X']], na.rm=TRUE),
             col=na_col, border=border, lend=1, ljoin=1)
        
        # add geno colors
        addgenorect(g1[1,,1], map[['X']], x_higher_left, x_lower_left,
                    col=col, swap_axes=swap_axes)
        addgenorect(g1[1,,2], map[['X']], x_higher_right, x_lower_right,
                    col=col, swap_axes=swap_axes)
        # borders
        rect(x_higher_left, min(map[['X']], na.rm=TRUE),
             x_lower_left, max(map[['X']], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
        rect(x_higher_right, min(map[['X']], na.rm=TRUE),
             x_lower_right, max(map[['X']], na.rm=TRUE),
             col=NULL, border=border, lend=1, ljoin=1)
      }
      box()
    }
  
  plot_onegeno_internal(geno, map, col=col, na_col=na_col,
                        swap_axes=swap_axes, border=border,
                        chrwidth=chrwidth, ...)
  
}

# add rectangles for the genotypes
addgenorect <- function(geno, map, x1, x2, col, swap_axes=FALSE) {
  intervals <- geno2intervals(geno, map)
  if(is.null(intervals) || nrow(intervals) < 1) return(NULL)
  
  for(i in seq_len(nrow(intervals))) {
    if(swap_axes) {
      rect(intervals[i,1], x1,
           intervals[i,2], x2,
           col=col[intervals[i,3]],
           border=NA, lend=1, ljoin=1)
    } else{
      rect(x1, intervals[i,1],
           x2, intervals[i,2],
           col=col[intervals[i,3]],
           border=NA, lend=1, ljoin=1)
    }
  }
}


# calculate the intervals between genos, taken from Karl's qtl2 code repo
# used in plot_onegeno_test
# @param geno (data.frame) - genotype data of first model to plot (qtl2 or haploqa)
# @param map (named list) - genetic map, should be the same for both models
#
# @output (dataframe object) - a data frame with calculated geno intervals
geno2intervals <- function(geno, map) {
  if(all(is.na(geno))) return(NULL)
  
  stopifnot(length(geno) == length(map))
  
  # drop missing values
  map <- map[!is.na(geno)]
  geno <- geno[!is.na(geno)]
  
  d <- diff(geno)
  xo_int <- which(d != 0)
  
  data.frame(lo=map[c(1,xo_int+1)],
             hi=map[c(xo_int, length(map))],
             geno=geno[c(xo_int, length(map))])

}



geno_align <- function(df, n_founders) {
  founder_codes <- LETTERS[seq(1,n_founders)]
  col_diff <- setdiff(founder_codes, colnames(df))
  dum_array <- array(0, dim = c(nrow(df), length(col_diff)))
  colnames(dum_array) <- col_diff
  df <- cbind(df, dum_array)
  df <- df[ ,founder_codes]
  
  return(df)
}


# flattens the genoprob object into one data frame
# data frame should have genoprob value columns with markers and chromosome metadata
# @param sample_result (object) - should be the output of the pipelines, an object that contains all relevant results
# @param sample_type (string) - type of sample, (CC, DO, F2, etc.)
#
# @output genoprob_all (data.frame) - flattened genoprob dataframe
get_genoprob_example <- function(sample_result, sample_type) {
  num_chr <- c((1:19),"X")
  
  genoprob_all <- list()
  for (chromosome in num_chr) {
    print(chromosome)
    chr_genoprob <- sample_result[['pr']][[chromosome]]
    
    chr_list <- list()
    for(i in dimnames(chr_genoprob)[[3]]) {
      m <- chr_genoprob[,,i] %>% as.data.frame()
      m$marker <- i
      m$sample_name <- rownames(m)
      chr_list[[i]] <- m
    }
    chr_df <- rbindlist(chr_list)
    genoprob_all[[chromosome]] <- chr_df
  }
  genoprob_all <- rbindlist(genoprob_all, fill = TRUE)
  return(genoprob_all)
}


ind_geno_comp <- function(sample_results, sample_name, sample_type) {
  num_chr <- c((1:19),"X")
  founder_lookup_table <- fread(file.path(root, 'founder_lookup_table.csv'))
  founder_codes_dict <- setNames(founder_lookup_table$founder_codes, founder_lookup_table$founder_id)
  
  geno_codes <- colnames(sample_results[['pr']]$X)
  founder_all_rev_lookup <- setNames(geno_codes, seq(1, length(geno_codes)))
  
  config <- fread(paste0(root, '/annotations_config.csv'))
  config_sample <- config[config$array_type == sample_type]
  qtl2_dir <- file.path(root, config_sample$qtl2_dir)
  ### sample folder
  sample_dir <- file.path(qtl2_dir, sample_name)
  
  df_cov <- fread(file.path(sample_dir, 'test_covar.csv')) %>% rename('sample_id' = 'id')
  data_dir <- paste0(root, '/', config_sample$data_dir)
  summary_df <- fread(paste0(sample_dir, '/', sample_type, '_', sample_name, '_summary.csv'))
  
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

