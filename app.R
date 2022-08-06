library(shiny)
library(qtl2)
library(dplyr)
library(data.table)
library(rstudioapi)
library(stringi)

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

#### RDS functions
sample_types <- c('CC', 'DO') # right now supports CC, DO, GigaMUGA, MiniMUGA
rds_loop(sample_types, config, results_dir) # check if their RDS files exist


## attributes from cc_test
total_ind <- 277 # take this automatically from cc_test?
num_chr <- c((1:19),"X")

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
founder_all_lookup <- setNames(all_map_codes, founder_all_codes)
founder_all_rev_lookup <- setNames(founder_all_codes, all_map_codes)


# result
df_all <- get_haplotypes(summary_df, data_dir)
df_test <- df_all %>% select(sample_id, snp_id, haplotype1, haplotype2, chromosome) %>% unique()
df_test[,c(3,4)] <- as.data.frame(apply(df_test[,c(3,4)], 2, function(x) founder_haplo_lookup[x]))
df_test$haplotype <- paste(df_test$haplotype2, df_test$haplotype1, sep='')

haploqa_maxmarg <- df_test %>% select(sample_id, snp_id, chromosome, haplotype)

haploqa_maxmarg[,4] <- lapply(haploqa_maxmarg[,4], function(col) {
  rev_col = stri_reverse(col)
  ifelse(rev_col %in% founder_all_codes, rev_col, col)
})

haploqa_maxmarg[,4] <- as.data.frame(apply(haploqa_maxmarg[,c(4)], 2, function(x) founder_all_lookup[x]))


##### csv data files
### find crossovers
ginf_haploqa <- list()
map_haploqa <- list()
for (i in num_chr) {
  #i<- 1
  print(i)
  map_df <- df_gmap %>% filter(chr == i) %>% arrange(pos) %>% as.data.frame()
  rownames(map_df) <- map_df$marker
  col_order <- map_df$marker
  map_df <- map_df %>% select(-c(marker, chr)) %>% t()
  map_haploqa[[i]] <- map_df
  
  ginf_df <- haploqa_maxmarg %>% filter(chromosome == i) %>% dcast(sample_id~snp_id, value.var = 'haplotype') %>% as.data.frame()
  rownames(ginf_df) <- ginf_df$sample_id
  ginf_df <- ginf_df[,col_order]
  ginf_haploqa[[i]] <- ginf_df
}

# add attributes as shown in cross
attr(ginf_haploqa, "crosstype") <- 'genail8'
attr(ginf_haploqa, "alleles") <- c("A", "B", "C", "D", "E", "F", "G", "H")
attr_chr <- setNames(num_chr == 'X', num_chr)
attr(ginf_haploqa, "is_x_chr") <- attr_chr
attr(ginf_haploqa, "class") <- c('viterbi', 'list')
attr(map_haploqa, "is_x_chr") <- attr_chr
pos_h <- locate_xo(ginf_haploqa, map_haploqa)

x_loc_diff <- list()
for (chr in num_chr) {
  #chr <- '1'
  print(chr)
  pos_h_chr <- pos_h[[chr]]
  pos_chr <- pos[[chr]]
  nmar_h <- sapply(pos_h_chr, length)
  nmar <- sapply(pos_chr, length) # none has the same length...
  #x_loc_diff[[chr]] <- list()
  for (i in seq(1:length(pos_chr))) {
    i <- 275
    x <- unlist(pos_h_chr[i])
    y <- unlist(pos_chr[i])
    if (length(x) > length(y)) {
      df <- apply(abs(outer(x, y, '-')), 2, min)
    }
    if (length(x) < length(y)) {
      df <- apply(abs(outer(y, x, '-')), 2, min)
    }
    x_loc_diff[[chr]][[i]]<- df
  }
}

## genocode comparison matrix
comp_matrix <- list()
for (chr in num_chr) {
  ## qtl2
  chr <- 1
  map_chr <- as.data.frame(map[[chr]]) %>% rename('pos' = 'map[[chr]]')
  map_chr$marker <- rownames(map_chr)
  # rowname is sample id
  df <- as.data.frame(ginf[[chr]])
  df_chr_qtl2 <- apply(df, 2, function(x) founder_all_rev_lookup[x]) # column level
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

}


### percent genomic difference
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

###### Visualizations
### genotype plots
df_haplo <- df_test %>% select(sample_id, snp_id, chromosome, haplotype1, haplotype2) %>% rename(marker = snp_id, chr = chromosome)
df_haplo[,c(4,5)] <- as.data.frame(apply(df_haplo[,c(4,5)], 2, function(x) founder_reverse[x])) ## heterozygosity?

df_chr_all <- merge(df_haplo, df_gmap, by = 'marker') # to align with the markers and samples
df_chr_all <- df_chr_all %>% select(marker, sample_id, pos, chr.x, haplotype1, haplotype2) %>% rename(chr = chr.x, mom = haplotype1, dad = haplotype2)


df_haplo_geno <- list()
df_haplo_map <- list()
for (i in num_chr) {
  #i<-1
  print(i)
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


### to output into local pdfs
## genotype
#pdf(file = paste0(results_dir, '/geno_plots_comp_cc.pdf'), width=20, height=10, bg="white")
#for (ind in seq(1:total_ind)) {
#  print(ind)
#  ind <- 1
#  par(mfrow = c(1, 2))
#  plot_onegeno(ph_geno, map, ind = ind, shift = TRUE, main = paste0('qtl2 - Geno-wide genotypes of individual ', ind))  # one individual's geno_wide genotype
#  plot_onegeno(df_haplo_geno, df_haplo_map, ind = ind, shift = TRUE, main = paste0('haploqa - Geno-wide genotypes of individual ', ind))  # one individual's geno_wide genotype
#}
#dev.off()
#
## genoprob
#pdf(file = paste0(results_dir, '/test_plots.pdf'))
#for (ind in seq(1:total_ind)) {
#  print(paste0('individual: ',ind))
#  for (chr in num_chr) {
#    print(chr)
#    plot_genoprob(pr, cc_test$gmap, ind = ind, chr = chr, main = paste0('Genotype Probabilities of individual ', ind, ' at chromosome ', chr))
#    
#  }
#}
#dev.off()


### to launch shiny app
ui <- fluidPage(
  titlePanel("Geno-wide Genotype comparison - qtl2 vs. HaploQA"),
  selectInput(inputId = "select", label = "Select model",
              choices = c("Diversity Outbred","Collaborative Cross"), plotOutput(outputId = "plot")),
  sidebarLayout(
    sidebarPanel(width = 3,
                 selectInput("individual", "Individual", choices = seq(1:total_ind), selected = 1),
                 uiOutput("chromosome")
    ),
    mainPanel(fluidRow(column(width = 6, plotOutput('genoplot_qtl2', width = '100%', height = '600px')), 
                       column(width = 6, plotOutput('genoplot_haploqa', width = '100%', height = '600px'))) 
      ))
)

server <- function(input, output, session) {
  
  observe({
    print(input$select)
    if(input$select=="Diversity Outbred"){
    individual <- as.numeric(input$individual)
    print(individual)
    output$genoplot_qtl2 <- renderPlot({
      plot_onegeno(ph_geno, map, ind = individual, shift = TRUE, main = paste0('qtl2 - Geno-wide genotypes of individual ', individual)) 
    })
    
    output$genoplot_haploqa <- renderPlot({
      plot_onegeno(df_haplo_geno, df_haplo_map, ind = individual, shift = TRUE, main = paste0('haploqa - Geno-wide genotypes of individual ', individual)) 
    })
  } else if (input$select=="Collaborative Cross") {
    individual <- as.numeric(input$individual)
    
    output$genoplot_qtl2 <- renderPlot({
      plot_onegeno(ph_geno_cc, map, ind = individual, shift = TRUE, main = paste0('qtl2 - Geno-wide genotypes of individual ', individual)) 
    })
    
    output$genoplot_haploqa <- renderPlot({
      plot_onegeno(df_haplo_geno, df_haplo_map, ind = individual, shift = TRUE, main = paste0('haploqa - Geno-wide genotypes of individual ', individual)) 
    })
    }
  })
}

shinyApp(ui, server)


### save the output of pr (calc_genoprob) in csv files?
#for (chr in num_chr) {
#  print(chr)
#  x <- eval(parse(text = deparse(substitute(chr))))
#  write.csv(pr[[chr]], paste0(results_dir, '/prob_chr_', x, '.csv'), row.names = F)
#}

