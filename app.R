library(rstudioapi)
library(shiny)


root <- dirname(getSourceEditorContext()$path)
source(paste0(root,"/input_data_prep_functions.R"))

## results directory
results_dir <- file.path(root, 'results')
rds_dir <- file.path(results_dir, 'RDS')

# phased geno files
phased_geno_comp_cc <- shiny_viz_input('CC', rds_dir) # contains phased geno of qtl2 AND haploqa
# contains: qtl2, haplo, map
phased_geno_comp_do <- shiny_viz_input('DO', rds_dir)
phased_geno_comp_bxd <- shiny_viz_input('BXD', rds_dir)
phased_geno_comp_f2 <- shiny_viz_input('F2', rds_dir)

# comparison files
truth_comp_fp <- file.path(results_dir, 'shiny_pct_csvs')

# minimuga data
minimuga_dir <- '/Users/linc/Documents/GitHub/HaploQA_qtl2_comp/minimuga_qtl2_genail/'

### add a legend that shows which side is which method
### to launch shiny app
ui <- fluidPage(
  titlePanel("Geno-wide Genotype comparison - qtl2 vs. HaploQA"),
  selectInput(inputId = "select", label = "Select model",
              choices = c("Diversity Outbred","Collaborative Cross", "BXD F3", "F2", "MiniMUGA"), plotOutput(outputId = "plot")),
  sidebarLayout(
    sidebarPanel(width = 3,
                 #htmlOutput("selectUI"),
                 #selectInput("individual", "Individual", ''),
                 conditionalPanel(
                   condition = "input.select == 'Diversity Outbred'",
                   selectInput("individual", "Individual", choices = seq(1:277), selected = 1)),
                 conditionalPanel(
                   condition = "input.select == 'Collaborative Cross'",
                   selectInput("individual1", "Individual", choices = seq(1:218), selected = 1)),
                 conditionalPanel(
                   condition = "input.select == 'BXD F3'",
                   selectInput("individual2", "Individual", choices = seq(1:71), selected = 1)),
                 conditionalPanel(
                   condition = "input.select == 'F2'",
                   selectInput("individual3", "Individual", choices = seq(1:34), selected = 1)),
                 conditionalPanel(
                   condition = "input.select == 'MiniMUGA'",
                   selectInput("individual4", "Individual", choices = seq(1:25), selected = 1)),
                 textOutput("text1")
    ),
    
    mainPanel(fluidRow(align = "center", htmlOutput("model_descipt")), 
              plotOutput('genoplot', width = '100%', height = '800px'), 
              div(style = "height:20px"),
              textOutput("col_descript"), 
              htmlOutput("col_descript1"), 
              htmlOutput("col_descript2"), 
              htmlOutput("col_descript3"), 
              div(style = "height:20px"),
              tableOutput('pct_qtl2_haploqa'))
  )
)

server <- function(input, output, session) {
  observe({
    print(input$select)
    if(input$select=="Diversity Outbred"){
      #individual <- as.numeric(input$individual)
      # print(individual)
      #output$selectUI <- renderUI({ 
      #selectInput("individual", "Select Individual", choices = seq(1:277))
      #})
      do_comp_truth <- fread(file.path(truth_comp_fp, 'do_truth_comp.csv'))
      
      individual <- as.numeric(input$individual)
      print(input$individual)
      
      output$genoplot <- renderPlot({
        plot_onegeno_test(phased_geno_comp_do[['qtl2']], phased_geno_comp_do[['haplo']], phased_geno_comp_do[['map']], ind = individual, 
                          shift = TRUE, main = paste0('Geno-wide genotypes of individual ', individual),
                          sub = 'Left - qtl2, Right - HaploQA')
        legend("bottomright", legend = c(names(qtl2::CCcolors)), title = 'Founder Colors', fill = qtl2::CCcolors)
      })
      
      # get the according comparison from csv
      sample_comp <- do_comp_truth[individual,]
      
      output$model_descipt <- renderText('</B>Left - qtl2, Right - HaploQA</B>')
      output$col_descript <- renderText('Column Descriptions')
      output$col_descript1 <- renderText('<B>qtl2_pct_diff</B>: percentage of markers that are different between the truth model and qtl2')
      output$col_descript2 <- renderText('<B>haploqa_pct_diff</B>: percentage of markers that are different between the truth model and haploqa')
      output$col_descript3 <- renderText('<B>haploqa_qtl2_pct_diff</B>: percentage of markers that are different between qtl2 and haploqa')
      
      output$pct_qtl2_haploqa <- renderTable({sample_comp})
      output$text1 <- NULL
      
    } else if (input$select=="Collaborative Cross") {
      cc_comp_truth <- fread(file.path(truth_comp_fp, 'cc_truth_comp.csv'))
      
      individual <- as.numeric(input$individual1)
      print(input$individual1)
      output$genoplot <- renderPlot({
        plot_onegeno_test(phased_geno_comp_cc[['qtl2']], phased_geno_comp_cc[['haplo']], phased_geno_comp_cc[['map']], ind = individual, 
                          shift = TRUE, main = paste0('qtl2 - Geno-wide genotypes of individual ', individual),
                          sub = 'Left - qtl2, Right - HaploQA') 
        legend("bottomright", legend = c(names(qtl2::CCcolors)), title = 'Founder Colors', fill = qtl2::CCcolors)
      })
      
      sample_comp <- cc_comp_truth[individual,]
      
      # text displays
      output$col_descript <- renderText('Column Descriptions')
      output$col_descript1 <- renderText('qtl2_pct_diff: percentage of markers that are different between the truth model and qtl2')
      output$col_descript2 <- renderText('haploqa_pct_diff: percentage of markers that are different between the truth model and haploqa')
      output$col_descript3 <- renderText('haploqa_qtl2_pct_diff: percentage of markers that are different between qtl2 and haploqa')
      
      output$pct_qtl2_haploqa <- renderTable({sample_comp})
      output$text1 <- NULL
      
    } else if (input$select=="BXD F3") {
      bxd_comp_truth <- fread(file.path(truth_comp_fp, 'bxd_truth_comp.csv'))
      individual <- as.numeric(input$individual2)
      print(input$individual2)
      output$genoplot <- renderPlot({
        plot_onegeno_test(phased_geno_comp_bxd[['qtl2']], phased_geno_comp_bxd[['haplo']], phased_geno_comp_bxd[['map']], ind = individual, 
                          shift = TRUE, main = paste0('qtl2 - Geno-wide genotypes of individual ', individual),
                          sub = 'Left - qtl2, Right - HaploQA') 
        #legend("bottomright", legend = c(names(qtl2::CCcolors)), title = 'Founder Colors', fill = qtl2::CCcolors)
      })
      
      sample_comp <- bxd_comp_truth[individual,]
      
      output$col_descript <- renderText('Column Descriptions')
      output$col_descript1 <- renderText('qtl2_pct_diff: percentage of markers that are different between the truth model and qtl2')
      output$col_descript2 <- renderText('haploqa_pct_diff: percentage of markers that are different between the truth model and haploqa')
      output$col_descript3 <- renderText('haploqa_qtl2_pct_diff: percentage of markers that are different between qtl2 and haploqa')
      
      output$pct_qtl2_haploqa <- renderTable({sample_comp})
      output$text1 <- NULL
      
    } else if (input$select=="F2") {
      individual <- as.numeric(input$individual3)
      print(input$individual3)
      output$genoplot <- renderPlot({
        plot_onegeno_test(phased_geno_comp_f2[['qtl2']], phased_geno_comp_f2[['haplo']], phased_geno_comp_f2[['map']], ind = individual, 
                          shift = TRUE, main = paste0('qtl2 - Geno-wide genotypes of individual ', individual),
                          sub = 'Left - qtl2, Right - HaploQA') 
        # legend("bottomright", legend = c(names(qtl2::CCcolors)), title = 'Founder Colors', fill = qtl2::CCcolors)
      })
      
      output$col_descript <- renderText('Column Descriptions')
      output$col_descript1 <- renderText('qtl2_pct_diff: percentage of markers that are different between the truth model and qtl2')
      output$col_descript2 <- renderText('haploqa_pct_diff: percentage of markers that are different between the truth model and haploqa')
      output$col_descript3 <- renderText('haploqa_qtl2_pct_diff: percentage of markers that are different between qtl2 and haploqa')
      
      output$pct_qtl2_haploqa <- NULL
      output$text1 <- NULL
      
    } else if (input$select=="MiniMUGA") {
      individual <- as.numeric(input$individual4)
      print(input$individual4)
      ind_name <- list.files(minimuga_dir)[individual]
      sample_mini_dir <- file.path(minimuga_dir, ind_name)
      #rds_fps <- list.files(path = sample_mini_dir, pattern = "\\.rds$")
      mini_qtl2 <- readRDS(file.path(sample_mini_dir, 'ph_geno.rds'))
      mini_haploqa <- readRDS(file.path(sample_mini_dir, 'ph_geno_haploqa.rds'))
      map_ind <- readRDS(file.path(sample_mini_dir, 'map.rds'))
      output$genoplot <- renderPlot({
        plot_onegeno_test(mini_qtl2, mini_haploqa, map_ind, ind = 1, 
                          shift = TRUE, main = paste0('qtl2 - Geno-wide genotypes of individual ', ind_name),
                          sub = 'Left - qtl2, Right - HaploQA') 
        #legend("bottomright", legend = c(names(qtl2::CCcolors)), title = 'Founder Colors', fill = qtl2::CCcolors)
      })
      
      output$col_descript <- NULL
      output$col_descript1 <- NULL
      output$col_descript2 <- NULL
      output$col_descript3 <- NULL
      
      output$pct_qtl2_haploqa <- NULL
      output$text1 <- renderText({
        ind_name
      })
      
    }
  })
}

shinyApp(ui, server)
