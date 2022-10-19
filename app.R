library(rstudioapi)
library(shiny)


root <- dirname(getSourceEditorContext()$path)
source(paste0(root,"/input_data_prep_functions.R"))

## results directory
results_dir <- file.path(root, 'results')
rds_dir <- file.path(results_dir, 'RDS')

total_ind <- 277 # different for CC and DO - fix it

phased_geno_comp_cc <- shiny_viz_input('CC', rds_dir) # contains phased geno of qtl2 AND haploqa
# contains: qtl2, haplo, map
phased_geno_comp_do <- shiny_viz_input('DO', rds_dir)

phased_geno_comp_bxd <- shiny_viz_input('BXD', rds_dir)

phased_geno_comp_f2 <- shiny_viz_input('F2', rds_dir)


minimuga_dir <- '/Users/linc/Documents/GitHub/HaploQA_qtl2_comp/minimuga_qtl2_genail/'


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
                   selectInput("individual1", "Individual", choices = seq(1:222), selected = 1)),
                 conditionalPanel(
                   condition = "input.select == 'BXD F3'",
                   selectInput("individual2", "Individual", choices = seq(1:71), selected = 1)),
                conditionalPanel(
                  condition = "input.select == 'F2'",
                  selectInput("individual3", "Individual", choices = seq(1:34), selected = 1)),
                conditionalPanel(
                  condition = "input.select == 'MiniMUGA'",
                  selectInput("individual4", "Individual", choices = seq(1:25), selected = 1))
    ),

    mainPanel(plotOutput('genoplot', width = '100%', height = '800px'), textOutput("text1"), textOutput('pct_qtl2_haploqa'))
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
      
      individual <- as.numeric(input$individual)
      print(input$individual)
      
      output$genoplot <- renderPlot({
        plot_onegeno_test(phased_geno_comp_do[['qtl2']], phased_geno_comp_do[['haplo']], phased_geno_comp_do[['map']], ind = individual, 
                          shift = TRUE, main = paste0('Geno-wide genotypes of individual ', individual),
                          sub = 'Left - qtl2, Right - HaploQA')
        legend("bottomright", legend = c(names(qtl2::CCcolors)), title = 'Founder Colors', fill = qtl2::CCcolors)
      })
      
      #output$pct_qtl2_haploqa <- ''
      output$text1 <- NULL
      
    } else if (input$select=="Collaborative Cross") {
      individual <- as.numeric(input$individual1)
      print(input$individual1)
      output$genoplot <- renderPlot({
        plot_onegeno_test(phased_geno_comp_cc[['qtl2']], phased_geno_comp_cc[['haplo']], phased_geno_comp_cc[['map']], ind = individual, 
                          shift = TRUE, main = paste0('qtl2 - Geno-wide genotypes of individual ', individual),
                          sub = 'Left - qtl2, Right - HaploQA') 
        legend("bottomright", legend = c(names(qtl2::CCcolors)), title = 'Founder Colors', fill = qtl2::CCcolors)
      })
      
      output$text1 <- NULL
      
    } else if (input$select=="BXD F3") {
      individual <- as.numeric(input$individual2)
      print(input$individual2)
      output$genoplot <- renderPlot({
        plot_onegeno_test(phased_geno_comp_bxd[['qtl2']], phased_geno_comp_bxd[['haplo']], phased_geno_comp_bxd[['map']], ind = individual, 
                          shift = TRUE, main = paste0('qtl2 - Geno-wide genotypes of individual ', individual),
                          sub = 'Left - qtl2, Right - HaploQA') 
        #legend("bottomright", legend = c(names(qtl2::CCcolors)), title = 'Founder Colors', fill = qtl2::CCcolors)
      })
      
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
      output$text1 <- renderPrint({
        ind_name
      })
      
    }
  })
}

shinyApp(ui, server)
