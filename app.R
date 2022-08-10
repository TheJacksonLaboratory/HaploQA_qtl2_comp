library(rstudioapi)
library(shiny)


root <- dirname(getSourceEditorContext()$path)
source(paste0(root,"/input_data_prep_functions.R"))

rds_dir <- file.path(results_dir, 'RDS')

total_ind <- 277 # different for CC and DO - fix it

phased_geno_comp_cc <- shiny_viz_input('CC', rds_dir) # contains phased geno of qtl2 AND haploqa
# contains: qtl2, haplo, map
phased_geno_comp_do <- shiny_viz_input('DO', rds_dir)

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
        plot_onegeno(phased_geno_comp_do[['qtl2']], phased_geno_comp_do[['map']], ind = individual, shift = TRUE, main = paste0('qtl2 - Geno-wide genotypes of individual ', individual)) 
      })
      
      output$genoplot_haploqa <- renderPlot({
        plot_onegeno(phased_geno_comp_do[['haplo']], phased_geno_comp_do[['map']], ind = individual, shift = TRUE, main = paste0('haploqa - Geno-wide genotypes of individual ', individual)) 
      })
    } else if (input$select=="Collaborative Cross") {
      individual <- as.numeric(input$individual)
      
      output$genoplot_qtl2 <- renderPlot({
        plot_onegeno(phased_geno_comp_cc[['qtl2']], phased_geno_comp_cc[['map']], ind = individual, shift = TRUE, main = paste0('qtl2 - Geno-wide genotypes of individual ', individual)) 
      })
      
      output$genoplot_haploqa <- renderPlot({
        plot_onegeno(phased_geno_comp_cc[['haplo']], phased_geno_comp_cc[['map']], ind = individual, shift = TRUE, main = paste0('haploqa - Geno-wide genotypes of individual ', individual)) 
      })
    }
  })
}

shinyApp(ui, server)
