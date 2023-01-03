#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)

options(shiny.maxRequestSize=30*1024^2)

# source("../functions/functions.R")

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  require(tidyverse)
  require(plotly)



  # read files --------------------------------------------------------------
  # get files containing csv of CQ values
  combined_df <- reactive({
    if (is.null(input$cq_files))
      return(NULL)
    inFile <- bind_rows(lapply(input$cq_files$datapath, read_csv)) %>% filter(Content %>% str_detect("Unkn"))
    return(inFile)
  })

  # get files containing metadata
  metadata_df <- reactive({
    if (is.null(input$file_input_md))
      return(NULL)
    inFile <- read_csv(input$file_input_md$datapath)
    return(inFile)
  })

  # get files containing efficiency
  efficiency_df <- reactive({
    if (is.null(input$file_input_efficiency))
      return(NULL)
    inFile <- read_csv(input$file_input_efficiency$datapath)
    return(inFile)
  })


  # select control data ----------------------------------------------------------
  # Select control genes
  output$gene_select_1 <- renderUI({
    if(is.null(input$cq_files))
      return()

    target_primer <- combined_df() %>% select(Target) %>% unique()

    # Create the dropdown menu
    selectInput("selected_ctrl_gene_1", "Choose control gene 1",
                choices  = target_primer,
                selected = target_primer[1])
  })

  output$gene_select_2 <- renderUI({
    if(is.null(input$cq_files))
      return()

    target_primer <- combined_df() %>% select(Target) %>% unique()

    # Create the dropdown menu
    selectInput("selected_ctrl_gene_2", "Choose control gene 2",
                choices  = target_primer,
                selected = target_primer[1])
  })

  # select control sample and replicate
  output$control_sample <- renderUI({
    if(is.null(input$cq_files))
      return()

    target_sample <- combined_df() %>% select(Sample) %>% unique()

    # Create the dropdown menu
    selectInput("selected_control_sample", "Choose control sample",
                choices  = target_sample,
                selected = target_sample[1])
  })

  output$control_sample_replicate <- renderUI({
    if(is.null(input$cq_files))
      return()

    target_sample <- combined_df() %>% select(Sample) %>% unique()
    rep <- combined_df() %>% filter(Sample == input$selected_control_sample) %>% select(`Biological Set Name`)

    # Create the dropdown menu
    selectInput("selected_control_sample_replicate", "Choose control sample replicate",
                choices  = rep,
                selected = rep[1])
  })

  # select samples to analyze
  output$genes_select <- renderUI({
    if(is.null(input$cq_files))
      return()

    target_sample <- combined_df() %>% pull(Sample) %>% unique()

    # Create the checkboxes and select them all by default
    checkboxGroupInput("selected_samples", "Choose samples to analyze",
                       choices  = target_sample)
  })




  # perform calculations ----------------------------------------------------



  # output combined df
  output$contents <- renderTable({
    combined_df()
  }, align = "c", width = "100%")



  # plot --------------------------------------------------------------------

})
