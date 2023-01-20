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
    inFile <- bind_rows(lapply(input$cq_files$datapath, read_csv)) %>%
      filter(Content %>% str_detect("Unkn")) %>%
      set_names(make.names(names(.))) %>%
      select(Target, Content, Sample, Biological.Set.Name, Cq.Mean, Cq.Std..Dev)
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
    inFile[[2]] <- (inFile[[2]]/100)+1
    colnames(inFile) <- c("Target", "efficiency")
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
    if(is.null(input$selected_control_sample))
      return()

    target_sample <- combined_df() %>% select(Sample) %>% unique()
    rep <- combined_df() %>% filter(Sample == input$selected_control_sample) %>% select(`Biological.Set.Name`)

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
  # create calculated dataframe
  calculated_df_list <- reactive({
    if (is.null(input$cq_files))
      return(NULL)
    if (is.null(input$file_input_efficiency))
      return(NULL)

    df <- combined_df()

    # filter samples
    df <- df %>% filter(Sample %in% input$selected_samples)

    # extract controls
    control_1 <- df %>% filter(Target == input$selected_ctrl_gene_1) %>% distinct() %>%
      select(Sample, Biological.Set.Name, mean_control_1 = Cq.Mean)
    control_2 <- df %>% filter(Target == input$selected_ctrl_gene_2) %>% distinct() %>%
      select(Sample, Biological.Set.Name, mean_control_2 = Cq.Mean)

    # extract cq for control sample on control genes
    control_1_ctrsample <- control_1 %>% filter(Sample == input$selected_control_sample, Biological.Set.Name == input$selected_control_sample_replicate) %>%
      select(Biological.Set.Name, mean_control_1_sample = mean_control_1, -Biological.Set.Name)
    control_2_ctrsample <- control_2 %>% filter(Sample == input$selected_control_sample, Biological.Set.Name == input$selected_control_sample_replicate) %>%
      select(Biological.Set.Name, mean_control_2_sample = mean_control_2, -Biological.Set.Name)

    # extract sample to use as control
    df_bio_control <- df %>% filter(Sample == input$selected_control_sample, Biological.Set.Name == input$selected_control_sample_replicate) %>%
      distinct() %>%
      select(Target, Biological.Set.Name, cq_control_sample = Cq.Mean, -Biological.Set.Name)

    # remove samples used as control
    df_samples <- df %>% mutate(is.control = case_when(
      (Sample == input$selected_control_sample & Biological.Set.Name == input$selected_control_sample_replicate) ~ "control",
      TRUE ~ "sample"
    )) %>% filter(is.control == "sample") %>% select(-is.control) %>% distinct()

    # remove erroneous data
    df_samples <- df_samples %>% filter(Cq.Std..Dev != 0)
    df_bio_control <- df_bio_control %>% filter(cq_control_sample != 0)

    # combine dataframes to one
    df_samples <- df_samples %>% left_join(control_1) %>% left_join(control_2) %>% left_join(df_bio_control) %>% left_join(efficiency_df()) %>%
      mutate(mean_control_1_sample = control_1_ctrsample$mean_control_1_sample, mean_control_2_sample = control_2_ctrsample$mean_control_2_sample)

    # efficiency of control genes
    efficiency_ctrl_gene_1 <- efficiency_df() %>% filter_at(1, all_vars(. == input$selected_ctrl_gene_1)) %>% .[[2]]
    efficiency_ctrl_gene_2 <- efficiency_df() %>% filter_at(1, all_vars(. == input$selected_ctrl_gene_2)) %>% .[[2]]

    # calculate RE and SD
    calculated_df <- df_samples %>% mutate(RE = ((efficiency)^(cq_control_sample - Cq.Mean))/
                                             sqrt((efficiency_ctrl_gene_1)^(mean_control_1_sample-mean_control_1) * (efficiency_ctrl_gene_2)^(mean_control_2_sample-mean_control_2)))
    calculated_df_summary <- calculated_df %>% group_by(Target, Sample) %>% summarise(mean_RE = mean(RE), sd = sd(RE)) %>% ungroup()

    # TODO combine metadata
    if (is.null(metadata_df())) {
      calculated_df_summary <- calculated_df_summary
    } else{
      calculated_df_summary <- left_join(calculated_df_summary, metadata_df())
      calculated_df_summary <- calculated_df_summary %>% rename(Target = "Name", Target_primer = "Target")
    }

    # output final df in list, 1st element for render table, other one for ggplot
    list(df, calculated_df_summary)
  })

  # output combined df
  output$contents <- renderTable({
    #combined_df()
    calculated_df_list() %>% .[[2]]
  }, align = "c", width = "100%")


  # download data -----------------------------------------------------------
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("dataset-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(calculated_df_list() %>% .[[2]], file)
    })



  # plot --------------------------------------------------------------------
  output$plot <- renderPlotly({
    if (is.null(input$cq_files))
      return(NULL)
    if (is.null(input$file_input_efficiency))
      return(NULL)

    df_final <- calculated_df_list() %>% .[[2]]

    # perform transformation
    # if (input$logfunc_data == "log2_1") {
    #   df_final <- df_final %>% mutate(mean_RE = log2(mean_RE), sd = log2(sd))
    # }else{
    #   df_final <- calculated_df_list() %>% .[[2]]
    # }

    # plot
    fig <- df_final %>% filter(!Target %in% c(input$selected_ctrl_gene_1, input$selected_ctrl_gene_2)) %>%
      ggplot(aes(x = Sample, y = mean_RE, fill = Sample)) +
      facet_wrap(~ Target, ncol = input$facet_cols) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_errorbar(aes(ymin=mean_RE-sd, ymax=mean_RE+sd), width=input$whiskersize_data, position=position_dodge(.9)) +
      cowplot::theme_cowplot(input$textsize_data) + theme(axis.text.x = element_blank()) + ylab("Relative expression (A.U.)") +
      cowplot::panel_border(color = "black")


    if (input$free_y == T) {
      fig <- df_final %>% filter(!Target %in% c(input$selected_ctrl_gene_1, input$selected_ctrl_gene_2)) %>%
        ggplot(aes(x = Sample, y = mean_RE, fill = Sample)) +
        facet_wrap(~ Target, ncol = input$facet_cols, scale = "free_y") +
        geom_bar(stat = "identity", position = "dodge") +
        geom_errorbar(aes(ymin=mean_RE-sd, ymax=mean_RE+sd), width=input$whiskersize_data, position=position_dodge(.9)) +
        cowplot::theme_cowplot(input$textsize_data) + theme(axis.text.x = element_blank()) + ylab("Relative expression (A.U.)") +
        cowplot::panel_border(color = "black")
    }

    print(ggplotly(fig
                   # , tooltip = c("combined_md")
                   ))

  })


})
