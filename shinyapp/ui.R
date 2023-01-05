#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)
library(shinyFiles)



# fileinput
# input for qpcr files containing cq data
file_input_expression <- fileInput("cq_files",
                                   "CSV files containing cq values",
                                   multiple = TRUE,
                                   accept=c('text/csv',
                                            'text/comma-separated-values,text/plain',
                                            '.csv'))

# metadata of files containing primer and ID
file_input_metadata <- fileInput("file_input_md", "Choose csv file of metadata: 1st col Target, 2nd col Name (optional)",
                                 accept = c(
                                   "text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")
)

# file containing efficiency
file_input_efficiency <- fileInput("file_input_efficiency", "Choose csv file containing efficiency (1st col target, 2nd col efficiency",
                                   accept = c(
                                     "text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv")
)



# logfunc for data transformation
logfunc <- radioButtons("logfunc_data", h4("Data transformation method"),
                        choices = list("log2(n+1)" = "log2_1",
                                       "none" = "none"),
                        selected = "none")


# number of column
facetcols <- numericInput("facet_cols", "Number of columns for facet", value = 3, min = 1, max = NA, step = 0.5)


# text size
textsize <- sliderInput("textsize_data", "Text size",
                        min = 1, max = 25, value = c(15))

whiskersize <- sliderInput("whiskersize_data", "Whisker width",
                           min = 0.01, max = 0.5, value = c("0.05"))

# metadata to show on hover
# metadata_hover <- sliderInput("metadata_hover_data", "Metadata to show upon hover",
#                               min = 1, max = 25, value = c(12))

# selection for control genes
control_gene_1 <- uiOutput("gene_select_1")
control_gene_2 <- uiOutput("gene_select_2")

# selection for control samples
ui_control_sample <- uiOutput("control_sample")
ui_control_sample_replicate <- uiOutput("control_sample_replicate")

# selection for data to show (checklist)
ui_genes_select <- uiOutput("genes_select")

# download button
dl_button <- downloadButton("downloadData", "Download normalized data")

# Define UI
shinyUI(pageWithSidebar(

  # Application title
  headerPanel("Visualization of gene expression levels"),

  sidebarPanel(
    file_input_expression,
    file_input_metadata,
    file_input_efficiency,
    control_gene_1,
    control_gene_2,
    ui_control_sample,
    ui_control_sample_replicate,
    ui_genes_select,
    # logfunc,
    whiskersize,
    facetcols,
    textsize,
    dl_button
  ),

  mainPanel(plotlyOutput("plot", height = "600px"),
            textOutput("tbl_dims"),
            div(style="height:600px; overflow:scroll; align:text-center",
                tableOutput("contents")
            )
  )
))
