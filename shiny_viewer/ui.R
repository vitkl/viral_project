library(shiny)

# Define UI for application that simulates gene expression and plots heatmap
shinyUI(pageWithSidebar(
    
    # Application title
    headerPanel("View the results of domain enrichment analysis"),
    
    # Sidebar with a slider input for number of observations
    sidebarPanel(
        textInput(inputId = "path", label = "path to *.RData file containing the results", value = "/Users/vitalii/Desktop/vitalii/viral_project/processed_data_files/what_we_find_VS_ELM_output_freq.RData"),
        selectInput(inputId = "enrich_plot_set", label = "Our prediction vs ELM domains: which results to compare?",
                    choices = "unavailable", multiple = T),
        radioButtons(inputId = "enrich_plot_type", label = "Our prediction vs ELM domains: plot type",
                     choices = list(pval = "pval", odds_ratio = "odds_ratio", count = "count"), selected = "count"),
        radioButtons(inputId = "bin2d_plot_set", label = "Characteristics of top protein-domain pairs: which results to look at?",
                     choices = "unavailable"),
        sliderInput(inputId = "N_pairs", label = "Characteristics of top protein-domain pairs: how many top pairs to choose",
                    min = 10, max = 5000, step = 10,
                    value = 250),
        radioButtons(inputId = "pval_plot_set", label = "Empirical pvalue plot: which results to look at?",
                     choices = "unavailable")
        ),
    
    # Show a plot of the generated distribution
    mainPanel(
        # plot heatmap
        tabsetPanel(
            tabPanel(title = "Our prediction vs ELM domains",
                     plotOutput("enrich_plot", height = "auto")
            ),
            tabPanel(title = "Characteristics of top protein-domain pairs:",
                     plotOutput("bin2d_plot", height = "auto")
            ),
            tabPanel(title = "Empirical pvalue plot",
                     plotOutput("pval_plot", height = "auto")
        )
    )
    
)
))

