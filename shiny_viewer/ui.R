library(shiny)

# Define UI for application that simulates gene expression and plots heatmap
shinyUI(pageWithSidebar(
    
    # Application title
    headerPanel("View the results of domain enrichment analysis"),
    
    # Sidebar with a slider input for number of observations
    sidebarPanel(
        textInput(inputId = "path", label = "path to *.RData file containing the results", value = "../processed_data_files/what_we_find_VS_ELM_clust11092017.RData"),
        selectInput(inputId = "enrich_plot_set", label = "Our prediction vs ELM domains: which results to compare?",
                    choices = "enrichment", selected = "enrichment", multiple = T),
        radioButtons(inputId = "enrich_plot_type", label = "Our prediction vs ELM domains: plot type",
                     choices = list(pval = "pval", odds_ratio = "odds_ratio", count = "count"), selected = "count"),
        textInput(inputId = "enrich_plot_args", label = "Our prediction vs ELM domains: plot args (like \"cex = 1.5\", separated by pipes)", value = "cex.lab = 2|cex.axis = 1.5"),
        textInput(inputId = "enrich_legend_args", label = "Our prediction vs ELM domains: legend args (like \"cex = 1.5\", separated by pipes)", value = "cex = 2"),
        numericInput(inputId = "leg_pos_x", label = "legend x position", value = 100, min = 0, max = 2000),
        checkboxInput(inputId = "show_known_domains", label = "Show known domains (when plotting count)?", value = FALSE),
        checkboxInput(inputId = "show_total_domains", label = "Show total domains found (when plotting count)?", value = FALSE),
        radioButtons(inputId = "bin2d_pval_plot_set", label = "Characteristics of top protein-domain pairs and pvalue distribution: which results to look at?",
                     choices = "unavailable"),
        textInput(inputId = "bin2d_plotname", label = "Characteristics of top protein-domain pairs: plot name", value = ""),
        selectInput(inputId = "bin2d_plot_rankby", label = "Characteristics of top protein-domain pairs: rank by column",
                    choices = "p.value", selected = "p.value", multiple = F),
        textInput(inputId = "bin2d_plot_filter", label = "Characteristics of top protein-domain pairs: filter criteria", value = "p.adjust(p.value, method = \"fdr\") < 0.05"),
        sliderInput(inputId = "N_pairs", label = "Characteristics of top protein-domain pairs: how many top pairs to choose",
                    min = 10, max = 5000, step = 10,
                    value = 250),
        textInput(inputId = "pval_plotname", label = "Empirical pvalue plot: plot name", value = "")
        ),
    
    # Show a plot of the generated distribution
    mainPanel(
        # plot heatmap
        tabsetPanel(
            tabPanel(title = "Our prediction vs ELM domains",
                     plotOutput("enrich_plot", height = "1000px") ###################### modify plot size here ###########
            ),
            tabPanel(title = "Characteristics of top protein-domain pairs",
                     plotOutput("bin2d_plot", height = "1200px") ###################### modify plot size here ###########
            ),
            tabPanel(title = "Empirical pvalue plot",
                     plotOutput("pval_plot", height = "600px") ###################### modify plot size here ###########
        )
    )
    
)
))

