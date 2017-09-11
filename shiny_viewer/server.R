library(shiny)

# Define server logic required to simulate gene expression and plot heatmap
shinyServer(function(input, output, session) {
    # identify missing packages
    packages = c("MItools", "ggplot2", "rtracklayer", "data.table", "GGally", "RColorBrewer")
    missing_packages = packages[!packages %in% names(installed.packages()[,"Package"])]
    if(length(missing_packages) > 0) stop(paste0("required packages are missing: ", paste(missing_packages, collapse = " ")))
    
    library(MItools)
    library(rtracklayer)
    library(ggplot2)
    library(GGally)
    library(RColorBrewer)
    library(R.utils)
    
    if(grepl("gz",input$path)) {
        path.gz = input$path
        gunzip(path, overwrite = T, remove = F)
        path = gsub("\\.gz$","", path.gz)
        } else path = input$path
    load(path)
    
})
    