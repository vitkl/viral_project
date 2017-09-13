library(shiny)
library(MItools)
library(rtracklayer)
library(ggplot2)
library(GGally)
library(RColorBrewer)
library(R.utils)

# Define server logic required to simulate gene expression and plot heatmap
shinyServer(function(input, output, session) {
    # identify missing packages
    packages = c("MItools", "ggplot2", "rtracklayer", "data.table", "GGally", "RColorBrewer")
    missing_packages = packages[!packages %in% names(installed.packages()[,"Package"])]
    if(length(missing_packages) > 0) stop(paste0("required packages are missing: ", paste(missing_packages, collapse = " ")))
    
    observe({
        if(grepl("gz",input$path)) {
            stop("gunzip *.RData first")
        } else path = input$path
        load(path)
        enrich_plot_set_choices = grep("enrichment",ls(), value = T)
        which_matrix = sapply(enrich_plot_set_choices, function(x) class(eval(parse(text = x)))) == "matrix"
        enrich_plot_set_choices = enrich_plot_set_choices[which_matrix]
        updateSelectInput(session, inputId = "enrich_plot_set", 
                          choices = enrich_plot_set_choices,
                          selected = enrich_plot_set_choices)
        
        which_EmpiricalPval = sapply(ls(), function(x) class(eval(parse(text = x)))) == "XYZinteration_XZEmpiricalPval"
        for_2d_bin = names(which_EmpiricalPval)[which_EmpiricalPval]
        updateRadioButtons(session, inputId = "bin2d_pval_plot_set", choices = for_2d_bin,
                           selected = for_2d_bin[1])
    })
    
    
    
    output$bin2d_plot <- renderPlot({
        PermutResult2D(res = eval(parse(text = input$bin2d_pval_plot_set)), N = input$N_pairs) +
            ggtitle(input$bin2d_plotname)
    })
    
    output$pval_plot <- renderPlot({
        plot(eval(parse(text = input$bin2d_pval_plot_set)), main = input$pval_plotname)
    })
    
    output$enrich_plot <- renderPlot({
        load(input$path)
        #which_function = sapply(ls(), function(x) class(eval(parse(text = x)))) == "function"
        #rm(list = names(which_function)[which_function])
        rm(plotEnrichment)
        
        results2plot = list()
        for (i in 1:(length(input$enrich_plot_set))) {
            results2plot = c(results2plot, list(eval(parse(text = input$enrich_plot_set[i]))))
        }
        
        if(is.null(input$enrich_plot_args)) plot_args = NULL else {
            plot_args = input$enrich_plot_args
            plot_args = unlist(strsplit(plot_args, "\\|"))
            }
        if(is.null(input$enrich_legend_args)) legend_args = NULL else {
            legend_args = input$enrich_legend_args
            legend_args = unlist(strsplit(legend_args, "\\|"))
            }
        par(mar = c(6,7,4,4))
        plotEnrichment(runningTestEnrichmentlist = results2plot,
                       random_domains = enrichmentRANDOM, 
                       domains_known_mapped = domains_known_mapped,
                       type = input$enrich_plot_type, plot_type = "l",
                       plot_args = plot_args, legend_args = legend_args
                       )
        # cex.lab = 2 cex = 2
    })
})
