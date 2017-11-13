library(shiny)
library(MItools)
library(rtracklayer)
library(ggplot2)
library(GGally)
library(RColorBrewer)
library(R.utils)
library(ROCR)

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
                          selected = c("enrichment", "enrichment_justfreq", "enrichmentFISHER_justodds", "enrichmentFISHER_justpval"))
        
        which_EmpiricalPval = sapply(ls(), function(x) class(eval(parse(text = x)))) == "XYZinteration_XZEmpiricalPval"
        for_2d_bin = names(which_EmpiricalPval)[which_EmpiricalPval]
        if(input$bin2d_pval_plot_set == "unavailable"){
            updateRadioButtons(session, inputId = "bin2d_pval_plot_set", choices = for_2d_bin,
                               selected = for_2d_bin[1])
        }
        if(input$ROCR_plot_set[1] == "unavailable"){
            updateSelectInput(session, inputId = "ROCR_plot_set", 
                              choices = for_2d_bin,
                              selected = c("res_count", "resJustFISHER"))
        }

        to_choose_rankby = isolate({input$bin2d_pval_plot_set %in% ls()})
        if(to_choose_rankby){
            choicesLISTrankby = colnames(eval(parse(text = input$bin2d_pval_plot_set))$data_with_pval)
            choicesLISTrankby = choicesLISTrankby[sapply(eval(parse(text = input$bin2d_pval_plot_set))$data_with_pval, class) == "numeric"]
            if(!"choicesLISTrankby" %in% ls()) choicesLISTrankby = "p.value"
            updateSelectInput(session, inputId = "bin2d_plot_rankby", 
                              choices = choicesLISTrankby,
                              selected = c("p.value"))
        }
    })
    
    
    
    output$bin2d_plot <- renderPlot({
        load(input$path)
        PermutResult2D(res = eval(parse(text = input$bin2d_pval_plot_set)),
                       N = input$N_pairs, rank_by = input$bin2d_plot_rankby,
                       filter = input$bin2d_plot_filter) +
            ggtitle(input$bin2d_plotname)
    })
    
    output$pval_plot <- renderPlot({
        load(input$path)
        plot(eval(parse(text = input$bin2d_pval_plot_set)), main = input$pval_plotname)
    })
    
    output$enrich_plot <- renderPlot({
        load(input$path)
        #which_function = sapply(ls(), function(x) class(eval(parse(text = x)))) == "function"
        #rm(list = names(which_function)[which_function])
        #rm(plotEnrichment)
        
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
        print(as.logical(input$show_total_domains))
        plotEnrichment(runningTestEnrichmentlist = results2plot,
                       random_domains = enrichmentRANDOM, 
                       domains_known_mapped = domains_known_mapped,
                       type = input$enrich_plot_type, plot_type = "l",
                       plot_args = plot_args, legend_args = legend_args,
                       leg_pos_x = input$leg_pos_x,
                       show_known_domains = as.logical(input$show_known_domains),
                       plot_total_domains_found = as.logical(input$show_total_domains)
                       )
        # cex.lab = 2 cex = 2
    })
    output$ROCR <- renderPlot({
        load(input$path)
        res_list = eval(parse(text = paste0("list(",paste0(input$ROCR_plot_set, collapse = ","),")")))
        predictor = lapply(res_list, function(XYZint){
            XYZint_table = unique(XYZint$data_with_pval[,c(XYZint$nodes$nodeX, XYZint$nodes$nodeZ, "p.value"), with = F])
            XYZint_table[, p.value := 1 - p.value]
            
            #XYZint_table = unique(XYZint$data_with_pval[,c(XYZint$nodes$nodeZ, "p.value"), with = F])
            #XYZint_table[, p.value := 1 - min(p.value), by = eval(XYZint$nodes$nodeZ)]
            #XYZint_table = unique(XYZint_table)
            
            #XYZint_table = XYZint_table[order(p.value, decreasing = F)[1:input$N_pairs]]
            XYZint_table$p.value
        })
        names(predictor) = input$ROCR_plot_set
        truth = lapply(res_list, function(XYZint){
            XYZint_table = unique(XYZint$data_with_pval[,c(XYZint$nodes$nodeX, XYZint$nodes$nodeZ, "p.value"), with = F])
            XYZint_table[, p.value := 1 - p.value]
            
            #XYZint_table = unique(XYZint$data_with_pval[,c(XYZint$nodes$nodeZ, "p.value"), with = F])
            #XYZint_table[, p.value := 1 - min(p.value), by = eval(XYZint$nodes$nodeZ)]
            #XYZint_table = unique(XYZint_table)
            
            #XYZint_table = XYZint_table[order(p.value, decreasing = F)[1:input$N_pairs]]
            eval(parse(text = paste0("XYZint_table[",XYZint$nodes$nodeZ," %in% domains_known_mapped, truth := 1]")))
            eval(parse(text = paste0("XYZint_table[!",XYZint$nodes$nodeZ," %in% domains_known_mapped, truth := 0]")))
            XYZint_table$truth
        })
        names(truth) = input$ROCR_plot_set

        pred <- ROCR::prediction(predictor, truth)
        perf <- ROCR::performance(pred, "tpr", "fpr")
        auc.perf = ROCR::performance(pred, measure = "auc")
        
        if(is.null(input$enrich_plot_args)) plot_args = NULL else {
            plot_args = input$enrich_plot_args
            plot_args = unlist(strsplit(plot_args, "\\|"))
        }
        par(mar = c(6,7,4,4))
        color = colorRampPalette(brewer.pal(7, "Dark2"))(length(input$ROCR_plot_set))
        eval(parse(text = paste0("plot(perf, colorize=F, lwd = 10,",paste0(plot_args, collapse = ","),")")))
        abline(a=0, b= 1)
        if(is.null(input$enrich_legend_args)) legend_args = NULL else {
            legend_args = input$enrich_legend_args
            legend_args = unlist(strsplit(legend_args, "\\|"))
        }
        eval(parse(text = paste0("text(x = 0.5, y = auc.perf@y.values, col = \"black\", labels = paste0(names(pred@predictions),\", AUC: \",signif(as.numeric(auc.perf@y.values), 3)),",paste0(legend_args, collapse = ","),")")))
        })
})


