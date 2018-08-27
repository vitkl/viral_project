---
title: "Comprehensive benchmarking of motif instances (all types) - only IntAct and BioPlex"
author: "Vitalii Kleshchevnikov"
date: "18/10/2017"
output: 
  html_document: 
    keep_md: yes
    toc: yes
---


```r
knitr::opts_chunk$set(fig.width = 12, fig.height = 8, warning = FALSE, message = FALSE)

packages = c("MItools", "RColorBrewer", "devtools", "grid","gridExtra", "ggplot2", "VennDiagram", "svglite", "GenomicRanges", "Biostrings", "RColorBrewer")
if(mean(packages %in% names(installed.packages()[,"Package"])) != 1){
    packages_to_install = packages[!packages %in% names(installed.packages()[,"Package"])]
    # specifying mirror is necessary for some Linux systems
    install.packages(packages_to_install, dependencies = T, repos = "http://mirrors.ebi.ac.uk/CRAN/")
    packages_to_install = packages[!packages %in% names(installed.packages()[,"Package"])]
    source("https://bioconductor.org/biocLite.R")
    biocLite(packages_to_install)
    devtools::install_github("vitkl/MItools", dependencies = T)
}
suppressPackageStartupMessages({
    library(MItools)
    library(grid)
    library(gridExtra)
    library(ggplot2)
    library(GenomicRanges)
    library(Biostrings)
    library(clustermq)
    library(RColorBrewer)
})
```

```
## Warning: replacing previous import 'IRanges::desc' by 'plyr::desc' when
## loading 'PSICQUIC'
```

```
## Warning: replacing previous import 'GenomicRanges::shift' by
## 'data.table::shift' when loading 'MItools'
```

```r
colors = RColorBrewer::brewer.pal(2, "Dark2")
```

```
## Warning in RColorBrewer::brewer.pal(2, "Dark2"): minimal value for n is 3, returning requested palette with 3 different levels
```

```r
motif_types = c("MOD", "LIG", "DOC") # "MOD", "LIG", "DOC" / "DEG", "CLV", "TRG"
MOD_LIG_DOC = T # if TRUE analyse only "MOD", "LIG", "DOC" motifs
normalise = F # if normalised Sig threshold of 0.3 becomes 1 (but this equalises Sig == 0.3 to true 1 - that is when motif was not found)
nice_colors = brewer.pal(6, "Dark2")

#ROC:
measure1 = "prec" # "tpr" or "prec"
measure2 = "rec" # "fpr" or "rec"
both_metrics_vs_cutoff = T
single_metric = c("auc","prbe")[1]
single_metric_name = c("Median AUC", "Prec-rec break-even")[1]

datasets = c("qslimfinder.Full_IntAct3.FALSE",
             "qslimfinder.BioPlex3cloudfixF.FALSE",
             "qslimfinder.all_viral_interaction3.FALSE")
descriptions = c("including all human data (IntAct)",
                 "including human data (BioPlex)",
                 "only viral-human data")
motif_setup_month = c("201802",
                      "201804",
                      "201802")

Sys.Date. = Sys.Date()
Sys.Date. = as.Date("2018-08-23")
Sys.Date.
```

```
## [1] "2018-08-23"
```

```r
filename = paste0("./processed_data_files/compr_benchmarking_venn_IntVidVir_human_dom_",format(Sys.Date., "%Y%m"),"_type_", measure1, measure2,".RData")
filename
```

```
## [1] "./processed_data_files/compr_benchmarking_venn_IntVidVir_human_dom_201808_type_precrec.RData"
```

```r
if(file.exists(filename)) load(filename)
```

## Overview

### 
1. filter_by_domain_data: NULL or "p.value < 0.5"
2. neg_set: "all_instances", "all_proteins"
3. motif_pval_cutoff: 1, precision == recall, precision > 0.5


```r
#####################################################################################
myBenchmarkMotifs = function(neg_set, filter_by_domain_data, motif_pval_cutoff,
                             dataset, datasets,
                             non_query_domain_results_obj = NULL, # res_count_all
                             query_domains_only = F, min_non_query_domain_support = 0, min_top_domain_support4motif_nq = 0,
                             return_all = F, merge_domain_data = F, motif_setup_obj2 = NULL, occurence_filt = NULL){
    
    if(isTRUE(grepl("all_viral_interaction",dataset))){
        non_query_domain_res_file = "./processed_data_files/domain_res_count_20171019.RData"
    } else non_query_domain_res_file = "./processed_data_files/predict_domain_human_clust20180819.RData"
    
    res = mBenchmarkMotifs(datasets = dataset,
                           descriptions = descriptions[datasets == dataset],
                           dir = "./",
                           motif_setup_months = motif_setup_month[datasets == dataset],
                           
                           domain_res_file = "./processed_data_files/domain_res_count_20171019.RData",
                           neg_set = neg_set,
                           domain_results_obj = "res_count", motif_input_obj = "forSLIMFinder_Ready",
                           motif_setup_obj2 = motif_setup_obj2, occurence_filt = occurence_filt,
                           one_from_cloud = T,
                           dbfile_main = "./data_files/instances9606.gff",
                           dburl_main = "http://elm.eu.org/instances.gff?q=None&taxon=Homo%20sapiens&instance_logic=",
                           dbfile_query = "./data_files/instances10239.gff",
                           dburl_query = "http://elm.eu.org/instances.gff?q=all&taxon=irus&instance_logic=",
                           query_res_query_only = T, motif_types = motif_types,
                           all_res_excl_query = T, merge_motif_variants = T,
                           seed = 21, N = 100, replace = T, within1sequence = T,
                           query_predictor_col = "Sig", all_predictor_col = "Sig", normalise = normalise,
                           minoverlap = 2, maxgap = 0,
                           minoverlap_redundant = 5,
                           merge_domain_data = merge_domain_data,
                           merge_by_occurence_mcols = c("query", "interacts_with"),
                           merge_by_domain_res_cols = c("IDs_interactor_viral", "IDs_interactor_human", "IDs_domain_human", "Taxid_interactor_human","Taxid_interactor_viral"),
                           merge_by_non_query_domain_res_cols = c("IDs_interactor_human_A", "IDs_interactor_human_B", "IDs_domain_human_B", "Taxid_interactor_human_A","Taxid_interactor_human_B"),
                           filter_by_domain_data = filter_by_domain_data, motif_pval_cutoff = motif_pval_cutoff,
                           select_predictor_per_range = max,
                           non_query_domain_res_file = non_query_domain_res_file,
                           non_query_domain_results_obj = non_query_domain_results_obj, 
                           non_query_domains_N = 0,
                           non_query_set_only = F,
                           query_domains_only = query_domains_only,
                           min_non_query_domain_support = min_non_query_domain_support,
                           min_top_domain_support4motif_nq = min_top_domain_support4motif_nq,
                           select_top_domain = F)
    
    DT = data.table(N_query_prot_with_known_instances = res[[1]]$N_query_prot_with_known_instances,
                    N_query_known_instances = res[[1]]$N_query_known_instances,
                    N_query_prot_with_known_instances_found = res[[1]]$N_query_prot_with_known_instances_found,
                    N_query_found_match_known_instances = res[[1]]$N_query_known_instances_found,
                    N_query_total_instances_found = res[[1]]$N_query_total_instances_found,
                    N_query_known_match_instances_found = res[[1]]$N_query_match_known_instances_found)
    if(return_all) list(res = res, DT = DT) else return(DT)
}
#####################################################################################
#library(lineprof)
#l <- lineprof(myBenchmarkMotifs(neg_set = "all_instances",
#                                   filter_by_domain_data = "p.value < 0.5",
#                                   motif_pval_cutoff = 1,
#                                   dataset = datasets[1], datasets,
#                                   non_query_domain_results_obj = "res_count_all", # res_count_all
#                                   query_domains_only = F, min_non_query_domain_support = 4,
#                                   return_all = F, merge_domain_data = T), interval = 0.1)
#shine(l) 
#sapply(colnames(benchmark), function(dataset) {
#    prec05cutoff = 1 - min(benchmark["ROC_prediction",dataset][[1]]@cutoffs[[1]][
#       benchmark["ROC_performance",dataset][[1]]@y.values[[1]] > 0.5], na.rm = T)
#    known_at_prec05cutoff = res_list_viral$qslimfinder.Full_IntAct3cloudfixF.FALSE$overlapping_GRanges_query[
#        res_list_viral$qslimfinder.Full_IntAct3cloudfixF.FALSE$overlapping_GRanges_query$Sig <= prec05cutoff]
#    unique(paste0(seqnames(known_at_prec05cutoff), "=", known_at_prec05cutoff$ID))
#})
#####################################################################################
plotVenn = function(res, ...){
    #grid.newpage()
    if(res$N_query_total_instances_found < res$N_query_known_instances) rotation.degree = 180 else rotation.degree = 0
    venn.plot <- VennDiagram::draw.pairwise.venn(
        area1 = res$N_query_total_instances_found,
        area2 = res$N_query_known_instances,
        cross.area = res$N_query_known_match_instances_found,
        category = c("discovered instances", "known instances"),
        fill = c("royalblue1", "seagreen1"),
        euler.d = F, scaled = F,
        lty = "blank",
        cat.pos = c(340, 175),
        cat.dist = c(0.05,0.05),
        ext.pos = 10,
        ext.dist = -0.05,
        ext.length = 0.8,
        ext.line.lwd = 2,
        ext.line.lty = "dashed",
        rotation.degree = rotation.degree,
        ind = F,
        ...
    );
    #grid.newpage()
    venn = gTree(children=venn.plot, name="venn")
    #grid.newpage()
    # since multiple matches are possible between datasets I need to edit numbers shown on a diagram
    venn = editGrob(venn, paste0(venn[["childrenOrder"]][7]),
                    label = paste0(res$N_query_found_match_known_instances,
                                   " / ",res$N_query_known_match_instances_found),
                    grep = T)
    #grid.newpage()
    if(rotation.degree == 0) discovered_ind = 5 else discovered_ind = 6
    venn = editGrob(venn, paste0(venn[["childrenOrder"]][discovered_ind]),
                    label = res$N_query_total_instances_found - res$N_query_found_match_known_instances,
                    grep = T)
    #grid.newpage()
    #grid.draw(venn)
    return(venn)
}
#####################################################################################
#plotVenn(res_list_viral)
#####################################################################################
plotAnnotVenn = function(res, main = "", sub = "",
                         main_fontsize = 22, sub_fontsize = 17,
                         main_lineheight = 1.4, sub_lineheight = 0.9,
                         cex = 1.8, cat.cex = 1.8,
                         fontfamily = "Arial", ...) { # "serif" "Arial"
    
    lay = matrix(c(1), nrow = 1, ncol = 1)
    annotated_venn = plotVenn(res, fontfamily = fontfamily, cat.fontfamily = fontfamily,
                              cex = cex, cat.cex = cat.cex)
    if(sub != "") {
        lay = matrix(c(1,2,2,2,2,2), nrow = 6, ncol = 1)
        sub_grob = textGrob(sub,
                            gp=gpar(fontsize=sub_fontsize, fontfamily = fontfamily,
                                    lineheight = sub_lineheight, ...))
        annotated_venn = gridExtra::arrangeGrob(sub_grob,
                                                annotated_venn, ncol=1,
                                                layout_matrix = lay)
    }
    
    if(main != ""){
        lay = matrix(c(1,2,2,2,2,2), nrow = 6, ncol = 1)
        main_grob = textGrob(main,
                             gp=gpar(fontsize=main_fontsize, fontfamily = fontfamily,
                                     lineheight = main_lineheight, ...))
        annotated_venn = gridExtra::arrangeGrob(main_grob,
                                                annotated_venn, ncol=1,
                                                layout_matrix = lay)
    } 
    
    #grid.newpage()
    #grid.draw(annotated_venn)
    return(annotated_venn)
}
#####################################################################################
#plotAnnotVenn(res_list_viral, main = "corrected p-value < 0.3")
#####################################################################################
mPlotAnnotVenn = function(filter_by_domain_data = NULL,
                          motif_pval_cutoffs = c(0.3, 0.05, 0.002),
                          motif_pval_tags = c("lenient cutoff", "optimal cutoff", "stringent cutoff"),
                          dataset = datasets[1], datasets = datasets, fontfamily = "Arial",
                          non_query_domain_results_obj = NULL, # res_count_all
                          query_domains_only = F, min_non_query_domain_support = 0, min_top_domain_support4motif_nq = 0,
                          return_all = F, merge_domain_data = F, motif_setup_obj2 = NULL, occurence_filt = NULL) {
    #function_4export = ls(envir = .GlobalEnv)
    #function_4export = function_4export[sapply(function_4export, function(func) class(eval(parse(text = func)))) == "function"]
    #export = lapply(function_4export, function(func) eval(parse(text = func), enclos = .GlobalEnv)) 
    #names(export) = function_4export
    #message(export)
    
    #annotated_venns = Q(fun = function(motif_pval_cutoff){
    annotated_venns = lapply(motif_pval_cutoffs, function(motif_pval_cutoff){
        res_instances = myBenchmarkMotifs(neg_set = "all_instances",
                                          filter_by_domain_data = filter_by_domain_data,
                                          motif_pval_cutoff = motif_pval_cutoff,
                                          dataset = dataset, datasets,
                                          non_query_domain_results_obj = non_query_domain_results_obj, # res_count_all
                                          query_domains_only = query_domains_only, min_non_query_domain_support = min_non_query_domain_support,
                                          min_top_domain_support4motif_nq = min_top_domain_support4motif_nq,
                                          return_all = return_all, merge_domain_data = merge_domain_data,
                                          motif_setup_obj2 = motif_setup_obj2, occurence_filt = occurence_filt)
        res_proteins = myBenchmarkMotifs(neg_set = "all_proteins",
                                         filter_by_domain_data = filter_by_domain_data,
                                         motif_pval_cutoff = motif_pval_cutoff,
                                         dataset = dataset, datasets,
                                         non_query_domain_results_obj = non_query_domain_results_obj, # res_count_all
                                         query_domains_only = query_domains_only, min_non_query_domain_support = min_non_query_domain_support,
                                         min_top_domain_support4motif_nq = min_top_domain_support4motif_nq,
                                         return_all = return_all, merge_domain_data = merge_domain_data,
                                         motif_setup_obj2 = motif_setup_obj2, occurence_filt = occurence_filt)
        annotated_venn = plotAnnotVenn(res_instances,
                                       main = paste0(motif_pval_tags[motif_pval_cutoffs == motif_pval_cutoff],
                                                     ": adj. p-value < ", signif(motif_pval_cutoff, 2)),
                                       sub = paste0(res_instances$N_query_prot_with_known_instances,
                                                    " proteins with known SLIMs (",
                                                    res_instances$N_query_prot_with_known_instances_found," proteins recovered)\n",
                                                    res_proteins$N_query_total_instances_found," SLIMs predicted in other proteins"),
                                       fontfamily = fontfamily)
    }) #, motif_pval_cutoffs,
    #export = c(list(filter_by_domain_data = filter_by_domain_data,
    #                motif_pval_tags = motif_pval_tags, dataset = dataset,
    #                datasets = datasets, fontfamily = fontfamily), export),
    ##export = list(myBenchmarkMotifs = myBenchmarkMotifs, plotAnnotVenn = plotAnnotVenn),
    #seed = 128965,
    #memory = 8000, n_jobs = 3)
    gridExtra::arrangeGrob(grobs = annotated_venns, ncol=1)
}
#####################################################################################
generateCutoffs = function(rocr_list, lenient_cutoff, stringent_prec){
    cutoffs = rocr_list[1,1][[1]]@cutoffs[[1]]
    prec = rocr_list[2,1][[1]]@y.values[[1]]
    rec = rocr_list[2,1][[1]]@x.values[[1]]
    motif_pval_cutoffs = c(lenient_cutoff,
                                   1 - min(cutoffs[prec >= rec & cutoffs != 0], na.rm = T),
                                   1 - min(cutoffs[prec >= stringent_prec & cutoffs != 0], na.rm = T))
    motif_pval_cutoffs[motif_pval_cutoffs == -Inf] = min(motif_pval_cutoffs[motif_pval_cutoffs != -Inf])
    motif_pval_cutoffs
}
#####################################################################################
figureVenn = function(filter_by_domain_data = NULL,
                      main = "",
                      motif_pval_tags = c("lenient cutoff", "optimal cutoff", "stringent cutoff"),
                      datasets = datasets,
                      stringent_prec = 0.5,
                      lenient_cutoff = 0.3,
                      motif_pval_cutoffs_IntAct = NULL,
                      motif_pval_cutoffs_BioPlex = NULL,
                      motif_pval_cutoffs_viral = NULL,
                      fontfamily = "Arial",
                      non_query_domain_results_obj = NULL, # res_count_all
                      query_domains_only = F,
                      min_non_query_domain_support = 0, min_top_domain_support4motif_nq = 0,
                      return_all = F, merge_domain_data = F, 
                      descriptions, viral_only_venn = T, ...){
    message(main)
    ################################################################################ IntAct
    res_list_IntAct = myBenchmarkMotifs(neg_set = "all_instances",
                                        filter_by_domain_data = filter_by_domain_data,
                                        motif_pval_cutoff = 1,
                                        dataset = datasets[1], datasets, return_all = T,
                                        non_query_domain_results_obj = non_query_domain_results_obj, # res_count_all
                                        query_domains_only = query_domains_only, min_non_query_domain_support = min_non_query_domain_support,
                                        min_top_domain_support4motif_nq = min_top_domain_support4motif_nq,
                                        merge_domain_data = merge_domain_data,
                                        motif_setup_obj2 = NULL, occurence_filt = NULL)
    rocr_list_IntAct = mBenchmarkMotifsROC(res_list_IntAct$res, measure1 = "prec", measure2 = "rec", ROCR_data = T)
    if(is.null(motif_pval_cutoffs_IntAct)){
        motif_pval_cutoffs_IntAct = generateCutoffs(rocr_list_IntAct, lenient_cutoff, stringent_prec)
        message(paste0("motif_pval_cutoffs_IntAct: ", paste0(motif_pval_cutoffs_IntAct, collapse = ", ")))
    }
    big_plot_IntAct = mPlotAnnotVenn(filter_by_domain_data = filter_by_domain_data,
                                     motif_pval_cutoffs = motif_pval_cutoffs_IntAct,
                                     motif_pval_tags = motif_pval_tags,
                                     dataset = datasets[1], datasets = datasets, fontfamily = fontfamily,
                                     non_query_domain_results_obj = non_query_domain_results_obj, # res_count_all
                                     query_domains_only = query_domains_only, min_non_query_domain_support = min_non_query_domain_support,
                                     min_top_domain_support4motif_nq = min_top_domain_support4motif_nq,
                                     return_all = return_all, merge_domain_data = merge_domain_data,
                                     motif_setup_obj2 = res_list_IntAct$res[[datasets[1]]]$motif_setup_obj2,
                                     occurence_filt = res_list_IntAct$res[[datasets[1]]]$occurence_filt)
    IntAct_grob = textGrob(descriptions[1],
                           gp=gpar(fontsize=25, fontfamily = fontfamily,
                                   lineheight = 1.4, ...))
    ################################################################################ BioPlex
    res_list_BioPlex = myBenchmarkMotifs(neg_set = "all_instances",
                                         filter_by_domain_data = filter_by_domain_data,
                                         motif_pval_cutoff = 1,
                                         dataset = datasets[2], datasets, return_all = T,
                                         non_query_domain_results_obj = non_query_domain_results_obj, # res_count_all
                                         query_domains_only = query_domains_only, min_non_query_domain_support = min_non_query_domain_support,
                                         min_top_domain_support4motif_nq = min_top_domain_support4motif_nq,
                                         merge_domain_data = merge_domain_data,
                                        motif_setup_obj2 = NULL, occurence_filt = NULL)
    rocr_list_BioPlex = mBenchmarkMotifsROC(res_list_BioPlex$res, measure1 = "prec", measure2 = "rec", ROCR_data = T)
    if(is.null(motif_pval_cutoffs_BioPlex)){
        motif_pval_cutoffs_BioPlex = generateCutoffs(rocr_list_BioPlex, lenient_cutoff, stringent_prec)
        message(paste0("motif_pval_cutoffs_BioPlex: ", paste0(motif_pval_cutoffs_BioPlex, collapse = ", ")))
    }
    big_plot_BioPlex = mPlotAnnotVenn(filter_by_domain_data = filter_by_domain_data,
                                      motif_pval_cutoffs = motif_pval_cutoffs_BioPlex,
                                      motif_pval_tags = motif_pval_tags,
                                      dataset = datasets[2], datasets = datasets, fontfamily = fontfamily,
                                      non_query_domain_results_obj = non_query_domain_results_obj, # res_count_all
                                      query_domains_only = query_domains_only, min_non_query_domain_support = min_non_query_domain_support,
                                      min_top_domain_support4motif_nq = min_top_domain_support4motif_nq,
                                      return_all = return_all, merge_domain_data = merge_domain_data,
                                     motif_setup_obj2 = res_list_BioPlex$res[[datasets[2]]]$motif_setup_obj2,
                                     occurence_filt = res_list_BioPlex$res[[datasets[2]]]$occurence_filt)
    BioPlex_grob = textGrob(descriptions[2],
                            gp=gpar(fontsize=25, fontfamily = fontfamily,
                                    lineheight = 1.4, ...))
    ################################################################################ viral only
    if(isTRUE(viral_only_venn)){
        res_list_viral = myBenchmarkMotifs(neg_set = "all_instances",
                                           filter_by_domain_data = filter_by_domain_data,
                                           motif_pval_cutoff = 1,
                                           dataset = datasets[3], datasets, return_all = T,
                                           non_query_domain_results_obj = non_query_domain_results_obj, # res_count_all
                                           query_domains_only = query_domains_only, min_non_query_domain_support = min_non_query_domain_support,
                                           min_top_domain_support4motif_nq = min_top_domain_support4motif_nq,
                                           merge_domain_data = merge_domain_data,
                                           motif_setup_obj2 = NULL, occurence_filt = NULL)
        rocr_list_viral = mBenchmarkMotifsROC(res_list_viral$res, measure1 = "prec", measure2 = "rec", ROCR_data = T)
        if(is.null(motif_pval_cutoffs_viral)){
            motif_pval_cutoffs_viral = generateCutoffs(rocr_list_viral, lenient_cutoff, stringent_prec)
            message(paste0("motif_pval_cutoffs_viral: ", paste0(motif_pval_cutoffs_viral, collapse = ", ")))
        } 
        
        big_plot_viral = mPlotAnnotVenn(filter_by_domain_data = filter_by_domain_data,
                                        motif_pval_cutoffs = motif_pval_cutoffs_viral,
                                        motif_pval_tags = motif_pval_tags,
                                        dataset = datasets[3], datasets = datasets, fontfamily = fontfamily,
                                        non_query_domain_results_obj = non_query_domain_results_obj, # res_count_all
                                        query_domains_only = query_domains_only, min_non_query_domain_support = min_non_query_domain_support,
                                        min_top_domain_support4motif_nq = min_top_domain_support4motif_nq,
                                        return_all = return_all, merge_domain_data = merge_domain_data,
                                        motif_setup_obj2 = res_list_viral$res[[datasets[3]]]$motif_setup_obj2,
                                        occurence_filt = res_list_viral$res[[datasets[3]]]$occurence_filt)
        viral_grob = textGrob(descriptions[3],
                              gp=gpar(fontsize=25, fontfamily = fontfamily,
                                      lineheight = 1.4, ...))
        big_plot_grobs = list(big_plot_IntAct, big_plot_BioPlex, big_plot_viral)
        big_text_grobs = list(IntAct_grob, BioPlex_grob, viral_grob)
    } else {
        res_list_viral = NULL
        rocr_list_viral = NULL
        big_plot_viral = NULL
        viral_grob = NULL
        motif_pval_cutoffs_viral = NULL
        big_plot_grobs = list(big_plot_IntAct, big_plot_BioPlex)
        big_text_grobs = list(IntAct_grob, BioPlex_grob)
    }
    big_plot = gridExtra::arrangeGrob(grobs = big_plot_grobs, nrow=1)
    big_text = gridExtra::arrangeGrob(grobs = big_text_grobs, nrow=1)
    big_plot_text = gridExtra::arrangeGrob(big_text, big_plot, ncol=1,
                                           layout_matrix = matrix(c(1,rep(2,6)),
                                                                  nrow = 7, ncol = 1))
    if(main != ""){
        main_grob = textGrob(main,
                             gp=gpar(fontsize=27, fontfamily = fontfamily,
                                     lineheight = 1, ...))
        big_plot_text = gridExtra::arrangeGrob(main_grob, big_plot_text, ncol=1,
                                               layout_matrix = matrix(c(1,rep(2,7)),
                                                                      nrow = 8, ncol = 1))
    }
    
    list(big_plot_text = big_plot_text,
         res_list_IntAct = res_list_IntAct, rocr_list_IntAct = rocr_list_IntAct, 
         res_list_BioPlex = res_list_BioPlex, rocr_list_BioPlex = rocr_list_BioPlex,
         res_list_viral = res_list_viral, rocr_list_viral = rocr_list_viral, 
         big_plot_IntAct = big_plot_IntAct, IntAct_grob = IntAct_grob,
         big_plot_BioPlex = big_plot_BioPlex, BioPlex_grob = BioPlex_grob,
         big_plot_viral = big_plot_viral, viral_grob = viral_grob,
         motif_pval_cutoffs_IntAct = motif_pval_cutoffs_IntAct,
         motif_pval_cutoffs_BioPlex = motif_pval_cutoffs_BioPlex,
         motif_pval_cutoffs_viral = motif_pval_cutoffs_viral)
}
#####################################################################################
replotFigureVenn = function(figure, main, fontfamily = "Arial", viral_only_venn = F, ...){
    IntAct_grob = textGrob("including all human data (IntAct)",
                           gp=gpar(fontsize=24, fontfamily = fontfamily,
                                   lineheight = 1, ...))
    BioPlex_grob = textGrob("including human data (BioPlex)",
                            gp=gpar(fontsize=24, fontfamily = fontfamily,
                                    lineheight = 1, ...))
    viral_grob = textGrob("only viral-human data",
                          gp=gpar(fontsize=24, fontfamily = fontfamily,
                                  lineheight = 1, ...))
    
    big_plot = gridExtra::arrangeGrob(rectGrob(gp=gpar(col = "white")),
                                      figure$big_plot_viral,
                                      rectGrob(gp=gpar(col = "white")),
                                      figure$big_plot_IntAct,
                                      rectGrob(gp=gpar(col = "white")),
                                      figure$big_plot_BioPlex, nrow=1,
                                      layout_matrix = matrix(c(1,rep(2,6),
                                                               3,rep(4,6),
                                                               5,rep(6,6)),
                                                             nrow = 1, ncol = 21))
    big_text = gridExtra::arrangeGrob(viral_grob, IntAct_grob, BioPlex_grob, nrow=1)
    big_plot_text = gridExtra::arrangeGrob(big_text, big_plot, ncol=1,
                                           layout_matrix = matrix(c(1,rep(2,7)),
                                                                  nrow = 8, ncol = 1))
    if(main != ""){
        main_grob = textGrob(main,
                             gp=gpar(fontsize=27, fontfamily = fontfamily,
                                     lineheight = 1, ...))
        big_plot_text = gridExtra::arrangeGrob(main_grob, big_plot_text, ncol=1,
                                               layout_matrix = matrix(c(1,rep(2,7)),
                                                                      nrow = 8, ncol = 1))
    }
    list(big_plot_text = big_plot_text,
         res_list_IntAct = figure$res_list_IntAct, rocr_list_IntAct = figure$rocr_list_IntAct, 
         res_list_BioPlex = figure$res_list_BioPlex, rocr_list_BioPlex = figure$rocr_list_BioPlex,
         res_list_viral = figure$res_list_viral, rocr_list_viral = figure$rocr_list_viral, 
         big_plot_IntAct = figure$big_plot_IntAct, IntAct_grob = IntAct_grob,
         big_plot_BioPlex = figure$big_plot_BioPlex, BioPlex_grob = BioPlex_grob,
         big_plot_viral = figure$big_plot_viral, viral_grob = viral_grob,
         motif_pval_cutoffs_IntAct = figure$motif_pval_cutoffs_IntAct,
         motif_pval_cutoffs_BioPlex = figure$motif_pval_cutoffs_BioPlex,
         motif_pval_cutoffs_viral = figure$motif_pval_cutoffs_viral)
}
```


```r
##################################################################################### "MOD", "LIG", "DOC"
# All motifs
start.time = Sys.time()
figureMODLIGDOC = figureVenn(filter_by_domain_data = NULL,
                             main = "human network searched for motifs present in viral proteins (not filtered by domain)",
                             motif_pval_tags = c("lenient cutoff", "optimal cutoff", "stringent cutoff"),
                             datasets = datasets,
                             motif_pval_cutoffs_viral = c(0.3, 0.1, 0.03),
                             non_query_domain_results_obj = NULL, # res_count_all
                             query_domains_only = F, min_non_query_domain_support = 0,
                             return_all = F, merge_domain_data = F, descriptions = descriptions)
figureMODLIGDOC2 = replotFigureVenn(figureMODLIGDOC,
                                    main = "")$big_plot_text
ggsave(filename = "Venn_not_filtered_MOD_LIG_DOC.svg", plot = figureMODLIGDOC2,
       device = "svg", width = 19.3, height = 11.60638, units = "in")
rm()
Sys.time() - start.time

# Filter by viral-human domain predictions
figureMODLIGDOC_filtered = figureVenn(filter_by_domain_data = "p.value < 0.5",
                                      main = "human network searched for motifs present in viral proteins (filtered by domain)",
                                      motif_pval_tags = c("lenient cutoff", "optimal cutoff", "stringent cutoff"),
                                      datasets = datasets,
                                      non_query_domain_results_obj = NULL, # res_count_all
                                      query_domains_only = F, min_non_query_domain_support = 0,
                                      return_all = F, merge_domain_data = F, descriptions = descriptions)
figureMODLIGDOC_filtered2 = replotFigureVenn(figureMODLIGDOC_filtered,
                                             main = "")$big_plot_text
ggsave(filename = "Venn_filtered_MOD_LIG_DOC.svg", plot = figureMODLIGDOC_filtered2,
       device = "svg", width = 19.3, height = 11.60638, units = "in")
save(list = ls(), file=filename)

########################################### -human domains
# Filter by human-human and viral-human domain predictions - some human proteins
start.time = Sys.time()
figureMODLIGDOC_1 = figureVenn(filter_by_domain_data = "p.value < 0.5",
                               main = "Filter by human-human and viral-human domain predictions (4 human inst. predict the same domain)",
                               motif_pval_tags = c("lenient cutoff", "optimal cutoff", "stringent cutoff"),
                               datasets = datasets,
                               non_query_domain_results_obj = "res_count_all", # res_count_all
                               query_domains_only = F, min_non_query_domain_support = 4, min_top_domain_support4motif_nq = 0,
                               return_all = F, merge_domain_data = T, descriptions = descriptions,
                               viral_only_venn = T)
figureMODLIGDOC_12 = replotFigureVenn(figureMODLIGDOC_1,
                                      main = "")$big_plot_text
ggsave(filename = "Venn_filtered_some_human_dom_MOD_LIG_DOC.svg", plot = figureMODLIGDOC_12,
       device = "svg", width = 19.3, height = 11.60638, units = "in")
save(list = ls(), file=filename)
Sys.time() - start.time

# Filter by human-human and viral-human domain predictions - many human proteins
figureMODLIGDOC_2 = figureVenn(filter_by_domain_data = "p.value < 0.5",
                               main = "Filter by human-human and viral-human domain predictions (8 human inst. predict the same domain)",
                               motif_pval_tags = c("lenient cutoff", "optimal cutoff", "stringent cutoff"),
                               datasets = datasets,
                               non_query_domain_results_obj = "res_count_all", # res_count_all
                               query_domains_only = F, min_non_query_domain_support = 8, min_top_domain_support4motif_nq = 0,
                               return_all = F, merge_domain_data = T, descriptions = descriptions,
                               viral_only_venn = T)
figureMODLIGDOC_22 = replotFigureVenn(figureMODLIGDOC_2,
                                      main = "")$big_plot_text
ggsave(filename = "Venn_filtered_many_human_dom_MOD_LIG_DOC.svg", plot = figureMODLIGDOC_22,
       device = "svg", width = 19.3, height = 11.60638, units = "in")
save(list = ls(), file=filename)

# Filter by human-human and viral-human domain predictions - both human and viral proteins predict domain
figureMODLIGDOC_3 = figureVenn(filter_by_domain_data = "p.value < 0.5",
                               main = "Filter by human-human and viral-human domain predictions (4 human inst. predict the same domain as viral instance)",
                               motif_pval_tags = c("lenient cutoff", "optimal cutoff", "stringent cutoff"),
                               datasets = datasets,
                               non_query_domain_results_obj = "res_count_all", # res_count_all
                               query_domains_only = T, min_non_query_domain_support = 4, min_top_domain_support4motif_nq = 0,
                               return_all = F, merge_domain_data = T, descriptions = descriptions,
                               viral_only_venn = T)
figureMODLIGDOC_32 = replotFigureVenn(figureMODLIGDOC_3,
                                      main = "")$big_plot_text
ggsave(filename = "Venn_filtered_some_humanviral_dom_MOD_LIG_DOC.svg", plot = figureMODLIGDOC_32,
       device = "svg", width = 19.3, height = 11.60638, units = "in")
save(list = ls(), file=filename)
```


```r
# Filter by human-human and viral-human domain predictions - (3 human have the same top-1 domain)
figureMODLIGDOC_4 = figureVenn(filter_by_domain_data = "p.value < 0.5",
                               main = "Filter by human-human and viral-human domain predictions (3 human protein have the same top-1 domain)",
                               motif_pval_tags = c("lenient cutoff", "optimal cutoff", "stringent cutoff"),
                               stringent_prec = 0.5,
                               datasets = datasets,
                               non_query_domain_results_obj = "res_count_all", # res_count_all
                               query_domains_only = F, min_non_query_domain_support = 0, min_top_domain_support4motif_nq = 3,
                               return_all = F, merge_domain_data = T, descriptions = descriptions,
                               viral_only_venn = T)
figureMODLIGDOC_42 = replotFigureVenn(figureMODLIGDOC_4,
                                      main = "")$big_plot_text
ggsave(filename = "Venn_filtered_some_human_top_dom_MOD_LIG_DOC.svg", plot = figureMODLIGDOC_42,
       device = "svg", width = 19.3, height = 11.60638, units = "in")

# Filter by human-human and viral-human domain predictions - (3 human have the same top-1 domain)
figureMODLIGDOC_5 = figureVenn(filter_by_domain_data = "p.value < 0.5",
                               main = "Filter by human-human and viral-human domain predictions (3 human protein have the same top-1 domain)",
                               motif_pval_tags = c("lenient cutoff", "optimal cutoff", "stringent cutoff"),
                               stringent_prec = 0.3,
                               datasets = datasets,
                               non_query_domain_results_obj = "res_count_all", # res_count_all
                               query_domains_only = F, min_non_query_domain_support = 0, min_top_domain_support4motif_nq = 3,
                               return_all = F, merge_domain_data = T, descriptions = descriptions,
                               viral_only_venn = T)
figureMODLIGDOC_52 = replotFigureVenn(figureMODLIGDOC_5,
                                      main = "")$big_plot_text
ggsave(filename = "Venn_str_prec_filtered_some_human_top_dom_MOD_LIG_DOC.svg", plot = figureMODLIGDOC_52,
       device = "svg", width = 19.3, height = 11.60638, units = "in")
#grid.newpage()
#grid.draw(figure_filtered2$big_plot_text)
save(list = ls(), file=filename)
```


```r
##################################################################################### ^ copy optimal and stringent motif_pval_cutoffs
##################################################################################### "MOD", "LIG", "DOC", "DEG", "CLV"
motif_types = c("MOD", "LIG", "DOC", "DEG", "CLV")

figureMODLIGDOC2DEGCLV = figureVenn(filter_by_domain_data = NULL,
                                    main = "human network searched for motifs present in viral proteins (not filtered by domain)",
                                    motif_pval_tags = c("lenient cutoff", "optimal cutoff", "stringent cutoff"),
                                    datasets = datasets,
                                    motif_pval_cutoffs_IntAct = figureMODLIGDOC$motif_pval_cutoffs_IntAct,
                                    motif_pval_cutoffs_BioPlex = figureMODLIGDOC$motif_pval_cutoffs_BioPlex,
                                    motif_pval_cutoffs_viral = figureMODLIGDOC$motif_pval_cutoffs_viral,
                                    non_query_domain_results_obj = NULL, # res_count_all
                                    query_domains_only = F, min_non_query_domain_support = 0,
                                    return_all = F, merge_domain_data = F)
figureMODLIGDOC2DEGCLV2 = replotFigureVenn(figureMODLIGDOC2DEGCLV,
                                           main = "")
ggsave(filename = "Venn_not_filtered_MOD_LIG_DOC_DEG_CLV.svg", plot = figureMODLIGDOC2DEGCLV2$big_plot_text,
       device = "svg", width = 19.3, height = 11.60638, units = "in")
#grid.newpage()
#grid.draw(figure2$big_plot_text)

figureMODLIGDOC2DEGCLV_filtered = figureVenn(filter_by_domain_data = "p.value < 0.5",
                                             main = "human network searched for motifs present in viral proteins (filtered by domain)",
                                             motif_pval_tags = c("lenient cutoff", "optimal cutoff", "stringent cutoff"),
                                             datasets = datasets,
                                             motif_pval_cutoffs_IntAct = figureMODLIGDOC2_filtered$motif_pval_cutoffs_IntAct,
                                             motif_pval_cutoffs_BioPlex = figureMODLIGDOC2_filtered$motif_pval_cutoffs_BioPlex,
                                             motif_pval_cutoffs_viral = figureMODLIGDOC2_filtered$motif_pval_cutoffs_viral,
                                             non_query_domain_results_obj = NULL, # res_count_all
                                             query_domains_only = F, min_non_query_domain_support = 0,
                                             return_all = F, merge_domain_data = F)
figureMODLIGDOC2DEGCLV_filtered2 = replotFigureVenn(figureMODLIGDOC2DEGCLV_filtered,
                                                    main = "")
ggsave(filename = "Venn_filtered_MOD_LIG_DOC_DEG_CLV.svg", plot = figureMODLIGDOC2DEGCLV_filtered2$big_plot_text,
       device = "svg", width = 19.3, height = 11.60638, units = "in")
#grid.newpage()
#grid.draw(figure_filtered2$big_plot_text)
##################################################################################### "MOD", "LIG", "DOC" "DEG", "CLV", "TRG"
motif_types = c("MOD", "LIG", "DOC", "DEG", "CLV", "TRG")

figureMODLIGDOC2DEGCLVTRG = figureVenn(filter_by_domain_data = NULL,
                                       main = "human network searched for motifs present in viral proteins (not filtered by domain)",
                                       motif_pval_tags = c("lenient cutoff", "optimal cutoff", "stringent cutoff"),
                                       datasets = datasets,
                                       motif_pval_cutoffs_IntAct = figureMODLIGDOC$motif_pval_cutoffs_IntAct,
                                       motif_pval_cutoffs_BioPlex = figureMODLIGDOC$motif_pval_cutoffs_BioPlex,
                                       motif_pval_cutoffs_viral = figureMODLIGDOC$motif_pval_cutoffs_viral,
                                       non_query_domain_results_obj = NULL, # res_count_all
                                       query_domains_only = F, min_non_query_domain_support = 0,
                                       return_all = F, merge_domain_data = F)
figureMODLIGDOC2DEGCLVTRG2 = replotFigureVenn(figureMODLIGDOC2DEGCLVTRG,
                                              main = "")
ggsave(filename = "Venn_not_filtered_MOD_LIG_DOC_DEG_CLV_TRG.svg", plot = figureMODLIGDOC2DEGCLVTRG2$big_plot_text,
       device = "svg", width = 19.3, height = 11.60638, units = "in")
#grid.newpage()
#grid.draw(figure2$big_plot_text)

figureMODLIGDOC2DEGCLVTRG_filtered = figureVenn(filter_by_domain_data = "p.value < 0.5",
                                                main = "human network searched for motifs present in viral proteins (filtered by domain)",
                                                motif_pval_tags = c("lenient cutoff", "optimal cutoff", "stringent cutoff"),
                                                datasets = datasets,
                                                motif_pval_cutoffs_IntAct = figureMODLIGDOC2_filtered$motif_pval_cutoffs_IntAct,
                                                motif_pval_cutoffs_BioPlex = figureMODLIGDOC2_filtered$motif_pval_cutoffs_BioPlex,
                                                motif_pval_cutoffs_viral = figureMODLIGDOC2_filtered$motif_pval_cutoffs_viral,
                                                non_query_domain_results_obj = NULL, # res_count_all
                                                query_domains_only = F, min_non_query_domain_support = 0,
                                                return_all = F, merge_domain_data = F)
figureMODLIGDOC2DEGCLVTRG_filtered2 = replotFigureVenn(figureMODLIGDOC2DEGCLVTRG_filtered,
                                                       main = "")
ggsave(filename = "Venn_filtered_MOD_LIG_DOC_DEG_CLV_TRG.svg", plot = figureMODLIGDOC2DEGCLVTRG_filtered2$big_plot_text,
       device = "svg", width = 19.3, height = 11.60638, units = "in")
#grid.newpage()
#grid.draw(figure_filtered2$big_plot_text)
#####################################################################################
```

# Create protein-domain-motif-protein network


```r
# choose motif types, domain and motif p-value threshold
motif_types = c("MOD", "LIG", "DOC")
# get motif & domains (stringent threshold)
IntAct = myBenchmarkMotifs(neg_set = "all_proteins",
                           filter_by_domain_data = "p.value < 0.5",
                           motif_pval_cutoff = 0.002,
                           dataset = datasets[1], datasets = datasets,
                           return_all = T, 
                           non_query_domain_results_obj = NULL, # res_count_all
                           query_domains_only = F, min_non_query_domain_support = 0, merge_domain_data = T)

figureMODLIGDOC_4$motif_pval_cutoffs_IntAct
```

```
## [1] 0.300 0.180 0.002
```

```r
figureMODLIGDOC_5$motif_pval_cutoffs_IntAct
```

```
## [1] 0.30 0.18 0.18
```

```r
IntAct_filt_hum_dom = myBenchmarkMotifs(neg_set = "all_proteins",
                                   filter_by_domain_data = "p.value < 0.5",
                                   motif_pval_cutoff = figureMODLIGDOC_4$motif_pval_cutoffs_IntAct[2],
                                   dataset = datasets[1], datasets = datasets,
                                   return_all = T, 
                                   non_query_domain_results_obj = "res_count_all", # res_count_all
                                   query_domains_only = F, min_non_query_domain_support = 0,
                                   min_top_domain_support4motif_nq = 3, merge_domain_data = T)
rocr_list_IntAct = mBenchmarkMotifsROC(IntAct_filt_hum_dom$res, measure1 = "prec", measure2 = "rec", ROCR_data = T, legend_args = "x = 0.2 | y = 0.9")
```

![](/nfs/research1/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/protein_domain_motif_protein_network-1.png)<!-- -->

```r
viral = myBenchmarkMotifs(neg_set = "all_proteins",
                          filter_by_domain_data = "p.value < 0.5",
                          motif_pval_cutoff = 0.04,
                          dataset = datasets[3], datasets = datasets,
                          return_all = T,
                          non_query_domain_results_obj = NULL, # res_count_all
                          query_domains_only = F, min_non_query_domain_support = 0, merge_domain_data = T)
save(list = ls(), file=filename)
```

# how many motif instance at a stringent threshold after domain filtering?


```r
bench_res_cyt = IntAct_filt_hum_dom$res$qslimfinder.Full_IntAct3.FALSE
instances = paste0(bench_res_cyt$occurence_query$query, "_",
                         bench_res_cyt$occurence_query$Pattern)
instances = instances[!is.na(bench_res_cyt$occurence_query$domain)]
instances = unique(instances)
proteins = sapply((strsplit(instances, "_")), function(x) x[1])
patterns = sapply((strsplit(instances, "_")), function(x) x[2])
unique(proteins)
```

```
##   [1] "B4URF7" "C5E522" "C5E526" "C5E527" "D1LN35" "E5LBT9" "F5HE15"
##   [8] "F5HFG5" "I6T1Z2" "I6TAH8" "K7Y1A2" "O39474" "O40939" "O56264"
##  [15] "O92837" "P03070" "P03087" "P03101" "P03126" "P03129" "P03177"
##  [22] "P03188" "P03209" "P03211" "P03220" "P03225" "P03243" "P03246"
##  [29] "P03255" "P03346" "P03366" "P03427" "P03428" "P03431" "P03433"
##  [36] "P03466" "P03495" "P03496" "P04012" "P04015" "P04296" "P04413"
##  [43] "P04487" "P04591" "P05919" "P06427" "P06428" "P06429" "P06430"
##  [50] "P06460" "P06462" "P06463" "P06464" "P06788" "P06821" "P06827"
##  [57] "P06930" "P08393" "P09992" "P0C1C6" "P0C1C7" "P0C213" "P0C739"
##  [64] "P0C746" "P0CK56" "P12418" "P13285" "P16717" "P17382" "P17386"
##  [71] "P21605" "P21698" "P21735" "P24772" "P24835" "P26555" "P27228"
##  [78] "P30119" "P31345" "P36780" "P50804" "Q01220" "Q04360" "Q05127"
##  [85] "Q05322" "Q0A2H0" "Q194T2" "Q1K9H2" "Q1K9H5" "Q2MG95" "Q2PJP0"
##  [92] "Q2PJP1" "Q2Q067" "Q5EP28" "Q67296" "Q69117" "Q6DP93" "Q6VGS8"
##  [99] "Q76S40" "Q77M19" "Q77Q36" "Q77UU1" "Q8AZK7" "Q997F2" "Q99AU3"
## [106] "Q9DGW5" "Q9QPN3" "Q9WMB5" "Q9WMX2" "Q9WPI5" "U5TQE9"
```

```r
unique(patterns)
```

```
##   [1] "KRKR"               "E..E..Q"            "RKR"               
##   [4] "Q..T.I"             "EVI"                "I..K.{1,2}R.{1,2}K"
##   [7] "N.N..T"             "E..D.I"             "L.R..K"            
##  [10] "KRK.E"              "R.A..E"             "KRK"               
##  [13] "L.D.D"              "EAL"                "R.RR"              
##  [16] "P..KR"              "G.G..T"             "S..AS"             
##  [19] "PA..P"              "TP[AS][AS]"         "KQK"               
##  [22] "E.EE"               "P[DE][DE]"          "SD.E"              
##  [25] "DSG"                "[ST].[LV]$"         "[DE]T.[ILV]$"      
##  [28] "E.E..I"             "P..QA"              "QA.P"              
##  [31] "EED"                "R..QE"              "E..[AGS]F"         
##  [34] "P.TP"               "[KR][KR][KR][KR]"   "SSG"               
##  [37] "APR.R"              "KR.R"               "Q..A..R"           
##  [40] "Q.EA"               "P[GS][AS]P"         "E..K.K"            
##  [43] "PV.A"               "I.KK"               "KEK.[DE]"          
##  [46] "E..E.L"             "EEL"                "T.SS"              
##  [49] "K.E..E"             "RF..[AGS]"          "K..RK.A"           
##  [52] "EE.D"               "[HK][KR].[HKR]"     "[HR][HK]R.R"       
##  [55] "[KR][KR].[KR]"      "P.S.G"              "AAG"               
##  [58] "PAP"                "E..K..K"            "DT..D"             
##  [61] "LKE.V"              "[DE]..P[LV]"        "STS"               
##  [64] "L.R..E"             "E..QR"              "[DE]..[KR].[HKR]"  
##  [67] "D..E.E"             "S.ED"               "E.DE"              
##  [70] "D.EE"               "A..A..G"            "Q..[FIMV][DE]"     
##  [73] "G.FV"               "EK.K"               "G[KR].RGR"         
##  [76] "LQ..R"              "GSS"                "S..T.P"            
##  [79] "SS..G"              "[DE]E..R"           "LS.P"              
##  [82] "R.K.E"              "KK..R"              "N.{0,1}GD"         
##  [85] "T.L$"               "D.VD"               "S..E[DE]..S"       
##  [88] "D.I..L"             "[AS]P..[AS].P"      "[DE].[ILV].[KR]"   
##  [91] "NFS"                "C..PG"              "Q..[HKR].[KR]"     
##  [94] "PP..D"              "KQ.{0,2}K.{1,2}L"   "EN.[DE]R"          
##  [97] "T..V.K"             "SRG.G..T"           "NN.N"              
## [100] "[ST].V$"            "NN.V"               "K..KR"             
## [103] "EE..Q"              "I.ED"               "[ST].[ILV]$"       
## [106] "A..S.D"             "EA..R..L"           "V.G..P.G"          
## [109] "V.S.{0,2}G"         "P..G..G"            "PA.P"              
## [112] "E.LR"               "G.G.G"              "[GS][ILM][AS].P"   
## [115] "Q..L..[DE]S"        "K..H.R..S"          "E.DD"              
## [118] "R.LQ"               "SEN"                "L..PP"             
## [121] "P..[LV]..P.L"       "S.SP[AS]"           "RG.R..G"           
## [124] "DED"                "P..PG..P"           "R..R..G"           
## [127] "P..PG"              "PG..P"              "L.Q.L.{0,1}R"      
## [130] "RK.A"               "E..R.Y"             "P..P.[HKR]"        
## [133] "AA.P"               "EE..R"              "D..P..A"
```

# do results overlap?


```r
IntAct_known = unique(paste0(seqnames(IntAct$res$qslimfinder.Full_IntAct3.FALSE$overlapping_GRanges_query), "_",
                             start(IntAct$res$qslimfinder.Full_IntAct3.FALSE$overlapping_GRanges_query), "_",
                             end(IntAct$res$qslimfinder.Full_IntAct3.FALSE$overlapping_GRanges_query)))
viral_known = unique(paste0(seqnames(viral$res$qslimfinder.all_viral_interaction3.FALSE$overlapping_GRanges_query), "_",
                            start(viral$res$qslimfinder.all_viral_interaction3.FALSE$overlapping_GRanges_query), "_",
                            end(viral$res$qslimfinder.all_viral_interaction3.FALSE$overlapping_GRanges_query)))
union_known = unique(c(IntAct_known, viral_known))

IntAct_only = union_known %in% IntAct_known
viral_only = union_known %in% viral_known

names(IntAct_only) = union_known
names(viral_only) = union_known
```


```r
# convert motif GRanges to simple network table
bechmark2Net = function(benchmarkMotifsResult = IntAct$res$qslimfinder.Full_IntAct3.FALSE,
                        top_domain_by = c(NA, "p.value", "combined_p.value", "domain_support4motif_nq", "top_domain_support4motif_nq")) {
    occurence_query = benchmarkMotifsResult$occurence_query
    occurence = benchmarkMotifsResult$occurence
    matching_known = benchmarkMotifsResult$overlapping_GRanges_query
    # check if some columns are in the result
    if(!"top_domain_support4motif_nq" %in% colnames(mcols(occurence_query))) occurence_query$top_domain_support4motif_nq = NA 
    if(!"domain_support4motif_nq" %in% colnames(mcols(occurence_query))) occurence_query$domain_support4motif_nq = NA
    if(!"domain_support4motif_q" %in% colnames(mcols(occurence_query))) occurence_query$domain_support4motif_q = NA 
    if(!"combined_p.value" %in% colnames(mcols(occurence_query))) occurence_query$combined_p.value = NA
    Net = data.table(Viral_protein = occurence_query$query,
                     Motif_Sig = occurence_query$Sig,
                     Domain_p.value = occurence_query$p.value,
                     top_domain_support4motif_nq = occurence_query$top_domain_support4motif_nq,
                     domain_support4motif_nq = occurence_query$domain_support4motif_nq,
                     domain_support4motif_q = occurence_query$domain_support4motif_q,
                     Domain_combined_p.value = occurence_query$combined_p.value,
                     Pattern = occurence_query$Pattern,
                     Human_recognition_domain = occurence_query$domain,
                     Human_protein = occurence_query$interacts_with,
                     unique_position = paste0(GenomicRanges::seqnames(occurence_query), "_",
                                                         GenomicRanges::start(occurence_query), "_",
                                                         GenomicRanges::end(occurence_query)))
    Net = unique(Net)
    if(!is.na(top_domain_by[1])){
        if(top_domain_by[1] == "p.value"){
            Net[, top_domain := Domain_p.value ==
                    min(Domain_p.value), by = unique_position]
        }
        if(top_domain_by[1] == "combined_p.value") {
            Net[, top_domain := Domain_combined_p.value ==
                    min(Domain_combined_p.value), by = unique_position]
        }
        if(top_domain_by[1] == "domain_support4motif_nq") {
            Net[, top_domain := domain_support4motif_nq ==
                    max(domain_support4motif_nq), by = unique_position]
        }
        if(top_domain_by[1] == "top_domain_support4motif_nq") {
            Net[, top_domain := top_domain_support4motif_nq ==
                    max(top_domain_support4motif_nq), by = unique_position]
        }
        Net = Net[top_domain == TRUE]
        Net[, top_domain := NULL]
    }
    
    matching_known = data.table(unique_position = paste0(GenomicRanges::seqnames(matching_known), "_",
                                                         GenomicRanges::start(matching_known), "_",
                                                         GenomicRanges::end(matching_known)),
                                matching_known = "matching_known")
    Net = merge(Net, matching_known, by = "unique_position", all.x = T, all.y = F, allow.cartesian=TRUE)
    Net[is.na(Net$matching_known), matching_known := ""]
    
    # add non-viral proteins
    occurence$unique_position = paste0(GenomicRanges::seqnames(occurence), "_",
                                                         GenomicRanges::start(occurence), "_",
                                                         GenomicRanges::end(occurence))
    occurence = unique(as.data.table(occurence)[query != prot_names, .(Pattern, query, interacts_with, unique_position)])
    occurence[, human_instances := paste0("when query interacts_with: ", unique(interacts_with),
                                          " supporting human_instances: ", paste0(unique(unique_position), collapse = " | ")),
              by = .(Pattern, query, interacts_with)]
    occurence[, unique_position := NULL]
    occurence = unique(occurence)
    occurence[, human_instances := paste0("query: ",unique(query), ", Pattern:", Pattern,
                                          ", interacts_with: (", paste0(unique(interacts_with), collapse = " | "),
                                          "), human_instances: (", paste0(unique(human_instances), collapse = " || "),")"),
              by = .(Pattern, query)]
    
    Net = merge(Net, unique(occurence), 
                          by.x = c("Pattern", "Viral_protein", "Human_protein"),
                          by.y = c("Pattern", "query", "interacts_with"),
                          all.x = T, all.y = F) 
    
    unique(Net)
}

IntAct_net = bechmark2Net(IntAct$res$qslimfinder.Full_IntAct3.FALSE)
IntAct_filt_hum_dom_net = bechmark2Net(IntAct_filt_hum_dom$res$qslimfinder.Full_IntAct3.FALSE,
                                       top_domain_by = "top_domain_support4motif_nq")
viral_net = bechmark2Net(viral$res$qslimfinder.all_viral_interaction3.FALSE)

# which known
IntAct_filt_hum_dom_net[matching_known == "matching_known"]
```

```
##          Pattern Viral_protein Human_protein unique_position Motif_Sig
##  1:      A..A..G        P06821        P50570    P06821_83_89     0.110
##  2:       E.E..I        P03129        P24278    P03129_33_38     0.002
##  3:          EED        P03129        P62191    P03129_34_36     0.071
##  4:         G.FV        P06821        Q8N448    P06821_89_92     0.043
##  5:       L.R..E        P06463        P17980  P06463_150_155     0.067
##  6:         SD.E        P03070        P06493  P03070_112_115     0.091
##  7:         SD.E        P03070        P06493  P03070_112_115     0.091
##  8: [DE]T.[ILV]$        P03126        Q14160  P03126_155_158     0.000
##  9: [DE]T.[ILV]$        P06463        Q14160  P06463_155_158     0.000
## 10:   [ST].[LV]$        P03126        Q12959  P03126_156_158     0.000
## 11:   [ST].[LV]$        P06463        Q12959  P06463_156_158     0.000
##     Domain_p.value top_domain_support4motif_nq domain_support4motif_nq
##  1:   1.764886e-03                           3                      16
##  2:   1.481115e-02                           3                       8
##  3:   4.752209e-03                           5                      17
##  4:   1.655336e-01                           3                       6
##  5:   7.646704e-06                           9                      14
##  6:   1.684202e-02                           5                      18
##  7:   1.684202e-02                           5                      19
##  8:   1.865102e-05                           6                      11
##  9:   1.591573e-04                           6                      11
## 10:   1.865102e-05                          13                      18
## 11:   1.591573e-04                          13                      18
##     domain_support4motif_q Domain_combined_p.value
##  1:                      2             0.008919226
##  2:                      2             0.006152263
##  3:                      1             0.004028852
##  4:                      1             0.065627800
##  5:                      3             0.006787246
##  6:                      1             0.006242212
##  7:                      1             0.006730409
##  8:                      8             0.002419502
##  9:                      8             0.002419502
## 10:                     10             0.002134549
## 11:                     10             0.002134549
##     Human_recognition_domain matching_known
##  1:                IPR027417 matching_known
##  2:                IPR013083 matching_known
##  3:                IPR027417 matching_known
##  4:                IPR013083 matching_known
##  5:                IPR027417 matching_known
##  6:                IPR000719 matching_known
##  7:                IPR011009 matching_known
##  8:                IPR001478 matching_known
##  9:                IPR001478 matching_known
## 10:                IPR001478 matching_known
## 11:                IPR001478 matching_known
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                human_instances
##  1:                                                                                                                                                                                                     query: P06821, Pattern:A..A..G, interacts_with: (P50570), human_instances: (when query interacts_with: P50570 supporting human_instances: A0JLT2_48_54 | O14576_571_577 | O14645_179_185 | O95817_529_535 | P00519_988_994 | P00519_923_929 | P03407_27_33 | P29353_436_442 | P29597_931_937 | P67809_315_321 | Q05193_805_811 | Q05516_312_318 | Q05516_135_141 | Q2MV58_40_46 | Q7Z3C6_591_597 | Q7Z3C6_679_685 | Q8NFF5-2_111_117 | Q8WV41_91_97 | Q9NX70_191_197 | Q9Y5X1_186_192)
##  2:                                                                                                                                                                                                                                                                                                                                                                                                      query: P03129, Pattern:E.E..I, interacts_with: (P24278), human_instances: (when query interacts_with: P24278 supporting human_instances: O75506_65_70 | P06821_6_11 | P11021_71_76 | P45984_344_349 | P55072_319_324 | Q7L273_66_71 | Q92609_16_21 | Q9BQD3_111_116 | Q9BRK4_606_611)
##  3:                                                                                                                                                         query: P03129, Pattern:EED, interacts_with: (P62191), human_instances: (when query interacts_with: P62191 supporting human_instances: I3L2W2_903_905 | O00487_261_263 | O95714_1829_1831 | P01106_416_418 | P17980_105_107 | P29692_166_168 | P46108-2_148_150 | P55036_373_375 | P55036_322_324 | Q05086-3_201_203 | Q05086-3_391_393 | Q13200_626_628 | Q13200_47_49 | Q13573_495_497 | Q15652_356_358 | Q15652_1665_1667 | Q15652_1560_1562 | Q16186_401_403 | Q96BD8_78_80 | Q96DT6_435_437 | Q9UBN6_375_377 | Q9UKE5_665_667)
##  4:                                                                                                                                                                                                                                                                                                                                                                                                                                                    query: P06821, Pattern:G.FV, interacts_with: (Q8N448), human_instances: (when query interacts_with: Q8N448 supporting human_instances: O43447_77_80 | P11802_29_32 | Q8WVD3_220_223 | Q92530_189_192 | Q9H832_118_121 | Q9UH92-3_27_30)
##  5:                                                                                                                                                                                                                   query: P06463, Pattern:L.R..E, interacts_with: (P17980), human_instances: (when query interacts_with: P17980 supporting human_instances: O00233_26_31 | O75832_220_225 | P03255_95_100 | P03255-2_95_100 | P35998_337_342 | P43686_327_332 | Q13573_132_137 | Q15008_81_86 | Q16401_380_385 | Q567U6_3_8 | Q5HYA8_809_814 | Q99459_545_550 | Q99459_267_272 | Q9BUZ4_295_300 | Q9H1R2-3_146_151 | Q9Y2J4_331_336 | Q9Y2J4_324_329 | Q9Y2J4-4_389_394 | Q9Y2J4-4_382_387)
##  6:                                                                                                                                                                              query: P03070, Pattern:SD.E, interacts_with: (P06493), human_instances: (when query interacts_with: P06493 supporting human_instances: O75164_616_619 | O75164_523_526 | O95067_11_14 | P01106_250_253 | P01106_348_351 | P19838_337_340 | P29590_565_568 | P30153_122_125 | P50991_158_161 | P51858_133_136 | P62136_2_5 | Q16543_13_16 | Q16667_15_18 | Q8N3C0_1850_1853 | Q8TDM6-4_3_6 | Q96EB6_162_165 | Q96EB6_659_662 | Q96FW1_16_19 | Q96PU4_643_646 | Q96RL1_46_49 | Q9H4X1-2_55_58 | Q9NQS7_214_217)
##  7:                                                                                                                                                                              query: P03070, Pattern:SD.E, interacts_with: (P06493), human_instances: (when query interacts_with: P06493 supporting human_instances: O75164_616_619 | O75164_523_526 | O95067_11_14 | P01106_250_253 | P01106_348_351 | P19838_337_340 | P29590_565_568 | P30153_122_125 | P50991_158_161 | P51858_133_136 | P62136_2_5 | Q16543_13_16 | Q16667_15_18 | Q8N3C0_1850_1853 | Q8TDM6-4_3_6 | Q96EB6_162_165 | Q96EB6_659_662 | Q96FW1_16_19 | Q96PU4_643_646 | Q96RL1_46_49 | Q9H4X1-2_55_58 | Q9NQS7_214_217)
##  8:                                                                                                                                                                                                                        query: P03126, Pattern:[DE]T.[ILV]$, interacts_with: (Q14160), human_instances: (when query interacts_with: Q14160 supporting human_instances: B7Z2Y1_449_452 | O60346_1714_1717 | P06427_146_149 | P06463_155_158 | P0C213_350_353 | P21735_155_158 | P22460_610_613 | P23508_826_829 | P24835_155_158 | P27228_146_149 | P33402_729_732 | P35222_778_781 | P50804_155_158 | P53778_364_367 | Q15311_652_655 | Q7Z628_593_596 | Q9NYB5_709_712 | Q9UDY2_1187_1190)
##  9:                                                                                                                                                                                                                        query: P06463, Pattern:[DE]T.[ILV]$, interacts_with: (Q14160), human_instances: (when query interacts_with: Q14160 supporting human_instances: B7Z2Y1_449_452 | O60346_1714_1717 | P03126_155_158 | P06427_146_149 | P0C213_350_353 | P21735_155_158 | P22460_610_613 | P23508_826_829 | P24835_155_158 | P27228_146_149 | P33402_729_732 | P35222_778_781 | P50804_155_158 | P53778_364_367 | Q15311_652_655 | Q7Z628_593_596 | Q9NYB5_709_712 | Q9UDY2_1187_1190)
## 10: query: P03126, Pattern:[ST].[LV]$, interacts_with: (Q12959), human_instances: (when query interacts_with: Q12959 supporting human_instances: B7Z2Y1_450_452 | O60333-3_1151_1153 | P06427_147_149 | P06463_156_158 | P0C213_351_353 | P16717_146_148 | P17386_147_149 | P21735_156_158 | P22459_651_653 | P22460_611_613 | P24835_156_158 | P26555_147_149 | P27228_147_149 | P50804_156_158 | P53778_365_367 | Q01814_1241_1243 | Q13224_1482_1484 | Q14524_2014_2016 | Q15303_1306_1308 | Q6VGS8_126_128 | Q7Z628_594_596 | Q96A65_972_974 | Q96GG9_257_259 | Q96PE1_1336_1338 | Q99569_1190_1192 | Q9NS75_344_346 | Q9NVW2_622_624 | Q9NYB5_710_712 | Q9P021_99_101 | Q9UQB3_1223_1225)
## 11: query: P06463, Pattern:[ST].[LV]$, interacts_with: (Q12959), human_instances: (when query interacts_with: Q12959 supporting human_instances: B7Z2Y1_450_452 | O60333-3_1151_1153 | P03126_156_158 | P06427_147_149 | P0C213_351_353 | P16717_146_148 | P17386_147_149 | P21735_156_158 | P22459_651_653 | P22460_611_613 | P24835_156_158 | P26555_147_149 | P27228_147_149 | P50804_156_158 | P53778_365_367 | Q01814_1241_1243 | Q13224_1482_1484 | Q14524_2014_2016 | Q15303_1306_1308 | Q6VGS8_126_128 | Q7Z628_594_596 | Q96A65_972_974 | Q96GG9_257_259 | Q96PE1_1336_1338 | Q99569_1190_1192 | Q9NS75_344_346 | Q9NVW2_622_624 | Q9NYB5_710_712 | Q9P021_99_101 | Q9UQB3_1223_1225)
```

```r
# human network missed
viral_net[unique_position == "P03129_21_26"]
```

```
##     Pattern Viral_protein Human_protein unique_position Motif_Sig
##  1:  DL.C.E        P03129        P28749    P03129_21_26     0.020
##  2:  DL.C.E        P03129        P28749    P03129_21_26     0.020
##  3:  DL.C.E        P03129        P28749    P03129_21_26     0.020
##  4:  DL.C.E        P03129        P28749    P03129_21_26     0.020
##  5:  DL.C.E        P03129        P28749    P03129_21_26     0.020
##  6:  DL.C.E        P03129        Q08999    P03129_21_26     0.013
##  7:  DL.C.E        P03129        Q08999    P03129_21_26     0.013
##  8:  DL.C.E        P03129        Q08999    P03129_21_26     0.013
##  9:  DL.C.E        P03129        Q08999    P03129_21_26     0.013
## 10:  DL.C.E        P03129        Q08999    P03129_21_26     0.013
##     Domain_p.value top_domain_support4motif_nq domain_support4motif_nq
##  1:    0.059906868                          NA                      NA
##  2:    0.059906868                          NA                      NA
##  3:    0.059906868                          NA                      NA
##  4:    0.004752209                          NA                      NA
##  5:    0.059906868                          NA                      NA
##  6:    0.059906868                          NA                      NA
##  7:    0.004752209                          NA                      NA
##  8:    0.059906868                          NA                      NA
##  9:    0.059906868                          NA                      NA
## 10:    0.059906868                          NA                      NA
##     domain_support4motif_q Domain_combined_p.value
##  1:                     NA                      NA
##  2:                     NA                      NA
##  3:                     NA                      NA
##  4:                     NA                      NA
##  5:                     NA                      NA
##  6:                     NA                      NA
##  7:                     NA                      NA
##  8:                     NA                      NA
##  9:                     NA                      NA
## 10:                     NA                      NA
##     Human_recognition_domain matching_known
##  1:                IPR024599 matching_known
##  2:                IPR002720 matching_known
##  3:                IPR002719 matching_known
##  4:                IPR013763 matching_known
##  5:                IPR015030 matching_known
##  6:                IPR002719 matching_known
##  7:                IPR013763 matching_known
##  8:                IPR002720 matching_known
##  9:                IPR015030 matching_known
## 10:                IPR024599 matching_known
##                                                                                                                                                                                                                                                                                                                                                                                      human_instances
##  1: query: P03129, Pattern:DL.C.E, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06430_26_31 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  2: query: P03129, Pattern:DL.C.E, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06430_26_31 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  3: query: P03129, Pattern:DL.C.E, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06430_26_31 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  4: query: P03129, Pattern:DL.C.E, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06430_26_31 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  5: query: P03129, Pattern:DL.C.E, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06430_26_31 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  6: query: P03129, Pattern:DL.C.E, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06430_26_31 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  7: query: P03129, Pattern:DL.C.E, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06430_26_31 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  8: query: P03129, Pattern:DL.C.E, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06430_26_31 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  9: query: P03129, Pattern:DL.C.E, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06430_26_31 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
## 10: query: P03129, Pattern:DL.C.E, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06430_26_31 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: A0MPS7_21_26 | P03255_121_126 | P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
```

```r
viral_net[unique_position == "P03129_32_36"]
```

```
##        Pattern Viral_protein Human_protein unique_position Motif_Sig
## 1: SEE.{0,2}ED        P03129        Q8N7W2    P03129_32_36     0.031
##    Domain_p.value top_domain_support4motif_nq domain_support4motif_nq
## 1:             NA                          NA                      NA
##    domain_support4motif_q Domain_combined_p.value Human_recognition_domain
## 1:                     NA                      NA                     <NA>
##    matching_known
## 1: matching_known
##                                                                                                                                                               human_instances
## 1: query: P03129, Pattern:SEE.{0,2}ED, interacts_with: (Q8N7W2), human_instances: (when query interacts_with: Q8N7W2 supporting human_instances: Q04360_57_61 | Q04360_57_63)
```

```r
viral_net[unique_position == "P03255_121_126"]
```

```
##     Pattern Viral_protein Human_protein unique_position Motif_Sig
##  1:  DL.CHE        P03255        P28749  P03255_121_126  0.000841
##  2:  DL.CHE        P03255        P28749  P03255_121_126  0.000841
##  3:  DL.CHE        P03255        P28749  P03255_121_126  0.000841
##  4:  DL.CHE        P03255        P28749  P03255_121_126  0.000841
##  5:  DL.CHE        P03255        P28749  P03255_121_126  0.000841
##  6:  DL.CHE        P03255        Q08999  P03255_121_126  0.000415
##  7:  DL.CHE        P03255        Q08999  P03255_121_126  0.000415
##  8:  DL.CHE        P03255        Q08999  P03255_121_126  0.000415
##  9:  DL.CHE        P03255        Q08999  P03255_121_126  0.000415
## 10:  DL.CHE        P03255        Q08999  P03255_121_126  0.000415
##     Domain_p.value top_domain_support4motif_nq domain_support4motif_nq
##  1:    0.018249255                          NA                      NA
##  2:    0.018249255                          NA                      NA
##  3:    0.018249255                          NA                      NA
##  4:    0.004741632                          NA                      NA
##  5:    0.018249255                          NA                      NA
##  6:    0.004741632                          NA                      NA
##  7:    0.018249255                          NA                      NA
##  8:    0.018249255                          NA                      NA
##  9:    0.018249255                          NA                      NA
## 10:    0.018249255                          NA                      NA
##     domain_support4motif_q Domain_combined_p.value
##  1:                     NA                      NA
##  2:                     NA                      NA
##  3:                     NA                      NA
##  4:                     NA                      NA
##  5:                     NA                      NA
##  6:                     NA                      NA
##  7:                     NA                      NA
##  8:                     NA                      NA
##  9:                     NA                      NA
## 10:                     NA                      NA
##     Human_recognition_domain matching_known
##  1:                IPR002720 matching_known
##  2:                IPR015030 matching_known
##  3:                IPR002719 matching_known
##  4:                IPR013763 matching_known
##  5:                IPR024599 matching_known
##  6:                IPR013763 matching_known
##  7:                IPR002719 matching_known
##  8:                IPR002720 matching_known
##  9:                IPR015030 matching_known
## 10:                IPR024599 matching_known
##                                                                                                                                                                                                                                                                                                       human_instances
##  1: query: P03255, Pattern:DL.CHE, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  2: query: P03255, Pattern:DL.CHE, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  3: query: P03255, Pattern:DL.CHE, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  4: query: P03255, Pattern:DL.CHE, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  5: query: P03255, Pattern:DL.CHE, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  6: query: P03255, Pattern:DL.CHE, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  7: query: P03255, Pattern:DL.CHE, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  8: query: P03255, Pattern:DL.CHE, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
##  9: query: P03255, Pattern:DL.CHE, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
## 10: query: P03255, Pattern:DL.CHE, interacts_with: (P28749 | Q08999), human_instances: (when query interacts_with: P28749 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109 || when query interacts_with: Q08999 supporting human_instances: P06788_24_29 | Q2PRM9_104_109 | Q9DUG7_104_109)
```

```r
# viral network missed
IntAct_net[unique_position == "P03428_736_739"]
```

```
##    Pattern Viral_protein Human_protein unique_position Motif_Sig
## 1:    KR.R        P03428        Q92997  P03428_736_739  0.000141
## 2:    KR.R        P03428        Q92997  P03428_736_739  0.000141
## 3:    KR.R        P03428        Q92997  P03428_736_739  0.000141
## 4:    KR.R        P03428        Q92997  P03428_736_739  0.000141
## 5:    KR.R        P03428        Q92997  P03428_736_739  0.000141
## 6:    KR.R        P03428        Q92997  P03428_736_739  0.000141
## 7:    KR.R        P03428        Q92997  P03428_736_739  0.000141
##    Domain_p.value top_domain_support4motif_nq domain_support4motif_nq
## 1:     0.19871977                          NA                      NA
## 2:     0.04574882                          NA                      NA
## 3:     0.19871977                          NA                      NA
## 4:     0.19871977                          NA                      NA
## 5:     0.19871977                          NA                      NA
## 6:     0.19871977                          NA                      NA
## 7:     0.08458042                          NA                      NA
##    domain_support4motif_q Domain_combined_p.value Human_recognition_domain
## 1:                     NA                      NA                IPR001158
## 2:                     NA                      NA                IPR001478
## 3:                     NA                      NA                IPR000591
## 4:                     NA                      NA                IPR011991
## 5:                     NA                      NA                IPR003351
## 6:                     NA                      NA                IPR024580
## 7:                     NA                      NA                IPR029071
##    matching_known
## 1:               
## 2:               
## 3:               
## 4:               
## 5:               
## 6:               
## 7:               
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              human_instances
## 1: query: P03428, Pattern:KR.R, interacts_with: (Q92997), human_instances: (when query interacts_with: Q92997 supporting human_instances: O00339_125_128 | O15169_414_417 | O15481_14_17 | O43159_131_134 | O43167_159_162 | O43167_161_164 | O43463_394_397 | O43474_416_419 | O75689_138_141 | O95251_171_174 | P03431_669_672 | P0DMU9_29_32 | P0DMU9_17_20 | P20719_196_199 | P23025_31_34 | P49759-3_178_181 | P61289_243_246 | P78337_146_149 | P86480_8_11 | Q00444_156_159 | Q14498_72_75 | Q14498_56_59 | Q15149_1951_1954 | Q15149_4262_4265 | Q3V6T2_839_842 | Q5S007_947_950 | Q5SQQ9-2_101_104 | Q66GS9_289_292 | Q67296_752_755 | Q67296_736_739 | Q6ZU80_141_144 | Q7Z7J5_155_158 | Q86YD7_157_160 | Q8IX06_226_229 | Q8N7W2-2_294_297 | Q8N8U2_105_108 | Q8N960_689_692 | Q8NE31_526_529 | Q8NHU0_17_20 | Q8NHU0_29_32 | Q8TAD8_257_260 | Q8TAD8_95_98 | Q8TBE0_169_172 | Q96C24_34_37 | Q96EB6_35_38 | Q96EG3_232_235 | Q96KM6_843_846 | Q96MU7_302_305 | Q96PV6-2_489_492 | Q96T60_138_141 | Q9BTC0_1214_1217 | Q9H0I2_113_116 | Q9HAJ7_87_90 | Q9HBE1-4_347_350 | Q9HBE1-4_264_267 | Q9P127_230_233 | Q9UBU8-2_136_139 | Q9UFB7_648_651 | Q9UK58_384_387)
## 2: query: P03428, Pattern:KR.R, interacts_with: (Q92997), human_instances: (when query interacts_with: Q92997 supporting human_instances: O00339_125_128 | O15169_414_417 | O15481_14_17 | O43159_131_134 | O43167_159_162 | O43167_161_164 | O43463_394_397 | O43474_416_419 | O75689_138_141 | O95251_171_174 | P03431_669_672 | P0DMU9_29_32 | P0DMU9_17_20 | P20719_196_199 | P23025_31_34 | P49759-3_178_181 | P61289_243_246 | P78337_146_149 | P86480_8_11 | Q00444_156_159 | Q14498_72_75 | Q14498_56_59 | Q15149_1951_1954 | Q15149_4262_4265 | Q3V6T2_839_842 | Q5S007_947_950 | Q5SQQ9-2_101_104 | Q66GS9_289_292 | Q67296_752_755 | Q67296_736_739 | Q6ZU80_141_144 | Q7Z7J5_155_158 | Q86YD7_157_160 | Q8IX06_226_229 | Q8N7W2-2_294_297 | Q8N8U2_105_108 | Q8N960_689_692 | Q8NE31_526_529 | Q8NHU0_17_20 | Q8NHU0_29_32 | Q8TAD8_257_260 | Q8TAD8_95_98 | Q8TBE0_169_172 | Q96C24_34_37 | Q96EB6_35_38 | Q96EG3_232_235 | Q96KM6_843_846 | Q96MU7_302_305 | Q96PV6-2_489_492 | Q96T60_138_141 | Q9BTC0_1214_1217 | Q9H0I2_113_116 | Q9HAJ7_87_90 | Q9HBE1-4_347_350 | Q9HBE1-4_264_267 | Q9P127_230_233 | Q9UBU8-2_136_139 | Q9UFB7_648_651 | Q9UK58_384_387)
## 3: query: P03428, Pattern:KR.R, interacts_with: (Q92997), human_instances: (when query interacts_with: Q92997 supporting human_instances: O00339_125_128 | O15169_414_417 | O15481_14_17 | O43159_131_134 | O43167_159_162 | O43167_161_164 | O43463_394_397 | O43474_416_419 | O75689_138_141 | O95251_171_174 | P03431_669_672 | P0DMU9_29_32 | P0DMU9_17_20 | P20719_196_199 | P23025_31_34 | P49759-3_178_181 | P61289_243_246 | P78337_146_149 | P86480_8_11 | Q00444_156_159 | Q14498_72_75 | Q14498_56_59 | Q15149_1951_1954 | Q15149_4262_4265 | Q3V6T2_839_842 | Q5S007_947_950 | Q5SQQ9-2_101_104 | Q66GS9_289_292 | Q67296_752_755 | Q67296_736_739 | Q6ZU80_141_144 | Q7Z7J5_155_158 | Q86YD7_157_160 | Q8IX06_226_229 | Q8N7W2-2_294_297 | Q8N8U2_105_108 | Q8N960_689_692 | Q8NE31_526_529 | Q8NHU0_17_20 | Q8NHU0_29_32 | Q8TAD8_257_260 | Q8TAD8_95_98 | Q8TBE0_169_172 | Q96C24_34_37 | Q96EB6_35_38 | Q96EG3_232_235 | Q96KM6_843_846 | Q96MU7_302_305 | Q96PV6-2_489_492 | Q96T60_138_141 | Q9BTC0_1214_1217 | Q9H0I2_113_116 | Q9HAJ7_87_90 | Q9HBE1-4_347_350 | Q9HBE1-4_264_267 | Q9P127_230_233 | Q9UBU8-2_136_139 | Q9UFB7_648_651 | Q9UK58_384_387)
## 4: query: P03428, Pattern:KR.R, interacts_with: (Q92997), human_instances: (when query interacts_with: Q92997 supporting human_instances: O00339_125_128 | O15169_414_417 | O15481_14_17 | O43159_131_134 | O43167_159_162 | O43167_161_164 | O43463_394_397 | O43474_416_419 | O75689_138_141 | O95251_171_174 | P03431_669_672 | P0DMU9_29_32 | P0DMU9_17_20 | P20719_196_199 | P23025_31_34 | P49759-3_178_181 | P61289_243_246 | P78337_146_149 | P86480_8_11 | Q00444_156_159 | Q14498_72_75 | Q14498_56_59 | Q15149_1951_1954 | Q15149_4262_4265 | Q3V6T2_839_842 | Q5S007_947_950 | Q5SQQ9-2_101_104 | Q66GS9_289_292 | Q67296_752_755 | Q67296_736_739 | Q6ZU80_141_144 | Q7Z7J5_155_158 | Q86YD7_157_160 | Q8IX06_226_229 | Q8N7W2-2_294_297 | Q8N8U2_105_108 | Q8N960_689_692 | Q8NE31_526_529 | Q8NHU0_17_20 | Q8NHU0_29_32 | Q8TAD8_257_260 | Q8TAD8_95_98 | Q8TBE0_169_172 | Q96C24_34_37 | Q96EB6_35_38 | Q96EG3_232_235 | Q96KM6_843_846 | Q96MU7_302_305 | Q96PV6-2_489_492 | Q96T60_138_141 | Q9BTC0_1214_1217 | Q9H0I2_113_116 | Q9HAJ7_87_90 | Q9HBE1-4_347_350 | Q9HBE1-4_264_267 | Q9P127_230_233 | Q9UBU8-2_136_139 | Q9UFB7_648_651 | Q9UK58_384_387)
## 5: query: P03428, Pattern:KR.R, interacts_with: (Q92997), human_instances: (when query interacts_with: Q92997 supporting human_instances: O00339_125_128 | O15169_414_417 | O15481_14_17 | O43159_131_134 | O43167_159_162 | O43167_161_164 | O43463_394_397 | O43474_416_419 | O75689_138_141 | O95251_171_174 | P03431_669_672 | P0DMU9_29_32 | P0DMU9_17_20 | P20719_196_199 | P23025_31_34 | P49759-3_178_181 | P61289_243_246 | P78337_146_149 | P86480_8_11 | Q00444_156_159 | Q14498_72_75 | Q14498_56_59 | Q15149_1951_1954 | Q15149_4262_4265 | Q3V6T2_839_842 | Q5S007_947_950 | Q5SQQ9-2_101_104 | Q66GS9_289_292 | Q67296_752_755 | Q67296_736_739 | Q6ZU80_141_144 | Q7Z7J5_155_158 | Q86YD7_157_160 | Q8IX06_226_229 | Q8N7W2-2_294_297 | Q8N8U2_105_108 | Q8N960_689_692 | Q8NE31_526_529 | Q8NHU0_17_20 | Q8NHU0_29_32 | Q8TAD8_257_260 | Q8TAD8_95_98 | Q8TBE0_169_172 | Q96C24_34_37 | Q96EB6_35_38 | Q96EG3_232_235 | Q96KM6_843_846 | Q96MU7_302_305 | Q96PV6-2_489_492 | Q96T60_138_141 | Q9BTC0_1214_1217 | Q9H0I2_113_116 | Q9HAJ7_87_90 | Q9HBE1-4_347_350 | Q9HBE1-4_264_267 | Q9P127_230_233 | Q9UBU8-2_136_139 | Q9UFB7_648_651 | Q9UK58_384_387)
## 6: query: P03428, Pattern:KR.R, interacts_with: (Q92997), human_instances: (when query interacts_with: Q92997 supporting human_instances: O00339_125_128 | O15169_414_417 | O15481_14_17 | O43159_131_134 | O43167_159_162 | O43167_161_164 | O43463_394_397 | O43474_416_419 | O75689_138_141 | O95251_171_174 | P03431_669_672 | P0DMU9_29_32 | P0DMU9_17_20 | P20719_196_199 | P23025_31_34 | P49759-3_178_181 | P61289_243_246 | P78337_146_149 | P86480_8_11 | Q00444_156_159 | Q14498_72_75 | Q14498_56_59 | Q15149_1951_1954 | Q15149_4262_4265 | Q3V6T2_839_842 | Q5S007_947_950 | Q5SQQ9-2_101_104 | Q66GS9_289_292 | Q67296_752_755 | Q67296_736_739 | Q6ZU80_141_144 | Q7Z7J5_155_158 | Q86YD7_157_160 | Q8IX06_226_229 | Q8N7W2-2_294_297 | Q8N8U2_105_108 | Q8N960_689_692 | Q8NE31_526_529 | Q8NHU0_17_20 | Q8NHU0_29_32 | Q8TAD8_257_260 | Q8TAD8_95_98 | Q8TBE0_169_172 | Q96C24_34_37 | Q96EB6_35_38 | Q96EG3_232_235 | Q96KM6_843_846 | Q96MU7_302_305 | Q96PV6-2_489_492 | Q96T60_138_141 | Q9BTC0_1214_1217 | Q9H0I2_113_116 | Q9HAJ7_87_90 | Q9HBE1-4_347_350 | Q9HBE1-4_264_267 | Q9P127_230_233 | Q9UBU8-2_136_139 | Q9UFB7_648_651 | Q9UK58_384_387)
## 7: query: P03428, Pattern:KR.R, interacts_with: (Q92997), human_instances: (when query interacts_with: Q92997 supporting human_instances: O00339_125_128 | O15169_414_417 | O15481_14_17 | O43159_131_134 | O43167_159_162 | O43167_161_164 | O43463_394_397 | O43474_416_419 | O75689_138_141 | O95251_171_174 | P03431_669_672 | P0DMU9_29_32 | P0DMU9_17_20 | P20719_196_199 | P23025_31_34 | P49759-3_178_181 | P61289_243_246 | P78337_146_149 | P86480_8_11 | Q00444_156_159 | Q14498_72_75 | Q14498_56_59 | Q15149_1951_1954 | Q15149_4262_4265 | Q3V6T2_839_842 | Q5S007_947_950 | Q5SQQ9-2_101_104 | Q66GS9_289_292 | Q67296_752_755 | Q67296_736_739 | Q6ZU80_141_144 | Q7Z7J5_155_158 | Q86YD7_157_160 | Q8IX06_226_229 | Q8N7W2-2_294_297 | Q8N8U2_105_108 | Q8N960_689_692 | Q8NE31_526_529 | Q8NHU0_17_20 | Q8NHU0_29_32 | Q8TAD8_257_260 | Q8TAD8_95_98 | Q8TBE0_169_172 | Q96C24_34_37 | Q96EB6_35_38 | Q96EG3_232_235 | Q96KM6_843_846 | Q96MU7_302_305 | Q96PV6-2_489_492 | Q96T60_138_141 | Q9BTC0_1214_1217 | Q9H0I2_113_116 | Q9HAJ7_87_90 | Q9HBE1-4_347_350 | Q9HBE1-4_264_267 | Q9P127_230_233 | Q9UBU8-2_136_139 | Q9UFB7_648_651 | Q9UK58_384_387)
```

```r
IntAct_net[unique_position == "P03428_752_755"]
```

```
##    Pattern Viral_protein Human_protein unique_position Motif_Sig
## 1:    KR.R        P03428        Q92997  P03428_752_755  0.000141
## 2:    KR.R        P03428        Q92997  P03428_752_755  0.000141
## 3:    KR.R        P03428        Q92997  P03428_752_755  0.000141
## 4:    KR.R        P03428        Q92997  P03428_752_755  0.000141
## 5:    KR.R        P03428        Q92997  P03428_752_755  0.000141
## 6:    KR.R        P03428        Q92997  P03428_752_755  0.000141
## 7:    KR.R        P03428        Q92997  P03428_752_755  0.000141
##    Domain_p.value top_domain_support4motif_nq domain_support4motif_nq
## 1:     0.19871977                          NA                      NA
## 2:     0.04574882                          NA                      NA
## 3:     0.19871977                          NA                      NA
## 4:     0.19871977                          NA                      NA
## 5:     0.19871977                          NA                      NA
## 6:     0.19871977                          NA                      NA
## 7:     0.08458042                          NA                      NA
##    domain_support4motif_q Domain_combined_p.value Human_recognition_domain
## 1:                     NA                      NA                IPR001158
## 2:                     NA                      NA                IPR001478
## 3:                     NA                      NA                IPR000591
## 4:                     NA                      NA                IPR011991
## 5:                     NA                      NA                IPR003351
## 6:                     NA                      NA                IPR024580
## 7:                     NA                      NA                IPR029071
##    matching_known
## 1:               
## 2:               
## 3:               
## 4:               
## 5:               
## 6:               
## 7:               
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              human_instances
## 1: query: P03428, Pattern:KR.R, interacts_with: (Q92997), human_instances: (when query interacts_with: Q92997 supporting human_instances: O00339_125_128 | O15169_414_417 | O15481_14_17 | O43159_131_134 | O43167_159_162 | O43167_161_164 | O43463_394_397 | O43474_416_419 | O75689_138_141 | O95251_171_174 | P03431_669_672 | P0DMU9_29_32 | P0DMU9_17_20 | P20719_196_199 | P23025_31_34 | P49759-3_178_181 | P61289_243_246 | P78337_146_149 | P86480_8_11 | Q00444_156_159 | Q14498_72_75 | Q14498_56_59 | Q15149_1951_1954 | Q15149_4262_4265 | Q3V6T2_839_842 | Q5S007_947_950 | Q5SQQ9-2_101_104 | Q66GS9_289_292 | Q67296_752_755 | Q67296_736_739 | Q6ZU80_141_144 | Q7Z7J5_155_158 | Q86YD7_157_160 | Q8IX06_226_229 | Q8N7W2-2_294_297 | Q8N8U2_105_108 | Q8N960_689_692 | Q8NE31_526_529 | Q8NHU0_17_20 | Q8NHU0_29_32 | Q8TAD8_257_260 | Q8TAD8_95_98 | Q8TBE0_169_172 | Q96C24_34_37 | Q96EB6_35_38 | Q96EG3_232_235 | Q96KM6_843_846 | Q96MU7_302_305 | Q96PV6-2_489_492 | Q96T60_138_141 | Q9BTC0_1214_1217 | Q9H0I2_113_116 | Q9HAJ7_87_90 | Q9HBE1-4_347_350 | Q9HBE1-4_264_267 | Q9P127_230_233 | Q9UBU8-2_136_139 | Q9UFB7_648_651 | Q9UK58_384_387)
## 2: query: P03428, Pattern:KR.R, interacts_with: (Q92997), human_instances: (when query interacts_with: Q92997 supporting human_instances: O00339_125_128 | O15169_414_417 | O15481_14_17 | O43159_131_134 | O43167_159_162 | O43167_161_164 | O43463_394_397 | O43474_416_419 | O75689_138_141 | O95251_171_174 | P03431_669_672 | P0DMU9_29_32 | P0DMU9_17_20 | P20719_196_199 | P23025_31_34 | P49759-3_178_181 | P61289_243_246 | P78337_146_149 | P86480_8_11 | Q00444_156_159 | Q14498_72_75 | Q14498_56_59 | Q15149_1951_1954 | Q15149_4262_4265 | Q3V6T2_839_842 | Q5S007_947_950 | Q5SQQ9-2_101_104 | Q66GS9_289_292 | Q67296_752_755 | Q67296_736_739 | Q6ZU80_141_144 | Q7Z7J5_155_158 | Q86YD7_157_160 | Q8IX06_226_229 | Q8N7W2-2_294_297 | Q8N8U2_105_108 | Q8N960_689_692 | Q8NE31_526_529 | Q8NHU0_17_20 | Q8NHU0_29_32 | Q8TAD8_257_260 | Q8TAD8_95_98 | Q8TBE0_169_172 | Q96C24_34_37 | Q96EB6_35_38 | Q96EG3_232_235 | Q96KM6_843_846 | Q96MU7_302_305 | Q96PV6-2_489_492 | Q96T60_138_141 | Q9BTC0_1214_1217 | Q9H0I2_113_116 | Q9HAJ7_87_90 | Q9HBE1-4_347_350 | Q9HBE1-4_264_267 | Q9P127_230_233 | Q9UBU8-2_136_139 | Q9UFB7_648_651 | Q9UK58_384_387)
## 3: query: P03428, Pattern:KR.R, interacts_with: (Q92997), human_instances: (when query interacts_with: Q92997 supporting human_instances: O00339_125_128 | O15169_414_417 | O15481_14_17 | O43159_131_134 | O43167_159_162 | O43167_161_164 | O43463_394_397 | O43474_416_419 | O75689_138_141 | O95251_171_174 | P03431_669_672 | P0DMU9_29_32 | P0DMU9_17_20 | P20719_196_199 | P23025_31_34 | P49759-3_178_181 | P61289_243_246 | P78337_146_149 | P86480_8_11 | Q00444_156_159 | Q14498_72_75 | Q14498_56_59 | Q15149_1951_1954 | Q15149_4262_4265 | Q3V6T2_839_842 | Q5S007_947_950 | Q5SQQ9-2_101_104 | Q66GS9_289_292 | Q67296_752_755 | Q67296_736_739 | Q6ZU80_141_144 | Q7Z7J5_155_158 | Q86YD7_157_160 | Q8IX06_226_229 | Q8N7W2-2_294_297 | Q8N8U2_105_108 | Q8N960_689_692 | Q8NE31_526_529 | Q8NHU0_17_20 | Q8NHU0_29_32 | Q8TAD8_257_260 | Q8TAD8_95_98 | Q8TBE0_169_172 | Q96C24_34_37 | Q96EB6_35_38 | Q96EG3_232_235 | Q96KM6_843_846 | Q96MU7_302_305 | Q96PV6-2_489_492 | Q96T60_138_141 | Q9BTC0_1214_1217 | Q9H0I2_113_116 | Q9HAJ7_87_90 | Q9HBE1-4_347_350 | Q9HBE1-4_264_267 | Q9P127_230_233 | Q9UBU8-2_136_139 | Q9UFB7_648_651 | Q9UK58_384_387)
## 4: query: P03428, Pattern:KR.R, interacts_with: (Q92997), human_instances: (when query interacts_with: Q92997 supporting human_instances: O00339_125_128 | O15169_414_417 | O15481_14_17 | O43159_131_134 | O43167_159_162 | O43167_161_164 | O43463_394_397 | O43474_416_419 | O75689_138_141 | O95251_171_174 | P03431_669_672 | P0DMU9_29_32 | P0DMU9_17_20 | P20719_196_199 | P23025_31_34 | P49759-3_178_181 | P61289_243_246 | P78337_146_149 | P86480_8_11 | Q00444_156_159 | Q14498_72_75 | Q14498_56_59 | Q15149_1951_1954 | Q15149_4262_4265 | Q3V6T2_839_842 | Q5S007_947_950 | Q5SQQ9-2_101_104 | Q66GS9_289_292 | Q67296_752_755 | Q67296_736_739 | Q6ZU80_141_144 | Q7Z7J5_155_158 | Q86YD7_157_160 | Q8IX06_226_229 | Q8N7W2-2_294_297 | Q8N8U2_105_108 | Q8N960_689_692 | Q8NE31_526_529 | Q8NHU0_17_20 | Q8NHU0_29_32 | Q8TAD8_257_260 | Q8TAD8_95_98 | Q8TBE0_169_172 | Q96C24_34_37 | Q96EB6_35_38 | Q96EG3_232_235 | Q96KM6_843_846 | Q96MU7_302_305 | Q96PV6-2_489_492 | Q96T60_138_141 | Q9BTC0_1214_1217 | Q9H0I2_113_116 | Q9HAJ7_87_90 | Q9HBE1-4_347_350 | Q9HBE1-4_264_267 | Q9P127_230_233 | Q9UBU8-2_136_139 | Q9UFB7_648_651 | Q9UK58_384_387)
## 5: query: P03428, Pattern:KR.R, interacts_with: (Q92997), human_instances: (when query interacts_with: Q92997 supporting human_instances: O00339_125_128 | O15169_414_417 | O15481_14_17 | O43159_131_134 | O43167_159_162 | O43167_161_164 | O43463_394_397 | O43474_416_419 | O75689_138_141 | O95251_171_174 | P03431_669_672 | P0DMU9_29_32 | P0DMU9_17_20 | P20719_196_199 | P23025_31_34 | P49759-3_178_181 | P61289_243_246 | P78337_146_149 | P86480_8_11 | Q00444_156_159 | Q14498_72_75 | Q14498_56_59 | Q15149_1951_1954 | Q15149_4262_4265 | Q3V6T2_839_842 | Q5S007_947_950 | Q5SQQ9-2_101_104 | Q66GS9_289_292 | Q67296_752_755 | Q67296_736_739 | Q6ZU80_141_144 | Q7Z7J5_155_158 | Q86YD7_157_160 | Q8IX06_226_229 | Q8N7W2-2_294_297 | Q8N8U2_105_108 | Q8N960_689_692 | Q8NE31_526_529 | Q8NHU0_17_20 | Q8NHU0_29_32 | Q8TAD8_257_260 | Q8TAD8_95_98 | Q8TBE0_169_172 | Q96C24_34_37 | Q96EB6_35_38 | Q96EG3_232_235 | Q96KM6_843_846 | Q96MU7_302_305 | Q96PV6-2_489_492 | Q96T60_138_141 | Q9BTC0_1214_1217 | Q9H0I2_113_116 | Q9HAJ7_87_90 | Q9HBE1-4_347_350 | Q9HBE1-4_264_267 | Q9P127_230_233 | Q9UBU8-2_136_139 | Q9UFB7_648_651 | Q9UK58_384_387)
## 6: query: P03428, Pattern:KR.R, interacts_with: (Q92997), human_instances: (when query interacts_with: Q92997 supporting human_instances: O00339_125_128 | O15169_414_417 | O15481_14_17 | O43159_131_134 | O43167_159_162 | O43167_161_164 | O43463_394_397 | O43474_416_419 | O75689_138_141 | O95251_171_174 | P03431_669_672 | P0DMU9_29_32 | P0DMU9_17_20 | P20719_196_199 | P23025_31_34 | P49759-3_178_181 | P61289_243_246 | P78337_146_149 | P86480_8_11 | Q00444_156_159 | Q14498_72_75 | Q14498_56_59 | Q15149_1951_1954 | Q15149_4262_4265 | Q3V6T2_839_842 | Q5S007_947_950 | Q5SQQ9-2_101_104 | Q66GS9_289_292 | Q67296_752_755 | Q67296_736_739 | Q6ZU80_141_144 | Q7Z7J5_155_158 | Q86YD7_157_160 | Q8IX06_226_229 | Q8N7W2-2_294_297 | Q8N8U2_105_108 | Q8N960_689_692 | Q8NE31_526_529 | Q8NHU0_17_20 | Q8NHU0_29_32 | Q8TAD8_257_260 | Q8TAD8_95_98 | Q8TBE0_169_172 | Q96C24_34_37 | Q96EB6_35_38 | Q96EG3_232_235 | Q96KM6_843_846 | Q96MU7_302_305 | Q96PV6-2_489_492 | Q96T60_138_141 | Q9BTC0_1214_1217 | Q9H0I2_113_116 | Q9HAJ7_87_90 | Q9HBE1-4_347_350 | Q9HBE1-4_264_267 | Q9P127_230_233 | Q9UBU8-2_136_139 | Q9UFB7_648_651 | Q9UK58_384_387)
## 7: query: P03428, Pattern:KR.R, interacts_with: (Q92997), human_instances: (when query interacts_with: Q92997 supporting human_instances: O00339_125_128 | O15169_414_417 | O15481_14_17 | O43159_131_134 | O43167_159_162 | O43167_161_164 | O43463_394_397 | O43474_416_419 | O75689_138_141 | O95251_171_174 | P03431_669_672 | P0DMU9_29_32 | P0DMU9_17_20 | P20719_196_199 | P23025_31_34 | P49759-3_178_181 | P61289_243_246 | P78337_146_149 | P86480_8_11 | Q00444_156_159 | Q14498_72_75 | Q14498_56_59 | Q15149_1951_1954 | Q15149_4262_4265 | Q3V6T2_839_842 | Q5S007_947_950 | Q5SQQ9-2_101_104 | Q66GS9_289_292 | Q67296_752_755 | Q67296_736_739 | Q6ZU80_141_144 | Q7Z7J5_155_158 | Q86YD7_157_160 | Q8IX06_226_229 | Q8N7W2-2_294_297 | Q8N8U2_105_108 | Q8N960_689_692 | Q8NE31_526_529 | Q8NHU0_17_20 | Q8NHU0_29_32 | Q8TAD8_257_260 | Q8TAD8_95_98 | Q8TBE0_169_172 | Q96C24_34_37 | Q96EB6_35_38 | Q96EG3_232_235 | Q96KM6_843_846 | Q96MU7_302_305 | Q96PV6-2_489_492 | Q96T60_138_141 | Q9BTC0_1214_1217 | Q9H0I2_113_116 | Q9HAJ7_87_90 | Q9HBE1-4_347_350 | Q9HBE1-4_264_267 | Q9P127_230_233 | Q9UBU8-2_136_139 | Q9UFB7_648_651 | Q9UK58_384_387)
```

```r
IntAct_net[unique_position == "P06463_148_155"]
```

```
##     Pattern Viral_protein Human_protein unique_position Motif_Sig
## 1: E..QR..E        P06463        Q99460  P06463_148_155  0.000853
## 2: E..QR..E        P06463        Q99460  P06463_148_155  0.000853
##    Domain_p.value top_domain_support4motif_nq domain_support4motif_nq
## 1:     0.01142691                          NA                      NA
## 2:     0.01142691                          NA                      NA
##    domain_support4motif_q Domain_combined_p.value Human_recognition_domain
## 1:                     NA                      NA                IPR011989
## 2:                     NA                      NA                IPR016024
##    matching_known
## 1: matching_known
## 2: matching_known
##                                                                                                                                                                                                                                    human_instances
## 1: query: P06463, Pattern:E..QR..E, interacts_with: (Q99460), human_instances: (when query interacts_with: Q99460 supporting human_instances: D3DR86_168_175 | P19532_410_417 | P55036_226_233 | Q99459_772_779 | Q9UBN6_312_319 | Q9UKE5_374_381)
## 2: query: P06463, Pattern:E..QR..E, interacts_with: (Q99460), human_instances: (when query interacts_with: Q99460 supporting human_instances: D3DR86_168_175 | P19532_410_417 | P55036_226_233 | Q99459_772_779 | Q9UBN6_312_319 | Q9UKE5_374_381)
```

### find ELM domains and plot p-values for those


```r
# read domain annotations
InterProScan_domains = readInterProGFF3("./processed_data_files/all_human_viral_protein_domains.gff3.gz")
# get InterProID to member database ID mapping
InterPro2memberDB = getInterPro2memberDB(InterProScan_domains)
InterPro2memberDB = InterPro2memberDB[complete.cases(InterPro2memberDB)]

# list of domains
if(!file.exists("./data_files/interactiondomains.tsv")) downloader::download("http://elm.eu.org/interactiondomains.tsv","./data_files/interactiondomains.tsv")
interactiondomains = fread("./data_files/interactiondomains.tsv")
interactiondomains[, pfam_id := `Interaction Domain Id`]

# list of interactions
if(!file.exists("./data_files/elm_interactions.tsv")) downloader::download("http://elm.eu.org/interactions/as_tsv","./data_files/elm_interactions.tsv")
elm_interactions = fread("./data_files/elm_interactions.tsv")
elm_interactions[, memberDBID := Domain]
elm_interactions = merge(elm_interactions, InterPro2memberDB,
                         by.x = "memberDBID", by.y = "memberDBID",
                         all.x = T, all.y = F)
elm_interactions[memberDBID %in% InterPro2memberDB$InterProID, InterProID := memberDBID]
# pick interactions where domain is in a human protein
elm_interactions = elm_interactions[grepl(9606,taxonomyDomain)]

domains_known = interactiondomains[, unique(pfam_id)]

InterProScan_domains = readInterProGFF3("./processed_data_files/all_human_viral_protein_domains.gff3.gz")
# get InterProID to member database ID mapping
InterPro2memberDB = getInterPro2memberDB(InterProScan_domains)
InterPro2memberDB = InterPro2memberDB[complete.cases(InterPro2memberDB)]
domains_known_mapped = unique(InterPro2memberDB[memberDBID %in% domains_known | InterProID %in% domains_known, InterProID])

dom_viral = R.utils::env(load("./processed_data_files/domain_res_count_20171019.RData"))
dom_viral = dom_viral$res_count$data_with_pval
dom_viral[, network := "predicted to bind\nviral proteins"]
dom_human_env = R.utils::env(load("./processed_data_files/predict_domain_human_clust20180819.RData"))
dom_human = dom_human_env$res_count_all$data_with_pval
dom_human[, network := "predicted to bind\nhuman proteins"]
dom_human_BioPlex = dom_human_env$res_count_BioPlex3$data_with_pval
dom_human_BioPlex[, network := "predicted to bind\nhuman proteins (BioPlex)"]
rm(dom_human_env)
##### / standardise column names
std_cols = c("protein_with_motif", "protein_with_domain", "domain", "Taxid_protein_with_motif","Taxid_protein_with_domain")
domain_res_cols = c("IDs_interactor_viral", "IDs_interactor_human", "IDs_domain_human", "Taxid_interactor_human","Taxid_interactor_viral")
for (col_name in domain_res_cols) {
    setnames(dom_viral, colnames(dom_viral),
             gsub(col_name,
                  std_cols[which(domain_res_cols == col_name)],
                  colnames(dom_viral)))
}
non_query_domain_res_cols = c("IDs_interactor_human_A", "IDs_interactor_human_B", "IDs_domain_human_B", "Taxid_interactor_human_A","Taxid_interactor_human_B")
for (col_name_nq in non_query_domain_res_cols) {
    setnames(dom_human, colnames(dom_human),
             gsub(col_name_nq,
                  std_cols[which(non_query_domain_res_cols == col_name_nq)],
                  colnames(dom_human)))
    setnames(dom_human_BioPlex, colnames(dom_human_BioPlex),
             gsub(col_name_nq,
                  std_cols[which(non_query_domain_res_cols == col_name_nq)],
                  colnames(dom_human_BioPlex)))
}
##### \
dom = rbind(dom_viral, dom_human, dom_human_BioPlex) # dom_human_BioPlex
dom[domain %in% domains_known_mapped, SLIM_binding := "yes"]
dom[!domain %in% domains_known_mapped, SLIM_binding := "no / not known"]
# how many human proteins have SLIM-binding domains
dom[SLIM_binding == "yes", uniqueN(protein_with_domain)]
```

```
## [1] 3837
```

```r
# how many viral proteins target human proteins that have SLIM-binding domains
dom[SLIM_binding == "yes", uniqueN(protein_with_motif)]
```

```
## [1] 17420
```

```r
dom2 = unique(dom[,.(protein_with_motif, domain, p.value, domain_type, SLIM_binding, network)])
# how many duplicated p-values?
dom2[,.N, by =.(protein_with_motif, domain, domain_type, SLIM_binding, network)][,table(N)]
```

```
## N
##      1 
## 863015
```

```r
# add known elm domains to interactions
dom_elm = merge(dom, elm_interactions,
            by.x = c("protein_with_motif", "protein_with_domain", "domain"),
            by.y = c("interactorElm", "interactorDomain", "InterProID"),
            all.x = T, all.y = F)
dom_elm[!is.na(Elm), Correct_SLiM_binding := "yes"]
dom_elm[is.na(Elm), Correct_SLiM_binding := "no"]

###### \ plot slim-binding domains
ggplot(dom2, aes(p.value, color = SLIM_binding, fill = SLIM_binding)) +
    geom_histogram(alpha = 0.3, aes(y = ..density..), position = "identity") +
    theme_bw() +
    facet_grid(network~., scales = "free_y")+
    theme(legend.position = "right",
          strip.text.y = element_text(angle = 0),
          strip.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          text = element_text(size = 14))+
  scale_color_manual(values=nice_colors)+
  scale_fill_manual(values=nice_colors)
```

![](/nfs/research1/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/ELM_domains-1.png)<!-- -->

```r
###### /
ks.test(dom2[network == "predicted to bind\nviral proteins" & 
                SLIM_binding == "no / not known", p.value],
        dom2[network == "predicted to bind\nviral proteins" &
                SLIM_binding == "yes", p.value], 
        alternative = "less")
```

```
## 
## 	Two-sample Kolmogorov-Smirnov test
## 
## data:  dom2[network == "predicted to bind\nviral proteins" & SLIM_binding ==  and dom2[network == "predicted to bind\nviral proteins" & SLIM_binding ==     "no / not known", p.value] and     "yes", p.value]
## D^- = 0.13548, p-value < 2.2e-16
## alternative hypothesis: the CDF of x lies below that of y
```

```r
ks.test(dom2[network ==  "predicted to bind\nhuman proteins" & 
                SLIM_binding == "no / not known", p.value],
        dom2[network == "predicted to bind\nhuman proteins" &
                SLIM_binding == "yes", p.value], 
        alternative = "less")
```

```
## 
## 	Two-sample Kolmogorov-Smirnov test
## 
## data:  dom2[network == "predicted to bind\nhuman proteins" & SLIM_binding ==  and dom2[network == "predicted to bind\nhuman proteins" & SLIM_binding ==     "no / not known", p.value] and     "yes", p.value]
## D^- = 0.13704, p-value < 2.2e-16
## alternative hypothesis: the CDF of x lies below that of y
```

```r
dom05 = dom2[p.value < 0.5]
dom05[, uniqueN(protein_with_motif)]
```

```
## [1] 13380
```

```r
dom05[, uniqueN(domain)]
```

```
## [1] 2407
```

```r
# how many motifs per protein
table(IntAct$res$qslimfinder.Full_IntAct3.FALSE$instances_query$ID, IntAct$res$qslimfinder.Full_IntAct3.FALSE$instances_query$Primary_Acc)[, table(IntAct$res$qslimfinder.Full_IntAct3.FALSE$instances_query$Primary_Acc) > 2]
```

```
##                     
##                      P03255 P27958 P69616
##   DOC_USP7_MATH_2         0      0      0
##   DOC_WW_Pin1_4           0      0      0
##   LIG_G3BP_FGDF_1         0      0      0
##   LIG_LIR_Gen_1           0      0      0
##   LIG_MYND_1              1      0      0
##   LIG_PDZ_Class_1         0      0      0
##   LIG_PTAP_UEV_1          0      0      1
##   LIG_Rb_LxCxE_1          1      0      0
##   LIG_Rb_pABgroove_1      1      0      0
##   LIG_SH3_3               0      0      2
##   LIG_TRAF2_2             0      0      0
##   LIG_WW_1                0      0      0
##   MOD_N-GLC_1             0      4      0
##   MOD_NMyristoyl          0      0      0
```
# Calculate and plot rank (by p-value) of SLIM-binding domains for each protein

This is more interpretable because a high p-value of 0.2 can still mean top-1 rank for proteins with low number of interactors. -Not very useful :(  


```r
setorder(dom2, protein_with_motif, p.value)
dom2[, rank := as.integer(frank(p.value, ties.method = "min", na.last = NA)), by = .(protein_with_motif)]
ggplot(dom2, aes(rank, p.value)) + # color = SLIM_binding, fill = SLIM_binding
    geom_bin2d(bins = 30, aes(colour = ..count.., fill = ..count..)) + theme_bw() +
    facet_grid(network~SLIM_binding, scales = "free_y") +
    ylim(0, 0.3)+
    theme(legend.position = "right",
          strip.text.y = element_text(angle = 0),
          strip.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          text = element_text(size = 14))
```

![](/nfs/research1/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/slim_domain_ranks-1.png)<!-- -->

# Calculate and plot p-values for known SLIM-mediated interactions

This is a more direct benchmark.


```r
# filter and select unique to have one row per protein_with_motif - domain pair
dom3 = unique(dom_elm[,.(protein_with_motif, domain, p.value, domain_type, Correct_SLiM_binding, network)])
dom3[, both_SLIM_interaction := uniqueN(Correct_SLiM_binding) >= 2,
     by = .(protein_with_motif, network)]
dom3 = dom3[both_SLIM_interaction == T]
# domain protein pairs in each network
dom3[, N_pairs := uniqueN(paste0(protein_with_motif, domain)), by = .(network, Correct_SLiM_binding)]
dom3[, network := paste0(network,"\n",
                         unique(N_pairs[Correct_SLiM_binding == "yes"]), 
                         " correct domains ",
                         "/ ", unique(N_pairs[Correct_SLiM_binding == "no"]),
                         " total"), by = .(network)]
ggplot(dom3, aes(p.value, color = Correct_SLiM_binding, fill = Correct_SLiM_binding)) +
    geom_histogram(alpha = 0.3, aes(y = ..density..), bins = 20, position = "identity") +
    theme_bw() +
    facet_grid(network~., scales = "free_y")+
    theme(legend.position = "none",
          strip.text.y = element_text(angle = 0),
          strip.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 14),
          text = element_text(size = 14))+
  scale_color_manual(values=nice_colors)+
  scale_fill_manual(values=nice_colors)
```

![](/nfs/research1/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/slim_interactions-1.png)<!-- -->

## Calculate precision - recall curves for known SLIM-mediated interactions - in progress
.    known      predicted  
TP:  PDZ        PDZ  
FN:  PDZ        no PDZ   
FP:  PDZ        ATP-ase, SH3 (2)  
TN:  not SH2    

First, look at the confusion matrix - in progress


```r
conf = merge(dom, elm_interactions,
            by.x = c("protein_with_motif", "protein_with_domain"),
            by.y = c("interactorElm", "interactorDomain"),
            all.x = T, all.y = F)
conf = unique(conf[,.(protein_with_motif, domain, InterProID, p.value)])
conf = conf[p.value < 0.5]
conf = dcast.data.table(conf, domain ~ InterProID,
                        value.var = "protein_with_motif",
                        fun.aggregate = uniqueN)
rownames = conf$domain
conf[, domain := NULL]
conf = as.matrix(conf)
rownames(conf) = rownames
conf = conf >= 1
conf = conf[colnames(conf),]
image(conf)
```


```r
# get readable names for proteins and domains 
readableNet = function(Net){
    proteins2gene_names = lapply(unique(c(Net$Viral_protein, Net$Human_protein)), function(uniprot_id) {
        suppressMessages({
            proteins2gene_names = fread(paste0("https://www.uniprot.org/uniprot/?query=id:", uniprot_id,
                                               "&format=tab&columns=id,genes,entry%20name,protein%20names,organism,organism-id"),
                                        stringsAsFactors = F)
        })
        proteins2gene_names
    })
    proteins2gene_names = Reduce(rbind, proteins2gene_names)
    setnames(proteins2gene_names, colnames(proteins2gene_names),
             c( "Entry", "Gene_names", "Entry_name", "Protein_names", "Organism", "Organism_ID"))
    proteins2gene_names[, Gene_names := gsub(" .*", "", Gene_names)] # pick only first gene name
    proteins2gene_names[, Protein_names := gsub(" \\(.*", "", Protein_names)] # pick only first protein name
    # fwrite(proteins2gene_names, "./results/IntAct_network_proteins2gene_names", sep = "\t")
    proteins2gene_names_human = proteins2gene_names[Organism_ID == 9606]
    proteins2gene_names_viral = proteins2gene_names[Organism_ID != 9606]
    
    proteins2gene_names_human[, c("Organism", "Organism_ID") := NULL]
    proteins2gene_names_viral[, c("Organism", "Organism_ID") := NULL]
    
    setnames(proteins2gene_names_human, colnames(proteins2gene_names_human),
             c("Human_protein", "Human_gene_names", "Human_entry_name", "Human_protein_names"))
    setnames(proteins2gene_names_viral, colnames(proteins2gene_names_viral),
             c("Viral_protein", "Viral_gene_names", "Viral_entry_name", "Viral_protein_names"))
    
    Net = merge(Net, proteins2gene_names_human, by = "Human_protein", all.x = T, all.y = F)
    Net = merge(Net, proteins2gene_names_viral, by = "Viral_protein", all.x = T, all.y = F)
    
    InterProEntryTypes = getInterProEntryTypes("./data_files/entry.list")
    InterProEntryTypes[, ENTRY_TYPE := NULL]
    setnames(InterProEntryTypes, colnames(InterProEntryTypes),
             c("Human_recognition_domain", "Human_recognition_domain_name"))
    Net = merge(Net, InterProEntryTypes, by = "Human_recognition_domain", all.x = T, all.y = F)
    
    # mark known SLIM-interacting domains
    Net[, in_ELM := Human_recognition_domain %in% domains_known_mapped]
    
    # check if some columns are in Net
    if(!"domain_support4motif_nq" %in% colnames(Net)) Net[, domain_support4motif_nq := NA]
    if(!"domain_support4motif_q" %in% colnames(Net)) Net[, domain_support4motif_q := NA]
    if(!"top_domain_support4motif_nq" %in% colnames(Net)) Net[, top_domain_support4motif_nq := NA]
    
    # reorder column names 
    Net = Net[, .(Pattern, Human_recognition_domain_name, Domain_p.value, 
                  top_domain_support4motif_nq, domain_support4motif_nq,
                  domain_support4motif_q, unique_position, 
                  Human_protein, in_ELM, matching_known, human_instances,
                  Viral_protein_names, Viral_entry_name, Human_protein_names, Human_entry_name,
                  Human_recognition_domain, Motif_Sig, Viral_protein,
                  Viral_gene_names, Human_gene_names)]
    setorder(Net, Viral_protein, Motif_Sig, Domain_p.value)
    unique(Net)
}

IntAct_net = readableNet(IntAct_net) 
IntAct_filt_hum_dom_net = readableNet(IntAct_filt_hum_dom_net) 
viral_net = readableNet(viral_net)

fwrite(IntAct_net, "./results/IntAct_net_motif_domain_readable.tsv", sep = "\t")
fwrite(IntAct_filt_hum_dom_net, "./results/IntAct_filt_hum_dom_net_motif_domain_readable.tsv", sep = "\t")
fwrite(viral_net, "./results/viral_net_motif_domain_readable.tsv", sep = "\t")

fwrite(unique(IntAct_net[,.(Pattern)]), "./results/IntAct_net_patterns.tsv", sep = "\t")
fwrite(unique(IntAct_filt_hum_dom_net[,.(Pattern)]), "./results/IntAct_filt_hum_dom_net_patterns.tsv", sep = "\t")
fwrite(unique(viral_net[,.(Pattern)]), "./results/viral_net_patterns.tsv", sep = "\t")
```

Armadillo like domains:  
Viral proteins 11.  
Viral motifs 14.  
Human proteins 4.   

PDZ domains:  
Viral proteins 14.  
Viral motifs 24.  
Human proteins 7.   

PDZ motifs:  
Viral proteins 16.  
Viral motifs 26.  
Human proteins 8.   



```r
motifSummary = function(benchmarkMotifsResult = IntAct$res$qslimfinder.Full_IntAct3.FALSE,
                        readable_net = IntAct_net,
                        motif = "E.V..G.{0,2}N.{0,1}Q",
                        viral_query = NULL,
                        human_seed = NULL,
                        results_dir = "./results/"){
    if(is.null(viral_query)) viral_query = readable_net$Viral_protein
    if(is.null(human_seed)) human_seed = readable_net$Human_protein
    if(is.null(motif)) motif = readable_net$Pattern
    
    readable_net_ind = readable_net$Pattern %in% motif & 
        readable_net$Viral_protein %in% viral_query & 
        readable_net$Human_protein %in% human_seed
    occurence_query_ind = benchmarkMotifsResult$occurence_query$Pattern %in% motif &
        benchmarkMotifsResult$occurence_query$query %in% viral_query &
        benchmarkMotifsResult$occurence_query$interacts_with %in% human_seed
    
    occurence_query = benchmarkMotifsResult$occurence_query[occurence_query_ind]
    occurence_query = data.table(UP = occurence_query$UP, UPNum = occurence_query$UPNum,
                                 unique_position = paste0(GenomicRanges::seqnames(occurence_query), "_",
                                                     GenomicRanges::start(occurence_query), "_",
                                                     GenomicRanges::end(occurence_query)),
                                 Occ = occurence_query$Occ, SeqNum = occurence_query$SeqNum, 
                                 domain_count_per_protein_with_motif = occurence_query$domain_count_per_protein_with_motif,
                                 top_domain_support4motif_nq = occurence_query$top_domain_support4motif_nq,
                                 domain_support4motif_nq = occurence_query$domain_support4motif_nq,
                                 protein_with_motif_degree = occurence_query$protein_with_motif_degree,
                                 domain = occurence_query$domain,
                                 Human_protein = occurence_query$interacts_with)
    occurence_query = unique(occurence_query)
    
    occurence_ind = benchmarkMotifsResult$occurence$Pattern %in% motif &
        benchmarkMotifsResult$occurence$query %in% viral_query &
        benchmarkMotifsResult$occurence$interacts_with %in% human_seed
    occurence_seqnames = as.character(GenomicRanges::seqnames(benchmarkMotifsResult$occurence[occurence_ind]))
    occurence_seqnames = unique(occurence_seqnames)
    
    summary = paste0(
        "viral query protein with motif: ",
        paste0(unique(paste0(readable_net[readable_net_ind]$unique_position, " = ",
                             readable_net[readable_net_ind]$Pattern, " = ",
                             readable_net[readable_net_ind]$Viral_protein_names, " (",
                             readable_net[readable_net_ind]$Viral_entry_name, ") ")), collapse = " | "), "\n",
        uniqueN(paste0(readable_net[readable_net_ind]$unique_position, "=", readable_net[readable_net_ind]$Viral_protein_names)),
        "\n\n",
        "viral targeted human proteins: ",
        paste0(unique(paste0(readable_net[readable_net_ind]$Human_protein, " = ", readable_net[readable_net_ind]$Human_protein_names)), collapse = " | "), "\n",
        uniqueN(paste0(readable_net[readable_net_ind]$Human_protein, "=", readable_net[readable_net_ind]$Human_protein_names)),
        "\n\n",
        "UP support / UP: ",
        paste0(unique(paste0("(", occurence_query$unique_position," in ", occurence_query$Human_protein," dataset) ",
            occurence_query$UP, "/",
            occurence_query$UPNum)), collapse = " | "),
        "\n\n",
        "Occurence / SeqNum: ",
        paste0(unique(paste0("(", occurence_query$unique_position," in ", occurence_query$Human_protein," dataset) ",
            occurence_query$Occ, "/",
            occurence_query$SeqNum)), collapse = " | "),
        "\n\n",
        "top-1 domain(s): ",
        paste0(paste0("(", occurence_query$unique_position,") ", occurence_query$domain), collapse = " | "),
        "\n\n",
        "domain_support4motif_nq / protein_with_motif_degree: ",
        paste0(paste0(
            occurence_query$domain_support4motif_nq, "/",
            occurence_query$protein_with_motif_degree), collapse = " | "),
        "\n\n",
        "all with motif: ",
        paste0(occurence_seqnames, collapse = " | "), "\n",
        "N ", uniqueN(occurence_seqnames),
        "\nN non-query with motif: ", uniqueN(occurence_seqnames) -
            uniqueN(readable_net[readable_net_ind]$Viral_protein),
        "\n\n"
    )
    
    readable_net2 = readable_net[readable_net_ind]
    # For Cytoscape: network table
    Nrows = nrow(readable_net2)
    readable_net3 = data.table(source = c(readable_net2$Human_protein,   readable_net2$Human_recognition_domain, readable_net2$Pattern),
                               target = c(readable_net2$Human_recognition_domain, readable_net2$Pattern,          readable_net2$Viral_protein),
                               source_readable = c(paste0(readable_net2$Human_protein_names, " (", readable_net2$Human_gene_names, ")"),  
                                                   readable_net2$Human_recognition_domain_name,
                                                   readable_net2$Pattern),
                               target_readable = c(readable_net2$Human_recognition_domain_name,
                                                   readable_net2$Pattern,
                                                   paste0(readable_net2$Viral_protein_names, " (", readable_net2$Viral_entry_name, ")")),
                               human_instances = c(readable_net2$human_instances, readable_net2$human_instances, readable_net2$human_instances),
                               top_domain_support4motif_nq =    c(rep(1, Nrows),  readable_net2$top_domain_support4motif_nq,   rep(1, Nrows)),
                               domain_support4motif_nq =    c(rep(1, Nrows),  readable_net2$domain_support4motif_nq,   rep(1, Nrows)),
                               type_source = c(rep("human_protein", Nrows), rep("human_domain", Nrows),  rep("viral_motif", Nrows)),
                               type_target = c(rep("human_domain", Nrows), rep("viral_motif", Nrows),  rep("viral_protein", Nrows)))
    # For Cytoscape: node table
    rescaled_motif_pval = readable_net2[, (1 - Motif_Sig/ max(readable_net$Motif_Sig)) *
                                            max(readable_net$top_domain_support4motif_nq)]
    rescaled_motif_pval[is.nan(rescaled_motif_pval)]
    network_types = unique(data.table(node = c(readable_net3$source, readable_net3$target),
                                      node_readable = c(readable_net3$source_readable, readable_net3$target_readable),
                                      type = c(readable_net3$type_source, readable_net3$type_target),
                                      p_value = c(rep("", Nrows), readable_net2$Domain_p.value, readable_net2$Motif_Sig,
                                                     readable_net2$Domain_p.value, readable_net2$Motif_Sig, rep("", Nrows)),
                                      N_human_instances_bind_domain = c(rep("", Nrows), readable_net2[, paste0("top: ", top_domain_support4motif_nq,", all: ",domain_support4motif_nq)], rep("", Nrows),
                                                     readable_net2[, paste0("top: ", top_domain_support4motif_nq,", all: ",domain_support4motif_nq)], rep("", Nrows), rep("", Nrows)),
                                      scale_by = c(rep("", Nrows), readable_net2$top_domain_support4motif_nq, rescaled_motif_pval,
                                                     readable_net2$top_domain_support4motif_nq, rescaled_motif_pval, rep("", Nrows)),
                                      scale_by2 = c(rep("", Nrows), readable_net2$domain_support4motif_nq, rescaled_motif_pval,
                                                     readable_net2$domain_support4motif_nq, rescaled_motif_pval, rep("", Nrows))))
    readable_net3_file = paste0(as.character(as.list(match.call())[["readable_net"]]), "_",
                                paste0(unique(readable_net[readable_net_ind]$unique_position), collapse = "_"), "_", 
                                paste0(unique(readable_net[readable_net_ind]$Human_protein), collapse = "_"))
    if(nchar(readable_net3_file) > 100) readable_net3_file = digest::digest(readable_net3_file, "sha1")
    if(!dir.exists(results_dir)) dir.create(results_dir, recursive = T)
    readable_net3_file = paste0(results_dir, readable_net3_file)
    readable_net3 = unique(readable_net3[complete.cases(readable_net3)])
    fwrite(readable_net3, 
           paste0(readable_net3_file, "_network"),
           sep = "\t")
    fwrite(network_types,
           paste0(readable_net3_file,"_node_table"),
           sep = "\t")
    
    cat(summary)
    list(network = readable_net3, node_table = network_types)
}
```

# files for cytoscape networks (for thesis and for paper)


```r
# use 
bench_res_cyt = IntAct_filt_hum_dom$res$qslimfinder.Full_IntAct3.FALSE
IntAct_net = bechmark2Net(bench_res_cyt,
                          top_domain_by = c(NA))
IntAct_net = readableNet(IntAct_net)

# all motifs and domains
all = motifSummary(bench_res_cyt, IntAct_net,
             motif = NULL, viral_query = NULL, human_seed = NULL,
             results_dir = "./results_w_human_dom/")
```

```
## viral query protein with motif: B4URF7_736_739 = KRKR = Polymerase basic protein 2 (B4URF7_I33A0)  | B4URF7_188_194 = E..E..Q = Polymerase basic protein 2 (B4URF7_I33A0)  | B4URF7_69_75 = E..E..Q = Polymerase basic protein 2 (B4URF7_I33A0)  | B4URF7_737_739 = RKR = Polymerase basic protein 2 (B4URF7_I33A0)  | C5E522_20_25 = Q..T.I = Nucleoprotein (C5E522_9INFA)  | C5E522_443_445 = EVI = Nucleoprotein (C5E522_9INFA)  | C5E526_205_214 = I..K.{1,2}R.{1,2}K = RNA-directed RNA polymerase catalytic subunit (C5E526_9INFA)  | C5E527_158_163 = E..D.I = Polymerase basic protein 2 (C5E527_9INFA)  | C5E527_100_105 = N.N..T = Polymerase basic protein 2 (C5E527_9INFA)  | C5E527_657_662 = N.N..T = Polymerase basic protein 2 (C5E527_9INFA)  | C5E527_737_739 = RKR = Polymerase basic protein 2 (C5E527_9INFA)  | D1LN35_36_41 = L.R..K = Non-structural protein 1 (D1LN35_9INFA)  | D1LN35_55_60 = R.A..E = Non-structural protein 1 (D1LN35_9INFA)  | D1LN35_219_223 = KRK.E = Non-structural protein 1 (D1LN35_9INFA)  | E5LBT9_339_342 = KRKR = Capsid scaffolding protein (E5LBT9_HHV8)  | E5LBT9_339_341 = KRK = Capsid scaffolding protein (E5LBT9_HHV8)  | F5HE15_93_96 = KRKR = Core gene UL26.5 family protein (F5HE15_HHV8)  | F5HFG5_268_272 = L.D.D = Membrane protein (F5HFG5_HHV8)  | I6T1Z2_35_38 = R.RR = Non-structural protein 1 (I6T1Z2_I68A0)  | I6T1Z2_75_77 = EAL = Non-structural protein 1 (I6T1Z2_I68A0)  | I6TAH8_20_25 = Q..T.I = Nucleoprotein (I6TAH8_I68A0)  | I6TAH8_95_99 = P..KR = Nucleoprotein (I6TAH8_I68A0)  | K7Y1A2_207_212 = G.G..T = Genome polyprotein (K7Y1A2_9HEPC)  | O39474_392_396 = S..AS = Non-structural 5A protein (O39474_9HEPC)  | O40939_246_250 = PA..P = ORF K10 (O40939_HHV8)  | O40939_383_387 = PA..P = ORF K10 (O40939_HHV8)  | O40939_94_98 = PA..P = ORF K10 (O40939_HHV8)  | O40939_296_299 = TP[AS][AS] = ORF K10 (O40939_HHV8)  | O56264_217_219 = KQK = Non-structural protein 1 (NS1_I97A1)  | O92837_203_205 = KRK = Minor capsid protein (O92837_SV40)  | P03070_656_658 = DSG = Large T antigen (LT_SV40)  | P03070_129_131 = KRK = Large T antigen (LT_SV40)  | P03070_126_130 = P..KR = Large T antigen (LT_SV40)  | P03070_112_115 = SD.E = Large T antigen (LT_SV40)  | P03070_261_264 = E.EE = Large T antigen (LT_SV40)  | P03070_259_261 = P[DE][DE] = Large T antigen (LT_SV40)  | P03087_5_7 = KRK = Major capsid protein VP1 (VP1_SV40)  | P03101_484_486 = KRK = Major capsid protein L1 (VL1_HPV16)  | P03101_499_501 = KRK = Major capsid protein L1 (VL1_HPV16)  | P03101_502_504 = KRK = Major capsid protein L1 (VL1_HPV16)  | P03126_156_158 = [ST].[LV]$ = Protein E6 (VE6_HPV16)  | P03126_155_158 = [DE]T.[ILV]$ = Protein E6 (VE6_HPV16)  | P03129_33_38 = E.E..I = Protein E7 (VE7_HPV16)  | P03129_34_36 = EED = Protein E7 (VE7_HPV16)  | P03129_41_45 = P..QA = Protein E7 (VE7_HPV16)  | P03129_44_47 = QA.P = Protein E7 (VE7_HPV16)  | P03177_16_20 = R..QE = Thymidine kinase (KITH_EBVB9)  | P03188_272_276 = E..[AGS]F = Envelope glycoprotein B (GB_EBVB9)  | P03209_569_572 = P.TP = Replication and transcription activator (RTA_EBVB9)  | P03211_458_461 = [KR][KR][KR][KR] = Epstein-Barr nuclear antigen 1 (EBNA1_EBVB9)  | P03211_635_638 = E.EE = Epstein-Barr nuclear antigen 1 (EBNA1_EBVB9)  | P03220_66_68 = SSG = Protein BGLF3 (UL95_EBVB9)  | P03225_124_126 = SSG = Protein BDLF2 (BDLF2_EBVB9)  | P03225_45_49 = APR.R = Protein BDLF2 (BDLF2_EBVB9)  | P03243_104_107 = KR.R = E1B 55 kDa protein (E1B55_ADE05)  | P03246_158_162 = R..QE = E1B protein, small T-antigen (E1BS_ADE05)  | P03246_153_159 = Q..A..R = E1B protein, small T-antigen (E1BS_ADE05)  | P03246_153_156 = Q.EA = E1B protein, small T-antigen (E1BS_ADE05)  | P03255_285_288 = KR.R = Early E1A protein (E1A_ADE05)  | P03346_142_145 = P[GS][AS]P = Gag polyprotein (GAG_HTLV2)  | P03346_94_97 = P[GS][AS]P = Gag polyprotein (GAG_HTLV2)  | P03366_107_112 = E..K.K = Gag-Pol polyprotein (POL_HV1B1)  | P03366_217_220 = PV.A = Gag-Pol polyprotein (POL_HV1B1)  | P03366_658_661 = PV.A = Gag-Pol polyprotein (POL_HV1B1)  | P03427_737_739 = RKR = Polymerase basic protein 2 (PB2_I33A0)  | P03427_187_191 = KEK.[DE] = Polymerase basic protein 2 (PB2_I33A0)  | P03427_30_33 = I.KK = Polymerase basic protein 2 (PB2_I33A0)  | P03428_736_739 = KR.R = Polymerase basic protein 2 (PB2_I34A1)  | P03428_752_755 = KR.R = Polymerase basic protein 2 (PB2_I34A1)  | P03428_188_193 = E..E.L = Polymerase basic protein 2 (PB2_I34A1)  | P03428_188_194 = E..E..Q = Polymerase basic protein 2 (PB2_I34A1)  | P03428_69_75 = E..E..Q = Polymerase basic protein 2 (PB2_I34A1)  | P03428_530_533 = T.SS = Polymerase basic protein 2 (PB2_I34A1)  | P03428_191_193 = EEL = Polymerase basic protein 2 (PB2_I34A1)  | P03431_669_672 = KR.R = RNA-directed RNA polymerase catalytic subunit (RDRP_I34A1)  | P03431_205_214 = I..K.{1,2}R.{1,2}K = RNA-directed RNA polymerase catalytic subunit (RDRP_I34A1)  | P03431_751_753 = EEL = RNA-directed RNA polymerase catalytic subunit (RDRP_I34A1)  | P03433_349_352 = E.EE = Polymerase acidic protein (PA_I34A1)  | P03466_20_25 = Q..T.I = Nucleoprotein (NCAP_I34A1)  | P03495_72_77 = E..E.L = Non-structural protein 1 (NS1_I72A2)  | P03495_216_220 = P..KR = Non-structural protein 1 (NS1_I72A2)  | P03495_200_204 = RF..[AGS] = Non-structural protein 1 (NS1_I72A2)  | P03495_217_223 = K..RK.A = Non-structural protein 1 (NS1_I72A2)  | P03495_219_221 = KRK = Non-structural protein 1 (NS1_I72A2)  | P03495_70_75 = K.E..E = Non-structural protein 1 (NS1_I72A2)  | P03496_71_74 = EE.D = Non-structural protein 1 (NS1_I34A1)  | P03496_72_77 = E..E.L = Non-structural protein 1 (NS1_I34A1)  | P03496_216_220 = P..KR = Non-structural protein 1 (NS1_I34A1)  | P03496_217_219 = KQK = Non-structural protein 1 (NS1_I34A1)  | P03496_83_87 = S..AS = Non-structural protein 1 (NS1_I34A1)  | P04012_493_495 = KRK = Major capsid protein L1 (VL1_HPV11)  | P04015_238_242 = [HR][HK]R.R = Regulatory protein E2 (VE2_HPV11)  | P04015_239_242 = [HK][KR].[HKR] = Regulatory protein E2 (VE2_HPV11)  | P04296_1170_1172 = RKR = Major DNA-binding protein (DNBI_HHV11)  | P04296_118_121 = R.RR = Major DNA-binding protein (DNBI_HHV11)  | P04413_145_148 = [KR][KR].[KR] = Serine/threonine-protein kinase US3 (US03_HHV11)  | P04413_199_203 = P.S.G = Serine/threonine-protein kinase US3 (US03_HHV11)  | P04487_7_9 = PAP = RNA-binding protein (RNB_HHV11)  | P04487_54_56 = AAG = RNA-binding protein (RNB_HHV11)  | P04591_107_113 = E..K..K = Gag polyprotein (GAG_HV1H2)  | P05919_51_53 = DSG = Protein Vpu (VPU_HV1H2)  | P06427_147_149 = [ST].[LV]$ = Protein E6 (VE6_HPV33)  | P06427_146_149 = [DE]T.[ILV]$ = Protein E6 (VE6_HPV33)  | P06428_11_15 = DT..D = Protein E6 (VE6_HPV08)  | P06429_8_12 = LKE.V = Protein E7 (VE7_HPV33)  | P06430_31_33 = EEL = Protein E7 (VE7_HPV08)  | P06460_2_6 = [DE]..P[LV] = Probable protein E5A (VE5A_HPV6B)  | P06462_7_9 = STS = Protein E6 (VE6_HPV6B)  | P06463_156_158 = [ST].[LV]$ = Protein E6 (VE6_HPV18)  | P06463_155_158 = [DE]T.[ILV]$ = Protein E6 (VE6_HPV18)  | P06463_148_152 = E..QR = Protein E6 (VE6_HPV18)  | P06463_150_155 = L.R..E = Protein E6 (VE6_HPV18)  | P06463_123_128 = [DE]..[KR].[HKR] = Protein E6 (VE6_HPV18)  | P06464_36_39 = E.DE = Protein E7 (VE7_HPV6B)  | P06464_32_35 = S.ED = Protein E7 (VE7_HPV6B)  | P06464_31_36 = D..E.E = Protein E7 (VE7_HPV6B)  | P06788_33_36 = D.EE = Protein E7 (VE7_HPV18)  | P06821_6_11 = E.E..I = Matrix protein 2 (M2_I34A1)  | P06821_81_85 = Q..[FIMV][DE] = Matrix protein 2 (M2_I34A1)  | P06821_89_92 = G.FV = Matrix protein 2 (M2_I34A1)  | P06821_83_89 = A..A..G = Matrix protein 2 (M2_I34A1)  | P06827_95_99 = P..KR = Nucleoprotein (NCAP_I72A2)  | P06930_13_16 = EK.K = Protein E6 (VE6_HPV05)  | P08393_339_344 = G[KR].RGR = E3 ubiquitin-protein ligase ICP0 (ICP0_HHV11)  | P09992_64_68 = LQ..R = Nucleoprotein (NCAP_LYCVA)  | P0C1C6_298_301 = E.EE = Protein W (W_HENDH)  | P0C1C7_243_246 = EE.D = Protein W (W_NIPAV)  | P0C1C7_320_324 = P..KR = Protein W (W_NIPAV)  | P0C1C7_268_271 = E.EE = Protein W (W_NIPAV)  | P0C213_351_353 = [ST].[LV]$ = Protein Tax-1 (TAX_HTL1F)  | P0C213_350_353 = [DE]T.[ILV]$ = Protein Tax-1 (TAX_HTL1F)  | P0C739_21_23 = GSS = Protein BNLF2a (BNL2A_EBVB9)  | P0C739_23_28 = S..T.P = Protein BNLF2a (BNL2A_EBVB9)  | P0C739_14_18 = SS..G = Protein BNLF2a (BNL2A_EBVB9)  | P0C746_111_115 = [DE]E..R = HTLV-1 basic zipper factor (HBZ_HTL1A)  | P0C746_133_137 = [DE]E..R = HTLV-1 basic zipper factor (HBZ_HTL1A)  | P0C746_54_58 = [DE]E..R = HTLV-1 basic zipper factor (HBZ_HTL1A)  | P0CK56_7_10 = LS.P = Uncharacterized protein BDLF4 (UL92_EBVB9)  | P0CK56_101_105 = R.K.E = Uncharacterized protein BDLF4 (UL92_EBVB9)  | P12418_109_113 = KK..R = Microtubule-associated protein mu-2 (MU2_REOVD)  | P13285_67_69 = N.{0,1}GD = Latent membrane protein 2 (LMP2_EBVB9)  | P16717_146_148 = [ST].[LV]$ = Putative fusion protein (VFUS_SHEVK)  | P16717_146_148 = T.L$ = Putative fusion protein (VFUS_SHEVK)  | P17382_46_49 = D.VD = Replication protein E1 (VE1_HPV31)  | P17386_147_149 = [ST].[LV]$ = Protein E6 (VE6_HPV31)  | P21605_93_100 = S..E[DE]..S = Protein E3 (E3_VACCW)  | P21698_139_144 = D.I..L = Non-structural protein S (NSS_RVFVZ)  | P21698_79_85 = [AS]P..[AS].P = Non-structural protein S (NSS_RVFVZ)  | P21698_134_138 = [DE].[ILV].[KR] = Non-structural protein S (NSS_RVFVZ)  | P21735_156_158 = [ST].[LV]$ = Protein E6 (VE6_HPV45)  | P21735_155_158 = [DE]T.[ILV]$ = Protein E6 (VE6_HPV45)  | P24772_4_6 = NFS = Protein B14 (B14_VACCW)  | P24835_156_158 = [ST].[LV]$ = Protein E6 (VE6_HPV39)  | P24835_155_158 = [DE]T.[ILV]$ = Protein E6 (VE6_HPV39)  | P26555_147_149 = [ST].[LV]$ = Protein E6 (VE6_HPV58)  | P27228_147_149 = [ST].[LV]$ = Protein E6 (VE6_HPV35)  | P27228_146_149 = [DE]T.[ILV]$ = Protein E6 (VE6_HPV35)  | P30119_4_8 = C..PG = Uncharacterized protein BTRF1 (BTRF1_EBVB9)  | P31345_736_739 = KRKR = Polymerase basic protein 2 (PB2_I75A3)  | P36780_203_207 = PP..D = Regulatory protein E2 (VE2_HPV09)  | P36780_277_282 = Q..[HKR].[KR] = Regulatory protein E2 (VE2_HPV09)  | P50804_156_158 = [ST].[LV]$ = Protein E6 (VE6_HPV70)  | P50804_155_158 = [DE]T.[ILV]$ = Protein E6 (VE6_HPV70)  | Q01220_101_106 = KQ.{0,2}K.{1,2}L = Protein A52 (A52_VACCW)  | Q01220_21_25 = EN.[DE]R = Protein A52 (A52_VACCW)  | Q01220_36_41 = T..V.K = Protein A52 (A52_VACCW)  | Q04360_185_192 = SRG.G..T = mRNA export factor ICP27 homolog (ICP27_EBVB9)  | Q05127_201_204 = NN.N = Polymerase cofactor VP35 (VP35_EBOZM)  | Q05322_233_235 = GSS = Membrane-associated protein VP24 (VP24_EBOZM)  | Q0A2H0_228_230 = [ST].V$ = Non-structural protein 1 (NS1_I59A0)  | Q194T2_189_192 = NN.V = Non-structural protein 1 (Q194T2_9INFA)  | Q1K9H2_20_25 = Q..T.I = Nucleoprotein (Q1K9H2_I33A0)  | Q1K9H2_4_8 = K..KR = Nucleoprotein (Q1K9H2_I33A0)  | Q1K9H5_232_237 = E..K.K = RNA-directed RNA polymerase catalytic subunit (Q1K9H5_I33A0)  | Q2MG95_260_262 = EEL = BVLF1 (Q2MG95_EBVB9)  | Q2MG95_260_264 = EE..Q = BVLF1 (Q2MG95_EBVB9)  | Q2PJP0_223_225 = [ST].[ILV]$ = Non-structural protein 1 (Q2PJP0_9INFA)  | Q2PJP0_223_225 = [ST].V$ = Non-structural protein 1 (Q2PJP0_9INFA)  | Q2PJP0_201_204 = I.ED = Non-structural protein 1 (Q2PJP0_9INFA)  | Q2PJP1_22_27 = A..S.D = Nuclear export protein (Q2PJP1_9INFA)  | Q2Q067_176_183 = EA..R..L = BZIP factor (Q2Q067_9DELA)  | Q5EP28_20_25 = Q..T.I = Nucleoprotein (Q5EP28_9INFA)  | Q67296_736_739 = KR.R = Polymerase basic protein 2 (PB2_I72A2)  | Q67296_752_755 = KR.R = Polymerase basic protein 2 (PB2_I72A2)  | Q67296_530_533 = T.SS = Polymerase basic protein 2 (PB2_I72A2)  | Q67296_191_193 = EEL = Polymerase basic protein 2 (PB2_I72A2)  | Q69117_21_28 = V.G..P.G = Tripartite terminase subunit 2 (Q69117_EBVG)  | Q69117_25_31 = P..G..G = Tripartite terminase subunit 2 (Q69117_EBVG)  | Q69117_8_12 = V.S.{0,2}G = Tripartite terminase subunit 2 (Q69117_EBVG)  | Q6DP93_223_225 = [ST].V$ = Non-structural protein 1 (Q6DP93_9INFA)  | Q6VGS8_126_128 = [ST].[LV]$ = 14.3 kDa protein (Q6VGS8_ADE05)  | Q76S40_201_204 = PA.P = Dihydrofolate reductase homolog (Q76S40_HHV8)  | Q77M19_229_237 = K..H.R..S = V protein (Q77M19_MEASW)  | Q77M19_148_150 = SEN = V protein (Q77M19_MEASW)  | Q77M19_98_101 = R.LQ = V protein (Q77M19_MEASW)  | Q77M19_192_195 = E.LR = V protein (Q77M19_MEASW)  | Q77M19_229_233 = KK..R = V protein (Q77M19_MEASW)  | Q77M19_132_139 = Q..L..[DE]S = V protein (Q77M19_MEASW)  | Q77M19_93_97 = L..PP = V protein (Q77M19_MEASW)  | Q77M19_75_79 = APR.R = V protein (Q77M19_MEASW)  | Q77M19_59_63 = [GS][ILM][AS].P = V protein (Q77M19_MEASW)  | Q77M19_85_88 = E.DD = V protein (Q77M19_MEASW)  | Q77M19_80_84 = G.G.G = V protein (Q77M19_MEASW)  | Q77Q36_8_16 = P..[LV]..P.L = viral cyclin homolog (VCYCL_HHV8P)  | Q77UU1_237_241 = [HR][HK]R.R = Regulatory protein E2 (Q77UU1_9PAPI)  | Q8AZK7_156_161 = L.Q.L.{0,1}R = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_156_162 = L.Q.L.{0,1}R = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_222_227 = L.Q.L.{0,1}R = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_222_228 = L.Q.L.{0,1}R = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_24_29 = L.Q.L.{0,1}R = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_24_30 = L.Q.L.{0,1}R = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_288_293 = L.Q.L.{0,1}R = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_288_294 = L.Q.L.{0,1}R = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_354_359 = L.Q.L.{0,1}R = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_354_360 = L.Q.L.{0,1}R = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_420_425 = L.Q.L.{0,1}R = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_420_426 = L.Q.L.{0,1}R = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_90_95 = L.Q.L.{0,1}R = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_90_96 = L.Q.L.{0,1}R = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_133_139 = R..R..G = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_199_205 = R..R..G = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_265_271 = R..R..G = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_331_337 = R..R..G = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_397_403 = R..R..G = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_67_73 = R..R..G = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_13_20 = P..PG..P = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_145_152 = P..PG..P = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_211_218 = P..PG..P = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_277_284 = P..PG..P = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_343_350 = P..PG..P = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_409_416 = P..PG..P = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_79_86 = P..PG..P = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_10_14 = P..PG = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_13_17 = P..PG = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_142_146 = P..PG = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_145_149 = P..PG = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_208_212 = P..PG = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_211_215 = P..PG = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_274_278 = P..PG = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_277_281 = P..PG = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_340_344 = P..PG = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_343_347 = P..PG = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_406_410 = P..PG = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_409_413 = P..PG = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_76_80 = P..PG = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_79_83 = P..PG = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_148_152 = PG..P = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_16_20 = PG..P = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_214_218 = PG..P = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_280_284 = PG..P = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_346_350 = PG..P = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_412_416 = PG..P = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_82_86 = PG..P = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_126_130 = S.SP[AS] = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_192_196 = S.SP[AS] = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_258_262 = S.SP[AS] = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_324_328 = S.SP[AS] = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_390_394 = S.SP[AS] = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_456_460 = S.SP[AS] = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_60_64 = S.SP[AS] = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_104_108 = R..QE = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_170_174 = R..QE = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_236_240 = R..QE = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_302_306 = R..QE = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_368_372 = R..QE = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_38_42 = R..QE = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_434_438 = R..QE = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_133_139 = RG.R..G = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_199_205 = RG.R..G = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_265_271 = RG.R..G = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_331_337 = RG.R..G = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_397_403 = RG.R..G = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_67_73 = RG.R..G = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q8AZK7_504_506 = DED = Epstein-Barr nuclear antigen leader protein (EBNA5_EBVB9)  | Q997F2_320_324 = P..KR = Non-structural protein V (V_NIPAV)  | Q99AU3_36_41 = L.R..K = Non-structural protein 1 (NS1_I18A0)  | Q99AU3_220_223 = RK.A = Non-structural protein 1 (NS1_I18A0)  | Q9DGW5_296_301 = E..R.Y = Oncoprotein MEQ (MEQ_GAHVM)  | Q9DGW5_90_94 = [DE]E..R = Oncoprotein MEQ (MEQ_GAHVM)  | Q9QPN3_72_77 = P..P.[HKR] = Protein Nef (NEF_HV1LA)  | Q9WMB5_464_467 = AA.P = Nucleoprotein (NCAP_MEASC)  | Q9WMX2_2273_2276 = E.LR = Genome polyprotein [Cleaved into: Core protein p21 (POLG_HCVCO)  | Q9WPI5_2_6 = EE..R = Nonstructural protein NS1 (Q9WPI5_9INFA)  | U5TQE9_57_63 = D..P..A = Ribonucleoside-diphosphate reductase large subunit (U5TQE9_HHV1) 
## 276
## 
## viral targeted human proteins: P52294 = Importin subunit alpha-5 | Q15051 = IQ calmodulin-binding motif-containing protein 1 | Q9NTJ3 = Structural maintenance of chromosomes protein 4 | O60812 = Heterogeneous nuclear ribonucleoprotein C-like 1 | Q6PKG0 = La-related protein 1 | Q6R327 = Rapamycin-insensitive companion of mTOR | Q13043 = Serine/threonine-protein kinase 4 | Q00839 = Heterogeneous nuclear ribonucleoprotein U | O60684 = Importin subunit alpha-7 | P52292 = Importin subunit alpha-1 | P38606 = V-type proton ATPase catalytic subunit A | Q14683 = Structural maintenance of chromosomes protein 1A | P07237 = Protein disulfide-isomerase | Q02539 = Histone H1.1 | Q14164 = Inhibitor of nuclear factor kappa-B kinase subunit epsilon | Q13451 = Peptidyl-prolyl cis-trans isomerase FKBP5 | Q12873 = Chromodomain-helicase-DNA-binding protein 3 | Q13042 = Cell division cycle protein 16 homolog | P42336 = Phosphatidylinositol 4,5-bisphosphate 3-kinase catalytic subunit alpha isoform | Q9UKB1 = F-box/WD repeat-containing protein 11 | Q9Y297 = F-box/WD repeat-containing protein 1A | O00505 = Importin subunit alpha-4 | P06493 = Cyclin-dependent kinase 1 | O43683 = Mitotic checkpoint serine/threonine-protein kinase BUB1 | Q12959 = Disks large homolog 1 | Q14160 = Protein scribble homolog | P24278 = Zinc finger and BTB domain-containing protein 25 | P62191 = 26S proteasome regulatory subunit 4 | P30153 = Serine/threonine-protein phosphatase 2A 65 kDa regulatory subunit A alpha isoform | P54646 = 5'-AMP-activated protein kinase catalytic subunit alpha-2 | Q9Y5M8 = Signal recognition particle receptor subunit beta | Q13077 = TNF receptor-associated factor 1 | O00268 = Transcription initiation factor TFIID subunit 4 | Q14192 = Four and a half LIM domains protein 2 | Q9Y230 = RuvB-like 2 | P29590 = Protein PML | P05556 = Integrin beta-1 | P51148 = Ras-related protein Rab-5C | Q53G59 = Kelch-like protein 12 | P33991 = DNA replication licensing factor MCM4 | Q99613 = Eukaryotic translation initiation factor 3 subunit C | P50990 = T-complex protein 1 subunit theta | Q92997 = Segment polarity protein dishevelled homolog DVL-3 | P45984 = Mitogen-activated protein kinase 9 | Q8N448 = Ligand of Numb protein X 2 | P04049 = RAF proto-oncogene serine/threonine-protein kinase | P62805 = Histone H4 | Q9H190 = Syntenin-2 | O95793 = Double-stranded RNA-binding protein Staufen homolog 1 | P27348 = 14-3-3 protein theta | Q07955 = Serine/arginine-rich splicing factor 1 | O60885 = Bromodomain-containing protein 4 | P09874 = Poly [ADP-ribose] polymerase 1 | P54132 = Bloom syndrome protein | P68104 = Elongation factor 1-alpha 1 | P78371 = T-complex protein 1 subunit beta | P26196 = Probable ATP-dependent RNA helicase DDX6 | Q96I25 = Splicing factor 45 | Q9UHD8 = Septin-9 | Q9H9D4 = Zinc finger protein 408 | Q9P0L0 = Vesicle-associated membrane protein-associated protein A | Q09472 = Histone acetyltransferase p300 | P17980 = 26S proteasome regulatory subunit 6A | P33993 = DNA replication licensing factor MCM7 | Q9NR12 = PDZ and LIM domain protein 7 | Q04864 = Proto-oncogene c-Rel | P85037 = Forkhead box protein K1 | P24941 = Cyclin-dependent kinase 2 | P61962 = DDB1- and CUL4-associated factor 7 | P50570 = Dynamin-2 | O15379 = Histone deacetylase 3 | P49368 = T-complex protein 1 subunit gamma | O95817 = BAG family molecular chaperone regulator 3 | P17544 = Cyclic AMP-dependent transcription factor ATF-7 | P53621 = Coatomer subunit alpha | Q71U36 = Tubulin alpha-1A chain | P37173 = TGF-beta receptor type-2 | P42574 = Caspase-3 | O60506 = Heterogeneous nuclear ribonucleoprotein Q | P48643 = T-complex protein 1 subunit epsilon | P50991 = T-complex protein 1 subunit delta | O14920 = Inhibitor of nuclear factor kappa-B kinase subunit beta | Q9NRL2 = Bromodomain adjacent to zinc finger domain protein 1A | Q01130 = Serine/arginine-rich splicing factor 2 | P38159 = RNA-binding motif protein, X chromosome | P55795 = Heterogeneous nuclear ribonucleoprotein H2 | P68431 = Histone H3.1 | P78352 = Disks large homolog 4 | Q92569 = Phosphatidylinositol 3-kinase regulatory subunit gamma | P20618 = Proteasome subunit beta type-1 | P15822 = Zinc finger protein 40 | P31249 = Homeobox protein Hox-D3 | P23458 = Tyrosine-protein kinase JAK1 | P63151 = Serine/threonine-protein phosphatase 2A 55 kDa regulatory subunit B alpha isoform | P41252 = Isoleucine--tRNA ligase, cytoplasmic | P62491 = Ras-related protein Rab-11A | Q16531 = DNA damage-binding protein 1 | Q09028 = Histone-binding protein RBBP4 | P11802 = Cyclin-dependent kinase 4 | Q9Y265 = RuvB-like 1 | O75533 = Splicing factor 3B subunit 1 | P47897 = Glutamine--tRNA ligase | P61978 = Heterogeneous nuclear ribonucleoprotein K | Q8N1F7 = Nuclear pore complex protein Nup93 | P15336 = Cyclic AMP-dependent transcription factor ATF-2 | P12931 = Proto-oncogene tyrosine-protein kinase Src | P11940 = Polyadenylate-binding protein 1 | P07910 = Heterogeneous nuclear ribonucleoproteins C1/C2 | Q9Y572 = Receptor-interacting serine/threonine-protein kinase 3
## 109
## 
## UP support / UP: (B4URF7_736_739 in P52294 dataset) 13/50 | (B4URF7_188_194 in Q15051 dataset) 27/122 | (B4URF7_69_75 in Q15051 dataset) 27/122 | (B4URF7_737_739 in Q9NTJ3 dataset) 12/26 | (C5E522_20_25 in O60812 dataset) 7/15 | (C5E522_443_445 in Q6PKG0 dataset) 12/38 | (C5E526_205_214 in Q6R327 dataset) 7/20 | (C5E527_100_105 in Q13043 dataset) 14/83 | (C5E527_657_662 in Q13043 dataset) 14/83 | (C5E527_158_163 in Q13043 dataset) 20/83 | (C5E527_737_739 in Q9NTJ3 dataset) 12/26 | (D1LN35_36_41 in Q00839 dataset) 34/165 | (D1LN35_219_223 in Q00839 dataset) 12/165 | (D1LN35_55_60 in Q00839 dataset) 30/165 | (E5LBT9_339_342 in O60684 dataset) 12/39 | (E5LBT9_339_341 in P52292 dataset) 24/84 | (E5LBT9_339_342 in P52294 dataset) 13/50 | (F5HE15_93_96 in O60684 dataset) 12/39 | (F5HE15_93_96 in P52294 dataset) 13/50 | (F5HFG5_268_272 in P38606 dataset) 13/57 | (I6T1Z2_75_77 in P07237 dataset) 16/53 | (I6T1Z2_35_38 in Q14683 dataset) 16/50 | (I6TAH8_20_25 in O60812 dataset) 7/15 | (I6TAH8_95_99 in Q02539 dataset) 9/21 | (K7Y1A2_207_212 in Q14164 dataset) 27/216 | (O39474_392_396 in Q13451 dataset) 21/54 | (O40939_383_387 in Q12873 dataset) 21/68 | (O40939_94_98 in Q12873 dataset) 21/68 | (O40939_246_250 in Q12873 dataset) 21/68 | (O40939_296_299 in Q13042 dataset) 8/33 | (O56264_217_219 in P42336 dataset) 10/35 | (O92837_203_205 in P52292 dataset) 24/84 | (P03070_261_264 in O00505 dataset) 27/62 | (P03070_129_131 in O00505 dataset) 21/62 | (P03070_259_261 in O43683 dataset) 32/56 | (P03070_129_131 in O60684 dataset) 21/39 | (P03070_112_115 in P06493 dataset) 19/54 | (P03070_129_131 in P52292 dataset) 24/84 | (P03070_129_131 in P52294 dataset) 21/50 | (P03070_126_130 in Q02539 dataset) 9/21 | (P03070_656_658 in Q9UKB1 dataset) 24/36 | (P03070_656_658 in Q9Y297 dataset) 29/56 | (P03087_5_7 in P52292 dataset) 24/84 | (P03101_499_501 in P52292 dataset) 24/84 | (P03101_484_486 in P52292 dataset) 24/84 | (P03101_502_504 in P52292 dataset) 24/84 | (P03126_156_158 in Q12959 dataset) 32/46 | (P03126_155_158 in Q14160 dataset) 23/56 | (P03129_33_38 in P24278 dataset) 11/27 | (P03129_41_45 in P30153 dataset) 21/112 | (P03129_44_47 in P54646 dataset) 19/74 | (P03129_34_36 in P62191 dataset) 20/67 | (P03177_16_20 in Q9Y5M8 dataset) 11/48 | (P03188_272_276 in Q13077 dataset) 32/135 | (P03209_569_572 in O00268 dataset) 11/16 | (P03211_635_638 in O00505 dataset) 27/62 | (P03211_458_461 in O60684 dataset) 19/39 | (P03211_458_461 in P52294 dataset) 21/50 | (P03220_66_68 in Q14192 dataset) 29/70 | (P03225_124_126 in Q14192 dataset) 29/70 | (P03225_45_49 in Q9Y230 dataset) 8/91 | (P03243_104_107 in P29590 dataset) 17/64 | (P03246_153_159 in P05556 dataset) 11/44 | (P03246_153_156 in P51148 dataset) 19/90 | (P03246_158_162 in Q9Y5M8 dataset) 11/48 | (P03255_285_288 in P29590 dataset) 17/64 | (P03346_94_97 in Q53G59 dataset) 16/36 | (P03346_142_145 in Q53G59 dataset) 16/36 | (P03366_107_112 in P33991 dataset) 21/42 | (P03366_658_661 in Q99613 dataset) 10/30 | (P03366_217_220 in Q99613 dataset) 10/30 | (P03427_737_739 in O00505 dataset) 22/62 | (P03427_30_33 in P50990 dataset) 12/42 | (P03427_187_191 in P62191 dataset) 11/67 | (P03428_188_193 in P45984 dataset) 23/51 | (P03428_188_194 in Q15051 dataset) 27/122 | (P03428_69_75 in Q15051 dataset) 27/122 | (P03428_191_193 in Q8N448 dataset) 20/59 | (P03428_752_755 in Q92997 dataset) 23/70 | (P03428_736_739 in Q92997 dataset) 23/70 | (P03428_530_533 in Q92997 dataset) 20/70 | (P03431_205_214 in Q6R327 dataset) 7/20 | (P03431_751_753 in Q8N448 dataset) 20/59 | (P03431_669_672 in Q92997 dataset) 23/70 | (P03433_349_352 in O00505 dataset) 27/62 | (P03466_20_25 in O60812 dataset) 7/15 | (P03495_70_75 in O95793 dataset) 16/46 | (P03495_200_204 in P04049 dataset) 21/70 | (P03495_72_77 in P45984 dataset) 23/51 | (P03495_217_223 in P62805 dataset) 7/62 | (P03495_216_220 in Q02539 dataset) 9/21 | (P03495_219_221 in Q9H190 dataset) 13/45 | (P03496_71_74 in P27348 dataset) 25/63 | (P03496_217_219 in P42336 dataset) 10/35 | (P03496_72_77 in P45984 dataset) 23/51 | (P03496_216_220 in Q02539 dataset) 9/21 | (P03496_83_87 in Q13451 dataset) 21/54 | (P04012_493_495 in P52292 dataset) 24/84 | (P04015_239_242 in O60885 dataset) 17/20 | (P04015_238_242 in Q07955 dataset) 8/24 | (P04296_1170_1172 in P09874 dataset) 21/65 | (P04296_118_121 in P54132 dataset) 12/33 | (P04413_145_148 in P68104 dataset) 58/129 | (P04413_199_203 in P78371 dataset) 18/48 | (P04487_54_56 in P26196 dataset) 20/60 | (P04487_7_9 in Q00839 dataset) 38/165 | (P04591_107_113 in Q96I25 dataset) 13/23 | (P05919_51_53 in Q9UKB1 dataset) 24/36 | (P05919_51_53 in Q9Y297 dataset) 29/56 | (P06427_147_149 in Q12959 dataset) 32/46 | (P06427_146_149 in Q14160 dataset) 23/56 | (P06428_11_15 in Q9UHD8 dataset) 7/18 | (P06429_8_12 in Q9H9D4 dataset) 6/33 | (P06430_31_33 in Q8N448 dataset) 20/59 | (P06460_2_6 in Q9P0L0 dataset) 40/119 | (P06462_7_9 in Q09472 dataset) 26/90 | (P06463_150_155 in P17980 dataset) 14/43 | (P06463_148_152 in P30153 dataset) 23/112 | (P06463_123_128 in P33993 dataset) 60/80 | (P06463_156_158 in Q12959 dataset) 32/46 | (P06463_155_158 in Q14160 dataset) 23/56 | (P06464_31_36 in P85037 dataset) 18/33 | (P06464_32_35 in Q04864 dataset) 23/134 | (P06464_36_39 in Q9NR12 dataset) 21/64 | (P06788_33_36 in P24941 dataset) 20/72 | (P06821_6_11 in P24278 dataset) 11/27 | (P06821_83_89 in P50570 dataset) 14/48 | (P06821_81_85 in P61962 dataset) 25/45 | (P06821_89_92 in Q8N448 dataset) 8/59 | (P06827_95_99 in O00505 dataset) 20/62 | (P06930_13_16 in Q9Y5M8 dataset) 10/48 | (P08393_339_344 in O15379 dataset) 5/38 | (P09992_64_68 in P49368 dataset) 17/72 | (P0C1C6_298_301 in O00505 dataset) 27/62 | (P0C1C7_320_324 in O00505 dataset) 20/62 | (P0C1C7_268_271 in O00505 dataset) 27/62 | (P0C1C7_243_246 in P27348 dataset) 25/63 | (P0C1C7_320_324 in Q02539 dataset) 9/21 | (P0C213_351_353 in Q12959 dataset) 32/46 | (P0C213_350_353 in Q14160 dataset) 23/56 | (P0C739_21_23 in O95817 dataset) 23/72 | (P0C739_23_28 in P51148 dataset) 20/90 | (P0C739_14_18 in Q9P0L0 dataset) 27/119 | (P0C746_111_115 in P17544 dataset) 11/13 | (P0C746_54_58 in P17544 dataset) 11/13 | (P0C746_133_137 in P17544 dataset) 11/13 | (P0CK56_7_10 in P53621 dataset) 14/39 | (P0CK56_101_105 in Q14192 dataset) 17/70 | (P12418_109_113 in Q71U36 dataset) 19/94 | (P13285_67_69 in P37173 dataset) 9/20 | (P16717_146_148 in Q12959 dataset) 32/46 | (P16717_146_148 in Q14160 dataset) 23/56 | (P17382_46_49 in P42574 dataset) 11/32 | (P17386_147_149 in Q12959 dataset) 32/46 | (P21605_93_100 in O60506 dataset) 16/100 | (P21698_139_144 in P48643 dataset) 12/62 | (P21698_79_85 in P50991 dataset) 12/50 | (P21698_134_138 in P78371 dataset) 39/48 | (P21735_156_158 in Q12959 dataset) 32/46 | (P21735_155_158 in Q14160 dataset) 23/56 | (P24772_4_6 in O14920 dataset) 9/50 | (P24835_156_158 in Q12959 dataset) 32/46 | (P24835_155_158 in Q14160 dataset) 23/56 | (P26555_147_149 in Q12959 dataset) 32/46 | (P27228_147_149 in Q12959 dataset) 32/46 | (P27228_146_149 in Q14160 dataset) 23/56 | (P30119_4_8 in P62191 dataset) 11/67 | (P31345_736_739 in O60684 dataset) 12/39 | (P31345_736_739 in P52294 dataset) 13/50 | (P36780_277_282 in Q01130 dataset) 14/20 | (P36780_203_207 in Q9NRL2 dataset) 10/18 | (P50804_156_158 in Q12959 dataset) 32/46 | (P50804_155_158 in Q14160 dataset) 23/56 | (Q01220_101_106 in P42336 dataset) 8/35 | (Q01220_21_25 in P42336 dataset) 6/35 | (Q01220_36_41 in P68104 dataset) 17/129 | (Q04360_185_192 in P38159 dataset) 4/47 | (Q05127_201_204 in P55795 dataset) 9/46 | (Q05322_233_235 in P68431 dataset) 16/37 | (Q0A2H0_228_230 in P78352 dataset) 7/16 | (Q194T2_189_192 in Q92569 dataset) 12/77 | (Q1K9H2_4_8 in O60684 dataset) 15/39 | (Q1K9H2_20_25 in O60812 dataset) 7/15 | (Q1K9H5_232_237 in P33991 dataset) 21/42 | (Q2MG95_260_264 in P20618 dataset) 19/61 | (Q2MG95_260_262 in P51148 dataset) 31/90 | (Q2PJP0_201_204 in P78352 dataset) 9/16 | (Q2PJP0_223_225 in P78352 dataset) 7/16 | (Q2PJP0_223_225 in Q14160 dataset) 32/56 | (Q2PJP1_22_27 in P15822 dataset) 12/43 | (Q2Q067_176_183 in P31249 dataset) 4/5 | (Q5EP28_20_25 in O60812 dataset) 7/15 | (Q67296_191_193 in Q8N448 dataset) 20/59 | (Q67296_752_755 in Q92997 dataset) 23/70 | (Q67296_530_533 in Q92997 dataset) 20/70 | (Q67296_736_739 in Q92997 dataset) 23/70 | (Q69117_21_28 in P23458 dataset) 5/30 | (Q69117_8_12 in P41252 dataset) 15/20 | (Q69117_25_31 in P63151 dataset) 15/42 | (Q6DP93_223_225 in P78352 dataset) 7/16 | (Q6VGS8_126_128 in Q12959 dataset) 32/46 | (Q76S40_201_204 in P63151 dataset) 17/42 | (Q77M19_192_195 in P50991 dataset) 19/50 | (Q77M19_80_84 in P52292 dataset) 24/84 | (Q77M19_229_233 in P52292 dataset) 20/84 | (Q77M19_59_63 in P53621 dataset) 11/39 | (Q77M19_132_139 in P53621 dataset) 7/39 | (Q77M19_229_237 in P62491 dataset) 5/35 | (Q77M19_85_88 in Q09028 dataset) 17/59 | (Q77M19_98_101 in Q16531 dataset) 18/54 | (Q77M19_148_150 in Q9P0L0 dataset) 24/119 | (Q77M19_93_97 in Q9P0L0 dataset) 25/119 | (Q77M19_75_79 in Q9Y230 dataset) 8/91 | (Q77Q36_8_16 in P11802 dataset) 9/56 | (Q77UU1_237_241 in Q07955 dataset) 8/24 | (Q8AZK7_456_460 in O60506 dataset) 16/100 | (Q8AZK7_199_205 in O60506 dataset) 8/100 | (Q8AZK7_504_506 in O60506 dataset) 25/100 | (Q8AZK7_60_64 in O60506 dataset) 16/100 | (Q8AZK7_343_350 in O60506 dataset) 10/100 | (Q8AZK7_192_196 in O60506 dataset) 16/100 | (Q8AZK7_211_218 in O60506 dataset) 10/100 | (Q8AZK7_13_20 in O60506 dataset) 10/100 | (Q8AZK7_390_394 in O60506 dataset) 16/100 | (Q8AZK7_409_416 in O60506 dataset) 10/100 | (Q8AZK7_277_284 in O60506 dataset) 10/100 | (Q8AZK7_126_130 in O60506 dataset) 16/100 | (Q8AZK7_258_262 in O60506 dataset) 16/100 | (Q8AZK7_145_152 in O60506 dataset) 10/100 | (Q8AZK7_324_328 in O60506 dataset) 16/100 | (Q8AZK7_79_86 in O60506 dataset) 10/100 | (Q8AZK7_133_139 in O60506 dataset) 8/100 | (Q8AZK7_265_271 in O60506 dataset) 8/100 | (Q8AZK7_67_73 in O60506 dataset) 8/100 | (Q8AZK7_397_403 in O60506 dataset) 8/100 | (Q8AZK7_331_337 in O60506 dataset) 8/100 | (Q8AZK7_265_271 in O75533 dataset) 15/47 | (Q8AZK7_199_205 in O75533 dataset) 15/47 | (Q8AZK7_397_403 in O75533 dataset) 15/47 | (Q8AZK7_67_73 in O75533 dataset) 15/47 | (Q8AZK7_331_337 in O75533 dataset) 15/47 | (Q8AZK7_133_139 in O75533 dataset) 15/47 | (Q8AZK7_10_14 in P47897 dataset) 20/64 | (Q8AZK7_145_149 in P47897 dataset) 20/64 | (Q8AZK7_79_83 in P47897 dataset) 20/64 | (Q8AZK7_343_347 in P47897 dataset) 20/64 | (Q8AZK7_142_146 in P47897 dataset) 20/64 | (Q8AZK7_13_17 in P47897 dataset) 20/64 | (Q8AZK7_208_212 in P47897 dataset) 20/64 | (Q8AZK7_211_215 in P47897 dataset) 20/64 | (Q8AZK7_340_344 in P47897 dataset) 20/64 | (Q8AZK7_277_281 in P47897 dataset) 20/64 | (Q8AZK7_76_80 in P47897 dataset) 20/64 | (Q8AZK7_274_278 in P47897 dataset) 20/64 | (Q8AZK7_409_413 in P47897 dataset) 20/64 | (Q8AZK7_406_410 in P47897 dataset) 20/64 | (Q8AZK7_346_350 in P61978 dataset) 19/55 | (Q8AZK7_412_416 in P61978 dataset) 19/55 | (Q8AZK7_280_284 in P61978 dataset) 19/55 | (Q8AZK7_82_86 in P61978 dataset) 19/55 | (Q8AZK7_148_152 in P61978 dataset) 19/55 | (Q8AZK7_16_20 in P61978 dataset) 19/55 | (Q8AZK7_214_218 in P61978 dataset) 19/55 | (Q8AZK7_38_42 in Q8N1F7 dataset) 11/40 | (Q8AZK7_104_108 in Q8N1F7 dataset) 11/40 | (Q8AZK7_434_438 in Q8N1F7 dataset) 11/40 | (Q8AZK7_170_174 in Q8N1F7 dataset) 11/40 | (Q8AZK7_368_372 in Q8N1F7 dataset) 11/40 | (Q8AZK7_236_240 in Q8N1F7 dataset) 11/40 | (Q8AZK7_302_306 in Q8N1F7 dataset) 11/40 | (Q8AZK7_354_359 in Q9Y265 dataset) 11/80 | (Q8AZK7_354_360 in Q9Y265 dataset) 11/80 | (Q8AZK7_156_162 in Q9Y265 dataset) 11/80 | (Q8AZK7_24_30 in Q9Y265 dataset) 11/80 | (Q8AZK7_420_425 in Q9Y265 dataset) 11/80 | (Q8AZK7_288_294 in Q9Y265 dataset) 11/80 | (Q8AZK7_222_227 in Q9Y265 dataset) 11/80 | (Q8AZK7_90_95 in Q9Y265 dataset) 11/80 | (Q8AZK7_24_29 in Q9Y265 dataset) 11/80 | (Q8AZK7_156_161 in Q9Y265 dataset) 11/80 | (Q8AZK7_288_293 in Q9Y265 dataset) 11/80 | (Q8AZK7_222_228 in Q9Y265 dataset) 11/80 | (Q8AZK7_420_426 in Q9Y265 dataset) 11/80 | (Q8AZK7_90_96 in Q9Y265 dataset) 11/80 | (Q997F2_320_324 in Q02539 dataset) 9/21 | (Q99AU3_220_223 in Q00839 dataset) 35/165 | (Q99AU3_36_41 in Q00839 dataset) 34/165 | (Q9DGW5_296_301 in P15336 dataset) 14/105 | (Q9DGW5_90_94 in P17544 dataset) 11/13 | (Q9QPN3_72_77 in P12931 dataset) 36/78 | (Q9WMB5_464_467 in P30153 dataset) 28/112 | (Q9WMX2_2273_2276 in P50991 dataset) 19/50 | (Q9WPI5_2_6 in P07910 dataset) 18/52 | (Q9WPI5_2_6 in P11940 dataset) 21/55 | (U5TQE9_57_63 in Q9Y572 dataset) 11/39
## 
## Occurence / SeqNum: (B4URF7_736_739 in P52294 dataset) 23/82 | (B4URF7_188_194 in Q15051 dataset) 54/195 | (B4URF7_69_75 in Q15051 dataset) 54/195 | (B4URF7_737_739 in Q9NTJ3 dataset) 23/39 | (C5E522_20_25 in O60812 dataset) 12/28 | (C5E522_443_445 in Q6PKG0 dataset) 16/69 | (C5E526_205_214 in Q6R327 dataset) 8/39 | (C5E527_100_105 in Q13043 dataset) 17/172 | (C5E527_657_662 in Q13043 dataset) 17/172 | (C5E527_158_163 in Q13043 dataset) 39/172 | (C5E527_737_739 in Q9NTJ3 dataset) 23/39 | (D1LN35_36_41 in Q00839 dataset) 55/303 | (D1LN35_219_223 in Q00839 dataset) 20/303 | (D1LN35_55_60 in Q00839 dataset) 51/303 | (E5LBT9_339_342 in O60684 dataset) 17/55 | (E5LBT9_339_341 in P52292 dataset) 41/125 | (E5LBT9_339_342 in P52294 dataset) 23/82 | (F5HE15_93_96 in O60684 dataset) 17/55 | (F5HE15_93_96 in P52294 dataset) 23/82 | (F5HFG5_268_272 in P38606 dataset) 14/81 | (I6T1Z2_75_77 in P07237 dataset) 18/71 | (I6T1Z2_35_38 in Q14683 dataset) 33/94 | (I6TAH8_20_25 in O60812 dataset) 12/28 | (I6TAH8_95_99 in Q02539 dataset) 11/31 | (K7Y1A2_207_212 in Q14164 dataset) 33/351 | (O39474_392_396 in Q13451 dataset) 31/116 | (O40939_383_387 in Q12873 dataset) 45/107 | (O40939_94_98 in Q12873 dataset) 45/107 | (O40939_246_250 in Q12873 dataset) 45/107 | (O40939_296_299 in Q13042 dataset) 9/42 | (O56264_217_219 in P42336 dataset) 12/50 | (O92837_203_205 in P52292 dataset) 41/125 | (P03070_261_264 in O00505 dataset) 37/84 | (P03070_129_131 in O00505 dataset) 30/84 | (P03070_259_261 in O43683 dataset) 63/81 | (P03070_129_131 in O60684 dataset) 36/55 | (P03070_112_115 in P06493 dataset) 25/90 | (P03070_129_131 in P52292 dataset) 41/125 | (P03070_129_131 in P52294 dataset) 42/82 | (P03070_126_130 in Q02539 dataset) 11/31 | (P03070_656_658 in Q9UKB1 dataset) 43/49 | (P03070_656_658 in Q9Y297 dataset) 55/85 | (P03087_5_7 in P52292 dataset) 41/125 | (P03101_499_501 in P52292 dataset) 41/125 | (P03101_484_486 in P52292 dataset) 41/125 | (P03101_502_504 in P52292 dataset) 41/125 | (P03126_156_158 in Q12959 dataset) 60/91 | (P03126_155_158 in Q14160 dataset) 32/102 | (P03129_33_38 in P24278 dataset) 11/36 | (P03129_41_45 in P30153 dataset) 32/218 | (P03129_44_47 in P54646 dataset) 24/136 | (P03129_34_36 in P62191 dataset) 29/82 | (P03177_16_20 in Q9Y5M8 dataset) 12/59 | (P03188_272_276 in Q13077 dataset) 52/261 | (P03209_569_572 in O00268 dataset) 14/20 | (P03211_635_638 in O00505 dataset) 37/84 | (P03211_458_461 in O60684 dataset) 45/55 | (P03211_458_461 in P52294 dataset) 62/82 | (P03220_66_68 in Q14192 dataset) 53/131 | (P03225_124_126 in Q14192 dataset) 53/131 | (P03225_45_49 in Q9Y230 dataset) 8/134 | (P03243_104_107 in P29590 dataset) 20/82 | (P03246_153_159 in P05556 dataset) 11/61 | (P03246_153_156 in P51148 dataset) 28/115 | (P03246_158_162 in Q9Y5M8 dataset) 12/59 | (P03255_285_288 in P29590 dataset) 20/82 | (P03346_94_97 in Q53G59 dataset) 57/69 | (P03346_142_145 in Q53G59 dataset) 57/69 | (P03366_107_112 in P33991 dataset) 34/65 | (P03366_658_661 in Q99613 dataset) 14/40 | (P03366_217_220 in Q99613 dataset) 14/40 | (P03427_737_739 in O00505 dataset) 32/84 | (P03427_30_33 in P50990 dataset) 13/94 | (P03427_187_191 in P62191 dataset) 12/82 | (P03428_188_193 in P45984 dataset) 53/82 | (P03428_188_194 in Q15051 dataset) 54/195 | (P03428_69_75 in Q15051 dataset) 54/195 | (P03428_191_193 in Q8N448 dataset) 26/84 | (P03428_752_755 in Q92997 dataset) 61/185 | (P03428_736_739 in Q92997 dataset) 61/185 | (P03428_530_533 in Q92997 dataset) 52/185 | (P03431_205_214 in Q6R327 dataset) 8/39 | (P03431_751_753 in Q8N448 dataset) 26/84 | (P03431_669_672 in Q92997 dataset) 61/185 | (P03433_349_352 in O00505 dataset) 37/84 | (P03466_20_25 in O60812 dataset) 12/28 | (P03495_70_75 in O95793 dataset) 26/81 | (P03495_200_204 in P04049 dataset) 25/120 | (P03495_72_77 in P45984 dataset) 53/82 | (P03495_217_223 in P62805 dataset) 11/126 | (P03495_216_220 in Q02539 dataset) 11/31 | (P03495_219_221 in Q9H190 dataset) 28/69 | (P03496_71_74 in P27348 dataset) 52/161 | (P03496_217_219 in P42336 dataset) 12/50 | (P03496_72_77 in P45984 dataset) 53/82 | (P03496_216_220 in Q02539 dataset) 11/31 | (P03496_83_87 in Q13451 dataset) 31/116 | (P04012_493_495 in P52292 dataset) 41/125 | (P04015_239_242 in O60885 dataset) 73/42 | (P04015_238_242 in Q07955 dataset) 12/68 | (P04296_1170_1172 in P09874 dataset) 40/102 | (P04296_118_121 in P54132 dataset) 17/39 | (P04413_145_148 in P68104 dataset) 160/179 | (P04413_199_203 in P78371 dataset) 27/110 | (P04487_54_56 in P26196 dataset) 32/108 | (P04487_7_9 in Q00839 dataset) 103/303 | (P04591_107_113 in Q96I25 dataset) 20/34 | (P05919_51_53 in Q9UKB1 dataset) 43/49 | (P05919_51_53 in Q9Y297 dataset) 55/85 | (P06427_147_149 in Q12959 dataset) 60/91 | (P06427_146_149 in Q14160 dataset) 32/102 | (P06428_11_15 in Q9UHD8 dataset) 9/21 | (P06429_8_12 in Q9H9D4 dataset) 6/75 | (P06430_31_33 in Q8N448 dataset) 26/84 | (P06460_2_6 in Q9P0L0 dataset) 76/180 | (P06462_7_9 in Q09472 dataset) 66/191 | (P06463_150_155 in P17980 dataset) 20/54 | (P06463_148_152 in P30153 dataset) 47/218 | (P06463_123_128 in P33993 dataset) 399/171 | (P06463_156_158 in Q12959 dataset) 60/91 | (P06463_155_158 in Q14160 dataset) 32/102 | (P06464_31_36 in P85037 dataset) 29/59 | (P06464_32_35 in Q04864 dataset) 36/201 | (P06464_36_39 in Q9NR12 dataset) 35/125 | (P06788_33_36 in P24941 dataset) 66/186 | (P06821_6_11 in P24278 dataset) 11/36 | (P06821_83_89 in P50570 dataset) 24/94 | (P06821_81_85 in P61962 dataset) 46/116 | (P06821_89_92 in Q8N448 dataset) 8/84 | (P06827_95_99 in O00505 dataset) 28/84 | (P06930_13_16 in Q9Y5M8 dataset) 12/59 | (P08393_339_344 in O15379 dataset) 9/78 | (P09992_64_68 in P49368 dataset) 25/139 | (P0C1C6_298_301 in O00505 dataset) 37/84 | (P0C1C7_320_324 in O00505 dataset) 28/84 | (P0C1C7_268_271 in O00505 dataset) 37/84 | (P0C1C7_243_246 in P27348 dataset) 52/161 | (P0C1C7_320_324 in Q02539 dataset) 11/31 | (P0C213_351_353 in Q12959 dataset) 60/91 | (P0C213_350_353 in Q14160 dataset) 32/102 | (P0C739_21_23 in O95817 dataset) 52/129 | (P0C739_23_28 in P51148 dataset) 23/115 | (P0C739_14_18 in Q9P0L0 dataset) 47/180 | (P0C746_111_115 in P17544 dataset) 26/29 | (P0C746_54_58 in P17544 dataset) 26/29 | (P0C746_133_137 in P17544 dataset) 26/29 | (P0CK56_7_10 in P53621 dataset) 17/58 | (P0CK56_101_105 in Q14192 dataset) 25/131 | (P12418_109_113 in Q71U36 dataset) 44/248 | (P13285_67_69 in P37173 dataset) 10/22 | (P16717_146_148 in Q12959 dataset) 60/91 | (P16717_146_148 in Q14160 dataset) 26/102 | (P17382_46_49 in P42574 dataset) 19/50 | (P17386_147_149 in Q12959 dataset) 60/91 | (P21605_93_100 in O60506 dataset) 21/158 | (P21698_139_144 in P48643 dataset) 14/111 | (P21698_79_85 in P50991 dataset) 15/103 | (P21698_134_138 in P78371 dataset) 117/110 | (P21735_156_158 in Q12959 dataset) 60/91 | (P21735_155_158 in Q14160 dataset) 32/102 | (P24772_4_6 in O14920 dataset) 9/88 | (P24835_156_158 in Q12959 dataset) 60/91 | (P24835_155_158 in Q14160 dataset) 32/102 | (P26555_147_149 in Q12959 dataset) 60/91 | (P27228_147_149 in Q12959 dataset) 60/91 | (P27228_146_149 in Q14160 dataset) 32/102 | (P30119_4_8 in P62191 dataset) 11/82 | (P31345_736_739 in O60684 dataset) 17/55 | (P31345_736_739 in P52294 dataset) 23/82 | (P36780_277_282 in Q01130 dataset) 50/59 | (P36780_203_207 in Q9NRL2 dataset) 10/30 | (P50804_156_158 in Q12959 dataset) 60/91 | (P50804_155_158 in Q14160 dataset) 32/102 | (Q01220_101_106 in P42336 dataset) 8/50 | (Q01220_21_25 in P42336 dataset) 7/50 | (Q01220_36_41 in P68104 dataset) 17/179 | (Q04360_185_192 in P38159 dataset) 7/86 | (Q05127_201_204 in P55795 dataset) 11/88 | (Q05322_233_235 in P68431 dataset) 39/84 | (Q0A2H0_228_230 in P78352 dataset) 14/25 | (Q194T2_189_192 in Q92569 dataset) 15/132 | (Q1K9H2_4_8 in O60684 dataset) 18/55 | (Q1K9H2_20_25 in O60812 dataset) 12/28 | (Q1K9H5_232_237 in P33991 dataset) 34/65 | (Q2MG95_260_264 in P20618 dataset) 30/94 | (Q2MG95_260_262 in P51148 dataset) 48/115 | (Q2PJP0_201_204 in P78352 dataset) 11/25 | (Q2PJP0_223_225 in P78352 dataset) 14/25 | (Q2PJP0_223_225 in Q14160 dataset) 52/102 | (Q2PJP1_22_27 in P15822 dataset) 15/56 | (Q2Q067_176_183 in P31249 dataset) 23/5 | (Q5EP28_20_25 in O60812 dataset) 12/28 | (Q67296_191_193 in Q8N448 dataset) 26/84 | (Q67296_752_755 in Q92997 dataset) 61/185 | (Q67296_530_533 in Q92997 dataset) 52/185 | (Q67296_736_739 in Q92997 dataset) 61/185 | (Q69117_21_28 in P23458 dataset) 5/40 | (Q69117_8_12 in P41252 dataset) 27/37 | (Q69117_25_31 in P63151 dataset) 18/55 | (Q6DP93_223_225 in P78352 dataset) 14/25 | (Q6VGS8_126_128 in Q12959 dataset) 60/91 | (Q76S40_201_204 in P63151 dataset) 22/55 | (Q77M19_192_195 in P50991 dataset) 24/103 | (Q77M19_80_84 in P52292 dataset) 57/125 | (Q77M19_229_233 in P52292 dataset) 24/125 | (Q77M19_59_63 in P53621 dataset) 13/58 | (Q77M19_132_139 in P53621 dataset) 7/58 | (Q77M19_229_237 in P62491 dataset) 6/52 | (Q77M19_85_88 in Q09028 dataset) 23/111 | (Q77M19_98_101 in Q16531 dataset) 28/104 | (Q77M19_148_150 in Q9P0L0 dataset) 31/180 | (Q77M19_93_97 in Q9P0L0 dataset) 26/180 | (Q77M19_75_79 in Q9Y230 dataset) 8/134 | (Q77Q36_8_16 in P11802 dataset) 11/107 | (Q77UU1_237_241 in Q07955 dataset) 12/68 | (Q8AZK7_456_460 in O60506 dataset) 27/158 | (Q8AZK7_199_205 in O60506 dataset) 15/158 | (Q8AZK7_504_506 in O60506 dataset) 50/158 | (Q8AZK7_60_64 in O60506 dataset) 27/158 | (Q8AZK7_343_350 in O60506 dataset) 36/158 | (Q8AZK7_192_196 in O60506 dataset) 27/158 | (Q8AZK7_211_218 in O60506 dataset) 36/158 | (Q8AZK7_13_20 in O60506 dataset) 36/158 | (Q8AZK7_390_394 in O60506 dataset) 27/158 | (Q8AZK7_409_416 in O60506 dataset) 36/158 | (Q8AZK7_277_284 in O60506 dataset) 36/158 | (Q8AZK7_126_130 in O60506 dataset) 27/158 | (Q8AZK7_258_262 in O60506 dataset) 27/158 | (Q8AZK7_145_152 in O60506 dataset) 36/158 | (Q8AZK7_324_328 in O60506 dataset) 27/158 | (Q8AZK7_79_86 in O60506 dataset) 36/158 | (Q8AZK7_133_139 in O60506 dataset) 15/158 | (Q8AZK7_265_271 in O60506 dataset) 15/158 | (Q8AZK7_67_73 in O60506 dataset) 15/158 | (Q8AZK7_397_403 in O60506 dataset) 15/158 | (Q8AZK7_331_337 in O60506 dataset) 15/158 | (Q8AZK7_265_271 in O75533 dataset) 20/67 | (Q8AZK7_199_205 in O75533 dataset) 20/67 | (Q8AZK7_397_403 in O75533 dataset) 20/67 | (Q8AZK7_67_73 in O75533 dataset) 20/67 | (Q8AZK7_331_337 in O75533 dataset) 20/67 | (Q8AZK7_133_139 in O75533 dataset) 20/67 | (Q8AZK7_10_14 in P47897 dataset) 55/98 | (Q8AZK7_145_149 in P47897 dataset) 55/98 | (Q8AZK7_79_83 in P47897 dataset) 55/98 | (Q8AZK7_343_347 in P47897 dataset) 55/98 | (Q8AZK7_142_146 in P47897 dataset) 55/98 | (Q8AZK7_13_17 in P47897 dataset) 55/98 | (Q8AZK7_208_212 in P47897 dataset) 55/98 | (Q8AZK7_211_215 in P47897 dataset) 55/98 | (Q8AZK7_340_344 in P47897 dataset) 55/98 | (Q8AZK7_277_281 in P47897 dataset) 55/98 | (Q8AZK7_76_80 in P47897 dataset) 55/98 | (Q8AZK7_274_278 in P47897 dataset) 55/98 | (Q8AZK7_409_413 in P47897 dataset) 55/98 | (Q8AZK7_406_410 in P47897 dataset) 55/98 | (Q8AZK7_346_350 in P61978 dataset) 71/137 | (Q8AZK7_412_416 in P61978 dataset) 71/137 | (Q8AZK7_280_284 in P61978 dataset) 71/137 | (Q8AZK7_82_86 in P61978 dataset) 71/137 | (Q8AZK7_148_152 in P61978 dataset) 71/137 | (Q8AZK7_16_20 in P61978 dataset) 71/137 | (Q8AZK7_214_218 in P61978 dataset) 71/137 | (Q8AZK7_38_42 in Q8N1F7 dataset) 17/46 | (Q8AZK7_104_108 in Q8N1F7 dataset) 17/46 | (Q8AZK7_434_438 in Q8N1F7 dataset) 17/46 | (Q8AZK7_170_174 in Q8N1F7 dataset) 17/46 | (Q8AZK7_368_372 in Q8N1F7 dataset) 17/46 | (Q8AZK7_236_240 in Q8N1F7 dataset) 17/46 | (Q8AZK7_302_306 in Q8N1F7 dataset) 17/46 | (Q8AZK7_354_359 in Q9Y265 dataset) 18/127 | (Q8AZK7_354_360 in Q9Y265 dataset) 18/127 | (Q8AZK7_156_162 in Q9Y265 dataset) 18/127 | (Q8AZK7_24_30 in Q9Y265 dataset) 18/127 | (Q8AZK7_420_425 in Q9Y265 dataset) 18/127 | (Q8AZK7_288_294 in Q9Y265 dataset) 18/127 | (Q8AZK7_222_227 in Q9Y265 dataset) 18/127 | (Q8AZK7_90_95 in Q9Y265 dataset) 18/127 | (Q8AZK7_24_29 in Q9Y265 dataset) 18/127 | (Q8AZK7_156_161 in Q9Y265 dataset) 18/127 | (Q8AZK7_288_293 in Q9Y265 dataset) 18/127 | (Q8AZK7_222_228 in Q9Y265 dataset) 18/127 | (Q8AZK7_420_426 in Q9Y265 dataset) 18/127 | (Q8AZK7_90_96 in Q9Y265 dataset) 18/127 | (Q997F2_320_324 in Q02539 dataset) 11/31 | (Q99AU3_220_223 in Q00839 dataset) 56/303 | (Q99AU3_36_41 in Q00839 dataset) 55/303 | (Q9DGW5_296_301 in P15336 dataset) 20/202 | (Q9DGW5_90_94 in P17544 dataset) 26/29 | (Q9QPN3_72_77 in P12931 dataset) 173/199 | (Q9WMB5_464_467 in P30153 dataset) 58/218 | (Q9WMX2_2273_2276 in P50991 dataset) 24/103 | (Q9WPI5_2_6 in P07910 dataset) 31/99 | (Q9WPI5_2_6 in P11940 dataset) 39/107 | (U5TQE9_57_63 in Q9Y572 dataset) 15/95
## 
## top-1 domain(s): (B4URF7_736_739) IPR016024 | (B4URF7_736_739) IPR011989 | (B4URF7_188_194) IPR011989 | (B4URF7_188_194) IPR027417 | (B4URF7_188_194) IPR016024 | (B4URF7_69_75) IPR011989 | (B4URF7_69_75) IPR027417 | (B4URF7_69_75) IPR016024 | (B4URF7_737_739) IPR027417 | (C5E522_20_25) IPR000504 | (C5E522_443_445) IPR011991 | (C5E526_205_214) IPR016024 | (C5E527_100_105) IPR011009 | (C5E527_100_105) IPR000719 | (C5E527_657_662) IPR011009 | (C5E527_657_662) IPR000719 | (C5E527_158_163) IPR011009 | (C5E527_158_163) IPR000719 | (C5E527_737_739) IPR027417 | (D1LN35_36_41) IPR027417 | (D1LN35_219_223) IPR027417 | (D1LN35_55_60) IPR027417 | (E5LBT9_339_342) IPR016024 | (E5LBT9_339_342) IPR011989 | (E5LBT9_339_341) IPR016024 | (E5LBT9_339_341) IPR011989 | (E5LBT9_339_342) IPR011989 | (E5LBT9_339_342) IPR016024 | (F5HE15_93_96) IPR016024 | (F5HE15_93_96) IPR011989 | (F5HE15_93_96) IPR011989 | (F5HE15_93_96) IPR016024 | (F5HFG5_268_272) IPR027417 | (I6T1Z2_75_77) IPR012336 | (I6T1Z2_35_38) IPR027417 | (I6TAH8_20_25) IPR000504 | (I6TAH8_95_99) IPR011991 | (K7Y1A2_207_212) IPR000719 | (K7Y1A2_207_212) IPR011009 | (O39474_392_396) IPR011990 | (O39474_392_396) IPR013026 | (O40939_383_387) IPR013083 | (O40939_383_387) IPR027417 | (O40939_94_98) IPR013083 | (O40939_94_98) IPR027417 | (O40939_246_250) IPR013083 | (O40939_246_250) IPR027417 | (O40939_296_299) IPR011990 | (O56264_217_219) IPR011009 | (O92837_203_205) IPR016024 | (O92837_203_205) IPR011989 | (P03070_261_264) IPR016024 | (P03070_261_264) IPR011989 | (P03070_129_131) IPR016024 | (P03070_129_131) IPR011989 | (P03070_259_261) IPR011009 | (P03070_259_261) IPR011990 | (P03070_259_261) IPR000719 | (P03070_129_131) IPR011989 | (P03070_129_131) IPR016024 | (P03070_112_115) IPR000719 | (P03070_112_115) IPR011009 | (P03070_129_131) IPR016024 | (P03070_129_131) IPR011989 | (P03070_129_131) IPR011989 | (P03070_129_131) IPR016024 | (P03070_126_130) IPR011991 | (P03070_656_658) IPR001810 | (P03070_656_658) IPR015943 | (P03070_656_658) IPR021977 | (P03070_656_658) IPR017986 | (P03070_656_658) IPR015943 | (P03070_656_658) IPR017986 | (P03070_656_658) IPR021977 | (P03070_656_658) IPR001810 | (P03087_5_7) IPR016024 | (P03087_5_7) IPR011989 | (P03101_499_501) IPR016024 | (P03101_499_501) IPR011989 | (P03101_484_486) IPR016024 | (P03101_484_486) IPR011989 | (P03101_502_504) IPR016024 | (P03101_502_504) IPR011989 | (P03126_156_158) IPR001478 | (P03126_156_158) IPR008144 | (P03126_156_158) IPR008145 | (P03126_156_158) IPR027417 | (P03126_155_158) IPR001478 | (P03129_33_38) IPR013083 | (P03129_41_45) IPR016024 | (P03129_44_47) IPR011009 | (P03129_44_47) IPR000719 | (P03129_34_36) IPR027417 | (P03177_16_20) IPR027417 | (P03188_272_276) IPR008974 | (P03209_569_572) IPR009072 | (P03211_635_638) IPR011989 | (P03211_635_638) IPR016024 | (P03211_458_461) IPR011989 | (P03211_458_461) IPR016024 | (P03211_458_461) IPR016024 | (P03211_458_461) IPR011989 | (P03220_66_68) IPR001781 | (P03225_124_126) IPR001781 | (P03225_45_49) IPR027417 | (P03243_104_107) IPR013083 | (P03246_153_159) IPR032695 | (P03246_153_156) IPR027417 | (P03246_153_156) IPR005225 | (P03246_158_162) IPR027417 | (P03255_285_288) IPR013083 | (P03346_94_97) IPR000210 | (P03346_142_145) IPR000210 | (P03366_107_112) IPR027925 | (P03366_107_112) IPR033762 | (P03366_107_112) IPR012340 | (P03366_107_112) IPR027417 | (P03366_107_112) IPR001208 | (P03366_658_661) IPR011991 | (P03366_217_220) IPR011991 | (P03427_737_739) IPR011989 | (P03427_737_739) IPR016024 | (P03427_30_33) IPR027409 | (P03427_30_33) IPR027413 | (P03427_187_191) IPR027417 | (P03428_188_193) IPR000719 | (P03428_188_193) IPR011009 | (P03428_188_194) IPR027417 | (P03428_188_194) IPR016024 | (P03428_188_194) IPR011989 | (P03428_69_75) IPR027417 | (P03428_69_75) IPR016024 | (P03428_69_75) IPR011989 | (P03428_191_193) IPR013083 | (P03428_191_193) IPR001478 | (P03428_752_755) IPR001478 | (P03428_736_739) IPR001478 | (P03428_530_533) IPR029071 | (P03428_530_533) IPR001478 | (P03428_530_533) IPR011991 | (P03431_205_214) IPR016024 | (P03431_751_753) IPR001478 | (P03431_751_753) IPR013083 | (P03431_669_672) IPR001478 | (P03433_349_352) IPR016024 | (P03433_349_352) IPR011989 | (P03466_20_25) IPR000504 | (P03495_70_75) IPR014720 | (P03495_200_204) IPR011009 | (P03495_200_204) IPR000719 | (P03495_72_77) IPR011009 | (P03495_72_77) IPR000719 | (P03495_217_223) IPR009072 | (P03495_216_220) IPR011991 | (P03495_219_221) IPR001478 | (P03496_71_74) IPR023410 | (P03496_217_219) IPR011009 | (P03496_72_77) IPR000719 | (P03496_72_77) IPR011009 | (P03496_216_220) IPR011991 | (P03496_83_87) IPR011990 | (P03496_83_87) IPR013026 | (P04012_493_495) IPR016024 | (P04012_493_495) IPR011989 | (P04015_239_242) IPR001487 | (P04015_238_242) IPR000504 | (P04296_1170_1172) IPR001357 | (P04296_118_121) IPR027417 | (P04413_145_148) IPR027417 | (P04413_199_203) IPR027409 | (P04413_199_203) IPR027413 | (P04487_54_56) IPR027417 | (P04487_7_9) IPR027417 | (P04591_107_113) IPR000504 | (P05919_51_53) IPR015943 | (P05919_51_53) IPR017986 | (P05919_51_53) IPR001810 | (P05919_51_53) IPR021977 | (P05919_51_53) IPR001810 | (P05919_51_53) IPR015943 | (P05919_51_53) IPR017986 | (P05919_51_53) IPR021977 | (P06427_147_149) IPR001478 | (P06427_146_149) IPR001478 | (P06428_11_15) IPR027417 | (P06429_8_12) IPR013087 | (P06430_31_33) IPR001478 | (P06430_31_33) IPR013083 | (P06460_2_6) IPR013783 | (P06462_7_9) IPR001487 | (P06463_150_155) IPR027417 | (P06463_148_152) IPR011989 | (P06463_148_152) IPR016024 | (P06463_123_128) IPR027417 | (P06463_156_158) IPR001452 | (P06463_156_158) IPR008144 | (P06463_156_158) IPR001478 | (P06463_156_158) IPR008145 | (P06463_156_158) IPR027417 | (P06463_155_158) IPR001478 | (P06464_31_36) IPR000253 | (P06464_31_36) IPR011991 | (P06464_31_36) IPR008984 | (P06464_32_35) IPR013783 | (P06464_36_39) IPR001478 | (P06788_33_36) IPR000719 | (P06788_33_36) IPR011009 | (P06821_6_11) IPR013083 | (P06821_83_89) IPR027417 | (P06821_81_85) IPR015943 | (P06821_81_85) IPR017986 | (P06821_89_92) IPR013083 | (P06827_95_99) IPR011989 | (P06827_95_99) IPR016024 | (P06930_13_16) IPR027417 | (P08393_339_344) IPR023801 | (P09992_64_68) IPR027409 | (P09992_64_68) IPR027413 | (P0C1C6_298_301) IPR011989 | (P0C1C6_298_301) IPR016024 | (P0C1C7_320_324) IPR016024 | (P0C1C7_320_324) IPR011989 | (P0C1C7_268_271) IPR016024 | (P0C1C7_268_271) IPR011989 | (P0C1C7_243_246) IPR023410 | (P0C1C7_320_324) IPR011991 | (P0C213_351_353) IPR001478 | (P0C213_350_353) IPR001478 | (P0C739_21_23) IPR001202 | (P0C739_23_28) IPR027417 | (P0C739_14_18) IPR013783 | (P0C746_111_115) IPR013087 | (P0C746_111_115) IPR004827 | (P0C746_54_58) IPR013087 | (P0C746_54_58) IPR004827 | (P0C746_133_137) IPR013087 | (P0C746_133_137) IPR004827 | (P0CK56_7_10) IPR015943 | (P0CK56_101_105) IPR001781 | (P12418_109_113) IPR003008 | (P12418_109_113) IPR018316 | (P12418_109_113) IPR008280 | (P13285_67_69) IPR011009 | (P16717_146_148) IPR001478 | (P16717_146_148) IPR001478 | (P17382_46_49) IPR015917 | (P17382_46_49) IPR029030 | (P17382_46_49) IPR001309 | (P17382_46_49) IPR002138 | (P17386_147_149) IPR027417 | (P17386_147_149) IPR001478 | (P17386_147_149) IPR008145 | (P17386_147_149) IPR008144 | (P21605_93_100) IPR000504 | (P21698_139_144) IPR027409 | (P21698_139_144) IPR027413 | (P21698_79_85) IPR027413 | (P21698_79_85) IPR027409 | (P21698_134_138) IPR027413 | (P21698_134_138) IPR027409 | (P21735_156_158) IPR001478 | (P21735_155_158) IPR001478 | (P24772_4_6) IPR000719 | (P24772_4_6) IPR011009 | (P24835_156_158) IPR001478 | (P24835_155_158) IPR001478 | (P26555_147_149) IPR027417 | (P27228_147_149) IPR001478 | (P27228_146_149) IPR001478 | (P30119_4_8) IPR027417 | (P31345_736_739) IPR011989 | (P31345_736_739) IPR016024 | (P31345_736_739) IPR016024 | (P31345_736_739) IPR011989 | (P36780_277_282) IPR000504 | (P36780_203_207) IPR013083 | (P50804_156_158) IPR001478 | (P50804_155_158) IPR001478 | (Q01220_101_106) IPR011009 | (Q01220_21_25) IPR011009 | (Q01220_36_41) IPR027417 | (Q04360_185_192) IPR000504 | (Q05127_201_204) IPR000504 | (Q05322_233_235) IPR009072 | (Q0A2H0_228_230) IPR001478 | (Q0A2H0_228_230) IPR001452 | (Q0A2H0_228_230) IPR027417 | (Q194T2_189_192) IPR000980 | (Q1K9H2_4_8) IPR016024 | (Q1K9H2_4_8) IPR011989 | (Q1K9H2_20_25) IPR000504 | (Q1K9H5_232_237) IPR027417 | (Q2MG95_260_264) IPR029055 | (Q2MG95_260_262) IPR027417 | (Q2MG95_260_262) IPR005225 | (Q2PJP0_201_204) IPR001478 | (Q2PJP0_223_225) IPR001452 | (Q2PJP0_223_225) IPR001478 | (Q2PJP0_223_225) IPR027417 | (Q2PJP0_223_225) IPR001478 | (Q2PJP1_22_27) IPR013083 | (Q2Q067_176_183) IPR009057 | (Q5EP28_20_25) IPR000504 | (Q67296_191_193) IPR001478 | (Q67296_191_193) IPR013083 | (Q67296_752_755) IPR001478 | (Q67296_530_533) IPR029071 | (Q67296_530_533) IPR001478 | (Q67296_530_533) IPR011991 | (Q67296_736_739) IPR001478 | (Q69117_21_28) IPR011009 | (Q69117_21_28) IPR000719 | (Q69117_8_12) IPR014729 | (Q69117_25_31) IPR017986 | (Q6DP93_223_225) IPR001478 | (Q6DP93_223_225) IPR001452 | (Q6DP93_223_225) IPR027417 | (Q6VGS8_126_128) IPR027417 | (Q6VGS8_126_128) IPR001452 | (Q76S40_201_204) IPR017986 | (Q76S40_201_204) IPR015943 | (Q77M19_192_195) IPR027409 | (Q77M19_192_195) IPR027413 | (Q77M19_80_84) IPR011989 | (Q77M19_80_84) IPR016024 | (Q77M19_229_233) IPR011989 | (Q77M19_229_233) IPR016024 | (Q77M19_59_63) IPR015943 | (Q77M19_132_139) IPR015943 | (Q77M19_229_237) IPR027417 | (Q77M19_85_88) IPR015943 | (Q77M19_85_88) IPR017986 | (Q77M19_85_88) IPR022052 | (Q77M19_98_101) IPR015943 | (Q77M19_98_101) IPR017986 | (Q77M19_148_150) IPR013783 | (Q77M19_93_97) IPR013783 | (Q77M19_75_79) IPR027417 | (Q77Q36_8_16) IPR000719 | (Q77Q36_8_16) IPR011009 | (Q77UU1_237_241) IPR000504 | (Q8AZK7_456_460) IPR000504 | (Q8AZK7_199_205) IPR000504 | (Q8AZK7_504_506) IPR000504 | (Q8AZK7_60_64) IPR000504 | (Q8AZK7_343_350) IPR000504 | (Q8AZK7_192_196) IPR000504 | (Q8AZK7_211_218) IPR000504 | (Q8AZK7_13_20) IPR000504 | (Q8AZK7_390_394) IPR000504 | (Q8AZK7_409_416) IPR000504 | (Q8AZK7_277_284) IPR000504 | (Q8AZK7_126_130) IPR000504 | (Q8AZK7_258_262) IPR000504 | (Q8AZK7_145_152) IPR000504 | (Q8AZK7_324_328) IPR000504 | (Q8AZK7_79_86) IPR000504 | (Q8AZK7_133_139) IPR000504 | (Q8AZK7_265_271) IPR000504 | (Q8AZK7_67_73) IPR000504 | (Q8AZK7_397_403) IPR000504 | (Q8AZK7_331_337) IPR000504 | (Q8AZK7_265_271) IPR011989 | (Q8AZK7_199_205) IPR011989 | (Q8AZK7_397_403) IPR011989 | (Q8AZK7_67_73) IPR011989 | (Q8AZK7_331_337) IPR011989 | (Q8AZK7_133_139) IPR011989 | (Q8AZK7_10_14) IPR014729 | (Q8AZK7_145_149) IPR014729 | (Q8AZK7_79_83) IPR014729 | (Q8AZK7_343_347) IPR014729 | (Q8AZK7_142_146) IPR014729 | (Q8AZK7_13_17) IPR014729 | (Q8AZK7_208_212) IPR014729 | (Q8AZK7_211_215) IPR014729 | (Q8AZK7_340_344) IPR014729 | (Q8AZK7_277_281) IPR014729 | (Q8AZK7_76_80) IPR014729 | (Q8AZK7_274_278) IPR014729 | (Q8AZK7_409_413) IPR014729 | (Q8AZK7_406_410) IPR014729 | (Q8AZK7_346_350) IPR004087 | (Q8AZK7_346_350) IPR004088 | (Q8AZK7_412_416) IPR004087 | (Q8AZK7_412_416) IPR004088 | (Q8AZK7_280_284) IPR004087 | (Q8AZK7_280_284) IPR004088 | (Q8AZK7_82_86) IPR004087 | (Q8AZK7_82_86) IPR004088 | (Q8AZK7_148_152) IPR004087 | (Q8AZK7_148_152) IPR004088 | (Q8AZK7_16_20) IPR004087 | (Q8AZK7_16_20) IPR004088 | (Q8AZK7_214_218) IPR004087 | (Q8AZK7_214_218) IPR004088 | (Q8AZK7_38_42) IPR011990 | (Q8AZK7_104_108) IPR011990 | (Q8AZK7_434_438) IPR011990 | (Q8AZK7_170_174) IPR011990 | (Q8AZK7_368_372) IPR011990 | (Q8AZK7_236_240) IPR011990 | (Q8AZK7_302_306) IPR011990 | (Q8AZK7_354_359) IPR027417 | (Q8AZK7_354_359) IPR003593 | (Q8AZK7_354_360) IPR027417 | (Q8AZK7_354_360) IPR003593 | (Q8AZK7_156_162) IPR027417 | (Q8AZK7_156_162) IPR003593 | (Q8AZK7_24_30) IPR027417 | (Q8AZK7_24_30) IPR003593 | (Q8AZK7_420_425) IPR027417 | (Q8AZK7_420_425) IPR003593 | (Q8AZK7_288_294) IPR027417 | (Q8AZK7_288_294) IPR003593 | (Q8AZK7_222_227) IPR027417 | (Q8AZK7_222_227) IPR003593 | (Q8AZK7_90_95) IPR027417 | (Q8AZK7_90_95) IPR003593 | (Q8AZK7_24_29) IPR027417 | (Q8AZK7_24_29) IPR003593 | (Q8AZK7_156_161) IPR027417 | (Q8AZK7_156_161) IPR003593 | (Q8AZK7_288_293) IPR027417 | (Q8AZK7_288_293) IPR003593 | (Q8AZK7_222_228) IPR027417 | (Q8AZK7_222_228) IPR003593 | (Q8AZK7_420_426) IPR027417 | (Q8AZK7_420_426) IPR003593 | (Q8AZK7_90_96) IPR027417 | (Q8AZK7_90_96) IPR003593 | (Q997F2_320_324) IPR011991 | (Q99AU3_220_223) IPR027417 | (Q99AU3_36_41) IPR027417 | (Q9DGW5_296_301) IPR004827 | (Q9DGW5_90_94) IPR004827 | (Q9DGW5_90_94) IPR013087 | (Q9QPN3_72_77) IPR001245 | (Q9QPN3_72_77) IPR000980 | (Q9QPN3_72_77) IPR001452 | (Q9QPN3_72_77) IPR000719 | (Q9QPN3_72_77) IPR011009 | (Q9WMB5_464_467) IPR016024 | (Q9WMB5_464_467) IPR011989 | (Q9WMX2_2273_2276) IPR027413 | (Q9WMX2_2273_2276) IPR027409 | (Q9WPI5_2_6) IPR000504 | (Q9WPI5_2_6) IPR000504 | (U5TQE9_57_63) IPR011009 | (U5TQE9_57_63) IPR000719
## 
## domain_support4motif_nq / protein_with_motif_degree: 12/210 | 12/210 | 19/210 | 28/210 | 21/210 | 19/210 | 28/210 | 21/210 | 13/210 | 6/77 | 11/77 | 6/124 | 13/140 | 13/140 | 13/140 | 13/140 | 26/140 | 26/140 | 13/140 | 31/42 | 11/42 | 31/42 | 12/6 | 12/6 | 22/6 | 22/6 | 12/6 | 12/6 | 12/12 | 12/12 | 12/12 | 12/12 | 8/25 | 6/138 | 17/138 | 6/106 | 4/106 | 29/3 | 29/3 | 19/3 | 17/3 | 21/32 | 23/32 | 21/32 | 23/32 | 21/32 | 23/32 | 8/32 | 7/108 | 22/4 | 22/4 | 19/63 | 19/63 | 13/63 | 13/63 | 34/63 | 16/63 | 34/63 | 21/63 | 23/63 | 18/63 | 19/63 | 22/63 | 22/63 | 17/63 | 17/63 | 4/63 | 21/63 | 25/63 | 21/63 | 23/63 | 29/63 | 27/63 | 23/63 | 26/63 | 22/4 | 22/4 | 22/4 | 22/4 | 22/4 | 22/4 | 22/4 | 22/4 | 18/32 | 6/32 | 6/32 | 8/32 | 11/32 | 8/97 | 14/97 | 18/97 | 18/97 | 17/97 | 9/43 | 26/10 | 8/13 | 19/51 | 19/51 | 18/51 | 20/51 | 23/51 | 23/51 | 21/61 | 21/78 | 6/78 | 13/13 | 8/95 | 21/95 | 19/95 | 9/95 | 13/30 | 9/6 | 9/6 | 22/38 | 22/38 | 22/38 | 23/38 | 22/38 | 9/38 | 9/38 | 13/168 | 13/168 | 12/168 | 12/168 | 10/168 | 19/166 | 19/166 | 28/166 | 21/166 | 19/166 | 28/166 | 21/166 | 19/166 | 15/166 | 8/166 | 35/166 | 35/166 | 24/166 | 28/166 | 20/166 | 6/133 | 8/133 | 15/133 | 35/133 | 19/48 | 19/48 | 6/143 | 11/60 | 20/60 | 20/60 | 19/60 | 19/60 | 6/60 | 4/60 | 15/60 | 33/212 | 7/212 | 19/212 | 19/212 | 4/212 | 19/212 | 17/212 | 22/3 | 22/3 | 11/43 | 9/43 | 14/46 | 11/46 | 54/29 | 22/29 | 22/29 | 23/145 | 49/145 | 13/11 | 25/55 | 23/55 | 21/55 | 21/55 | 26/55 | 29/55 | 27/55 | 23/55 | 18/4 | 11/4 | 6/111 | 5/52 | 8/51 | 15/51 | 23/276 | 18/12 | 14/81 | 23/81 | 24/81 | 98/81 | 7/81 | 6/81 | 18/81 | 6/81 | 8/81 | 11/81 | 14/48 | 18/48 | 15/48 | 14/48 | 18/48 | 34/33 | 35/33 | 8/101 | 16/101 | 29/101 | 28/101 | 6/101 | 13/9 | 13/9 | 9/20 | 6/18 | 16/11 | 16/11 | 19/53 | 19/53 | 13/29 | 13/29 | 19/29 | 19/29 | 33/29 | 4/29 | 18/3 | 11/3 | 11/252 | 16/252 | 16/252 | 11/14 | 8/14 | 11/14 | 8/14 | 11/14 | 8/14 | 13/202 | 12/202 | 23/9 | 23/9 | 23/9 | 8/74 | 18/2 | 11/2 | 5/2 | 5/2 | 5/2 | 5/2 | 8/6 | 18/6 | 6/6 | 6/6 | 11/42 | 11/26 | 11/26 | 10/26 | 10/26 | 60/26 | 60/26 | 18/2 | 11/2 | 8/16 | 8/16 | 18/2 | 11/2 | 8/2 | 18/2 | 11/2 | 6/76 | 12/3 | 12/3 | 12/3 | 12/3 | 26/89 | 8/89 | 18/3 | 11/3 | 5/21 | 6/21 | 15/21 | 5/10 | 7/60 | 21/90 | 6/41 | 7/41 | 5/41 | 10/49 | 16/112 | 15/112 | 6/112 | 23/192 | 16/196 | 32/196 | 26/196 | 6/42 | 7/42 | 6/42 | 5/42 | 19/42 | 8/36 | 3/10 | 6/81 | 8/25 | 15/25 | 35/25 | 24/25 | 28/25 | 20/25 | 35/25 | 4/77 | 4/77 | 14/77 | 11/77 | 6/32 | 7/32 | 5/32 | 8/26 | 7/26 | 14/8 | 13/8 | 20/438 | 20/438 | 19/438 | 20/438 | 16/438 | 16/438 | 11/438 | 6/438 | 4/438 | 15/438 | 14/438 | 13/438 | 20/438 | 20/438 | 11/438 | 11/438 | 6/438 | 7/4 | 7/4 | 9/39 | 12/150 | 5/150 | 17/150 | 12/150 | 12/150 | 12/150 | 12/150 | 12/150 | 12/150 | 12/150 | 12/150 | 12/150 | 12/150 | 12/150 | 12/150 | 12/150 | 5/150 | 5/150 | 5/150 | 5/150 | 5/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 14/150 | 14/150 | 14/150 | 14/150 | 14/150 | 14/150 | 14/150 | 14/150 | 14/150 | 14/150 | 14/150 | 14/150 | 14/150 | 14/150 | 20/150 | 20/150 | 20/150 | 20/150 | 20/150 | 20/150 | 20/150 | 20/150 | 20/150 | 20/150 | 20/150 | 20/150 | 20/150 | 20/150 | 8/150 | 8/150 | 8/150 | 8/150 | 8/150 | 8/150 | 8/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 10/150 | 4/34 | 37/24 | 31/24 | 11/13 | 8/13 | 11/13 | 77/5 | 79/5 | 82/5 | 80/5 | 82/5 | 23/13 | 24/13 | 20/61 | 20/61 | 23/61 | 22/61 | 12/2 | 12/2
## 
## all with motif: A0A024R8A9 | A0A024R8L2 | A0AUL9 | A0JLT2 | A3KN83 | A6ND36 | A6NHQ4 | A6NK89 | A7KAX9-2 | A8E1C4 | A8K1F4 | A8MTZ0 | A8MW92 | B2RXH8 | B4DE84 | B4URF7 | B7Z2Y1 | B8PZP3 | C5E522 | C5E526 | C5E527 | C9JA69 | C9JIZ3 | D0UZS0 | D1LN35 | D3DR86 | E5LBT9 | E9PK67 | F5HE15 | F5HFG5 | G3V1X1 | H2A1D8 | I3L2W2 | I6L996 | I6T1Z2 | I6TAH8 | J3KNL6 | K7Y1A2 | O00159 | O00167-2 | O00221 | O00233 | O00254 | O00273 | O00303 | O00358 | O00418 | O00422 | O00444 | O00487 | O00506 | O00560 | O00567 | O00571 | O00716 | O00757 | O14490 | O14503 | O14530 | O14576 | O14579 | O14639 | O14641 | O14645 | O14654 | O14713 | O14744 | O14746 | O14776 | O14862 | O14893 | O14908 | O14920 | O14964 | O14974 | O15027 | O15042 | O15078 | O15160 | O15169 | O15198 | O15209 | O15212 | O15226 | O15259 | O15264 | O15287 | O15294 | O15297 | O15315 | O15371 | O15381 | O15399 | O15455 | O15481 | O15535 | O39474 | O40939 | O43147 | O43148 | O43150 | O43159 | O43166 | O43167 | O43172 | O43189 | O43264 | O43281 | O43296 | O43299 | O43318-2 | O43395 | O43447 | O43463 | O43491 | O43504 | O43524 | O43542 | O43633 | O43684 | O43684-2 | O43741 | O43772 | O43795 | O43815 | O43847 | O43852 | O43896 | O43900 | O43909 | O43918 | O55779 | O56264 | O60216 | O60232 | O60264 | O60271 | O60318 | O60333 | O60333-3 | O60341 | O60346 | O60356 | O60493 | O60502 | O60504 | O60506 | O60563 | O60566 | O60684 | O60716 | O60739 | O60828 | O60884 | O60934 | O75031 | O75044 | O75052-3 | O75128 | O75161 | O75164 | O75190 | O75223 | O75365 | O75376 | O75381 | O75385 | O75398 | O75400 | O75400-2 | O75410 | O75420 | O75449 | O75496 | O75506 | O75525 | O75530 | O75533 | O75554 | O75604 | O75643 | O75665 | O75683 | O75689 | O75694 | O75815 | O75821 | O75832 | O75888 | O75909 | O75914 | O75923 | O75955 | O76031 | O76041 | O76071 | O76075 | O76081 | O92837 | O94760 | O94776 | O94782 | O94817 | O94905 | O94966 | O95067 | O95070 | O95071 | O95157 | O95159 | O95166 | O95182 | O95229 | O95251 | O95299 | O95365 | O95402 | O95405 | O95425 | O95425-2 | O95453 | O95487 | O95503 | O95602 | O95677 | O95684 | O95696 | O95714 | O95758 | O95782-2 | O95816 | O95817 | O95835 | O95848 | O95863 | O95983 | O95995 | O95999 | O96017 | P00519 | P00533 | P01023 | P01100 | P01106 | P01112 | P01116 | P02545 | P02675 | P02751 | P03070 | P03087 | P03101 | P03107 | P03120 | P03126 | P03129 | P03177 | P03179 | P03188 | P03209 | P03211 | P03220 | P03225 | P03243 | P03246 | P03255 | P03255-2 | P03346 | P03366 | P03372 | P03407 | P03427 | P03428 | P03431 | P03433 | P03466 | P03495 | P03496 | P04004 | P04012 | P04015 | P04075 | P04183 | P04296 | P04406 | P04413 | P04487 | P04591 | P04604 | P04626 | P04637 | P04792 | P05023 | P05067 | P05090 | P05412 | P05783 | P05787 | P05919 | P06401 | P06422 | P06423 | P06427 | P06428 | P06429 | P06430 | P06460 | P06462 | P06463 | P06464 | P06748 | P06753-2 | P06756 | P06788 | P06790 | P06821 | P06827 | P06921 | P06930 | P07814 | P07900 | P07949 | P07951-2 | P07954 | P08047 | P08069 | P08238 | P08240 | P08393 | P08473 | P08727 | P09012 | P09430 | P09466 | P09496-2 | P09497-2 | P09619 | P09651 | P09758 | P09874 | P09992 | P0C1C6 | P0C1C7 | P0C213 | P0C739 | P0C746 | P0CG47 | P0CK56 | P0DMU9 | P10275 | P10600 | P10636 | P10636-2 | P10636-6 | P10644 | P10909 | P11021 | P11142 | P11171 | P11274 | P11387 | P11388 | P11586 | P11717 | P11802 | P12004 | P12418 | P12532 | P12814 | P12830 | P12956 | P13010 | P13285 | P13473 | P13569 | P13639 | P13994 | P14136 | P14373 | P14618 | P14866 | P14923 | P15336 | P15408 | P15586 | P15822 | P15880 | P15884 | P15884-3 | P15923 | P15924 | P15927 | P16234 | P16383 | P16401 | P16615 | P16717 | P17028 | P17382 | P17386 | P17542 | P17844 | P17980 | P17987 | P18031 | P18847 | P18848 | P18887 | P19105 | P19320 | P19338 | P19404 | P19474 | P19838 | P20073 | P20226 | P20333 | P20700 | P20719 | P20774 | P20810 | P20930 | P20936 | P20963 | P21281 | P21333 | P21605 | P21698 | P21735 | P21796 | P22459 | P22460 | P22681 | P22736 | P23025 | P23258 | P23396 | P23508 | P23511-2 | P23759 | P23760 | P24278 | P24772 | P24830 | P24835 | P24863 | P25054 | P25205 | P25391 | P25445 | P25490 | P25685 | P25786 | P25789 | P25963 | P26367 | P26368 | P26373 | P26555 | P27228 | P27658 | P27695 | P27708 | P27958 | P27986 | P28289 | P28290 | P28340 | P28749 | P29317 | P29353 | P29401 | P29590 | P29597 | P29692 | P30119 | P30153 | P30304 | P30305 | P30480 | P30530 | P30622-2 | P31040 | P31150 | P31273 | P31345 | P31629 | P31689 | P31749 | P31930 | P31943 | P31946 | P31947 | P31947-2 | P31948 | P32121 | P32314 | P33240 | P33402 | P33991 | P33992 | P33993 | P33993-2 | P34969 | P35221 | P35222 | P35232 | P35250 | P35268 | P35372 | P35579 | P35610 | P35638 | P35968 | P35998 | P36578 | P36778 | P36780 | P37198 | P37840 | P38398 | P38432 | P38919 | P38935 | P38936 | P39019 | P39023 | P39687 | P40337 | P40429 | P40692 | P40818 | P41182 | P41208 | P41218 | P41226 | P42226 | P42285 | P42345 | P42566 | P42568 | P42574 | P42684 | P42695 | P42858 | P43115-12 | P43243 | P43246 | P43351 | P43686 | P43699 | P45973 | P45983 | P45984 | P46013 | P46063 | P46087 | P46100 | P46108 | P46108-2 | P46109 | P46527 | P46531 | P46778 | P46779 | P46781 | P46782 | P46937 | P46940 | P47756 | P47756-2 | P47914 | P47928 | P48050 | P48380 | P48643 | P49321 | P49321-2 | P49327 | P49336 | P49368 | P49407 | P49427 | P49450 | P49593 | P49662 | P49711 | P49736 | P49755 | P49757 | P49759-3 | P49761 | P49790 | P49792 | P49815 | P49848 | P49916 | P49918 | P50570 | P50616 | P50804 | P50991 | P51116 | P51532 | P51610 | P51617 | P51665 | P51798 | P51858 | P52272 | P52732 | P52735 | P52907 | P53350 | P53355 | P53396 | P53567 | P53582 | P53778 | P53804 | P53990 | P53990-2 | P54253 | P54646 | P54725 | P54819 | P55036 | P55072 | P55209 | P55212 | P55884 | P56945 | P57060 | P57078 | P57682 | P57796-2 | P58012 | P58340 | P58753 | P60520 | P60842 | P60953 | P61244 | P61254 | P61289 | P61326 | P61353 | P61421 | P61962 | P61978 | P62136 | P62241 | P62269 | P62277 | P62280 | P62424 | P62491 | P62633 | P62714 | P62750 | P62753 | P62805 | P62829 | P62841 | P62910 | P62917 | P62995 | P63000 | P63010 | P63104 | P63167 | P63208 | P63261 | P63279 | P67809 | P67936 | P68400 | P68431 | P69284 | P78318 | P78337 | P78345 | P78347 | P78423 | P78424 | P78527 | P78536 | P78563-4 | P82914 | P82930 | P82933 | P83731 | P83881 | P83916 | P84243 | P85037 | P98082 | P98170 | P98175 | P98179 | Q00005 | Q00005-7 | Q00335 | Q00341 | Q00403 | Q00444 | Q00534 | Q00535 | Q00537 | Q00610-2 | Q00613 | Q00653 | Q00839 | Q00987 | Q01082 | Q01167 | Q01196 | Q01201 | Q01220 | Q01780 | Q01804 | Q01814 | Q01844 | Q01968 | Q01973 | Q02078 | Q02543 | Q02750 | Q03252 | Q03518 | Q03701 | Q04206 | Q04206-2 | Q04323 | Q04360 | Q04637 | Q04864 | Q05086-3 | Q05127 | Q05193 | Q05322 | Q05397 | Q05516 | Q05639 | Q05996 | Q05BL1 | Q06190 | Q06547 | Q06609 | Q07065 | Q07666 | Q07817 | Q07889 | Q07890 | Q08050 | Q08117 | Q08117-2 | Q08170 | Q08379 | Q08945 | Q08999 | Q09161 | Q09472 | Q09666 | Q0A2H0 | Q0HD54 | Q0VD86 | Q12772 | Q12797 | Q12830 | Q12834 | Q12841 | Q12873 | Q12879 | Q12888 | Q12905 | Q12906 | Q12923 | Q12933 | Q12952 | Q12955 | Q12968 | Q12996 | Q13023 | Q13033 | Q13043 | Q13064 | Q13085 | Q13094 | Q13098 | Q13112 | Q13127 | Q13136 | Q13137 | Q13148 | Q13153 | Q13162 | Q13163 | Q13177 | Q13185 | Q13191 | Q13200 | Q13224 | Q13233 | Q13255 | Q13263 | Q13315 | Q13330 | Q13332 | Q13362 | Q13371 | Q13395 | Q13404 | Q13415 | Q13418 | Q13428 | Q13435 | Q13444 | Q13451 | Q13469 | Q13472 | Q13480 | Q13501 | Q13509 | Q13526 | Q13541 | Q13547 | Q13555 | Q13561 | Q13573 | Q13595 | Q13601 | Q13615 | Q13619 | Q13620 | Q13625 | Q13627 | Q13761 | Q13796 | Q13813 | Q13823 | Q13905 | Q14008 | Q14011 | Q14103 | Q14103-4 | Q14118 | Q14126 | Q14141 | Q14149 | Q14160 | Q14164 | Q14185 | Q14186 | Q14190 | Q14191 | Q14203 | Q14204 | Q14247 | Q14257 | Q14289 | Q14494 | Q14498 | Q14500 | Q14511 | Q14524 | Q14566 | Q14573 | Q14653 | Q14669 | Q14674 | Q14676 | Q14683 | Q14684 | Q14738 | Q14781 | Q14807 | Q14839 | Q14849 | Q14966 | Q14BN4 | Q14D33 | Q15006 | Q15008 | Q15013 | Q15018 | Q15022 | Q15036 | Q15051 | Q15052 | Q15075 | Q15149 | Q15154 | Q15257 | Q15287 | Q15303 | Q15311 | Q15397 | Q15424 | Q15428 | Q15468 | Q15477 | Q15554 | Q15562 | Q15637 | Q15645 | Q15648 | Q15652 | Q15653 | Q15672 | Q15691 | Q15696 | Q15714 | Q15742 | Q15796 | Q15797 | Q15811 | Q15834 | Q15906 | Q15910 | Q16186 | Q16254 | Q16401 | Q16512 | Q16513 | Q16514 | Q16531 | Q16543 | Q16594 | Q16637 | Q16643 | Q16644 | Q16650 | Q16659 | Q16665 | Q16666-2 | Q16667 | Q16695 | Q16795 | Q16891 | Q16891-2 | Q194T2 | Q1K9H2 | Q1K9H5 | Q1KMD3 | Q20MH3 | Q29RF7 | Q2HR82 | Q2KJY2 | Q2MG95 | Q2MV58 | Q2PJP0 | Q2PJP1 | Q2Q067 | Q3SYB3 | Q3V6T2 | Q49A26-4 | Q49AN0 | Q4FZB7 | Q4G0J3 | Q4G0X9 | Q4G0X9-5 | Q4LE39 | Q4TT70 | Q53EL6 | Q53FC7 | Q53GS9 | Q53HY9 | Q53TQ3 | Q53TS8 | Q567U6 | Q59FP4 | Q5D862 | Q5EP28 | Q5EP37 | Q5FWF4 | Q5HYA8 | Q5HYN5 | Q5HYW2 | Q5JR59-3 | Q5JTC6 | Q5JTH9 | Q5JU00 | Q5JVS0 | Q5PSV4 | Q5S007 | Q5SQQ9-2 | Q5SW79 | Q5T3F8 | Q5T4F4 | Q5T4S7 | Q5T5P2-6 | Q5T7B8 | Q5T7W7 | Q5T8D3-2 | Q5TA45 | Q5TAQ9 | Q5TBB1 | Q5TCZ1 | Q5U5Q3 | Q5UIP0 | Q5VSL9 | Q5VU43 | Q5VWN6 | Q5VY09 | Q5VZL5 | Q60I27 | Q63ZY3 | Q66GS9 | Q66K89 | Q66LE6 | Q67296 | Q676U5 | Q68CZ1 | Q68CZ2 | Q68D86 | Q69117 | Q69YQ0 | Q6A162 | Q6BDS2 | Q6DD88 | Q6DN90 | Q6DP93 | Q6FHJ7 | Q6IBT1 | Q6IBW4 | Q6IN84 | Q6IN85 | Q6IQ16 | Q6IQ23 | Q6K0P9 | Q6NSI4 | Q6NUS6 | Q6NXE6 | Q6NYC1 | Q6NZI2 | Q6P1J9 | Q6P2E9 | Q6P2Q9 | Q6P597 | Q6P5Z2 | Q6P9B9 | Q6PCD5 | Q6PI77 | Q6PIV2 | Q6PK04 | Q6PKG0 | Q6PNE5 | Q6SPF0 | Q6T4R5-2 | Q6UN15 | Q6UVJ0 | Q6VAB6 | Q6VGS8 | Q6W2J9 | Q6WCQ1 | Q6XUX3 | Q6ZNH5 | Q6ZQN5 | Q6ZRS2 | Q6ZRY4 | Q6ZS30 | Q6ZTQ3 | Q6ZU15 | Q6ZU80 | Q6ZUT1 | Q6ZUT1-2 | Q6ZW61 | Q71F23 | Q76N32 | Q76S40 | Q77M19 | Q77Q36 | Q77UU1 | Q77UV9 | Q7KZF4 | Q7KZI7 | Q7KZS0 | Q7L0Y3 | Q7L273 | Q7L2E3 | Q7L4I2 | Q7L590 | Q7L7X3 | Q7L804 | Q7Z2E3 | Q7Z2T5 | Q7Z333 | Q7Z3C6 | Q7Z401 | Q7Z4S6 | Q7Z5V6-2 | Q7Z614 | Q7Z628 | Q7Z6G8 | Q7Z6Z7 | Q7Z7A1 | Q7Z7C8-2 | Q7Z7J5 | Q86SE5-3 | Q86SE9 | Q86TG7 | Q86TH3 | Q86UK5 | Q86UR1-2 | Q86V38 | Q86V81 | Q86VF2-5 | Q86VM9 | Q86W54 | Q86W92 | Q86WG5 | Q86WH2 | Q86X10 | Q86X95 | Q86Y26 | Q86YD7 | Q86YS7 | Q86Z20 | Q88489 | Q8AZK7 | Q8HWS3 | Q8IUD2 | Q8IUH5 | Q8IVD9 | Q8IVT2 | Q8IVT4 | Q8IVT5 | Q8IW19 | Q8IW35 | Q8IWL3 | Q8IWZ6 | Q8IWZ8 | Q8IX03 | Q8IX06 | Q8IX15-3 | Q8IXK2 | Q8IY81 | Q8IY92 | Q8IYA8 | Q8IZD4 | Q8IZD9 | Q8IZP0-5 | Q8IZQ1 | Q8IZQ1-2 | Q8IZU0 | Q8IZY2 | Q8N157 | Q8N163 | Q8N2Q7 | Q8N2W9 | Q8N3C0 | Q8N3P4 | Q8N3V7-2 | Q8N402 | Q8N448 | Q8N4C6 | Q8N4V1 | Q8N531 | Q8N5L8 | Q8N6L0 | Q8N6M6 | Q8N6Y0 | Q8N720 | Q8N7W2-2 | Q8N8A6 | Q8N8B7-2 | Q8N8S7 | Q8N960 | Q8N9N5 | Q8N9N8 | Q8N9T8 | Q8NA19-2 | Q8NB46 | Q8NBS9 | Q8NC44 | Q8NC51 | Q8NCD3 | Q8NCN4 | Q8ND90 | Q8NE01 | Q8NE31 | Q8NEM7-2 | Q8NEU8 | Q8NFF5-2 | Q8NFJ9 | Q8NFW5 | Q8NG31 | Q8NG31-2 | Q8NHU0 | Q8NHV4 | Q8NI38 | Q8TAD8 | Q8TAE8 | Q8TAM1 | Q8TAU3 | Q8TB24 | Q8TBB1 | Q8TBE0 | Q8TBK6 | Q8TBZ8 | Q8TC59 | Q8TCG1 | Q8TCT1 | Q8TCU6 | Q8TD10 | Q8TDM6-4 | Q8TDR0 | Q8TDY4 | Q8TES7 | Q8TEW0 | Q8TEX9 | Q8TF72 | Q8TF76 | Q8WU02 | Q8WUA4 | Q8WUH2 | Q8WV41 | Q8WV92 | Q8WVD3 | Q8WVM7 | Q8WVM8 | Q8WW22 | Q8WWM7 | Q8WWW0 | Q8WWW0-2 | Q8WWY3 | Q8WXI2 | Q8WXI9 | Q8WXX7 | Q8WY64 | Q8WYP5 | Q8WZ74 | Q90VU7 | Q92466 | Q92508 | Q92530 | Q92542 | Q92597 | Q92609 | Q92614 | Q92630 | Q92636 | Q92685 | Q92688 | Q92731 | Q92769 | Q92793 | Q92830 | Q92833 | Q92844 | Q92845 | Q92900 | Q92918 | Q92934 | Q92956 | Q92974 | Q92988 | Q92993 | Q92997 | Q93009 | Q93050 | Q93074 | Q969G5 | Q969M3 | Q969R5 | Q969U6 | Q96A65 | Q96AE4 | Q96AQ6 | Q96AY2 | Q96BD8 | Q96BE0 | Q96BK5 | Q96BT7 | Q96BY7 | Q96C10 | Q96C36 | Q96C86 | Q96CA5 | Q96CS3 | Q96CV9 | Q96DE5 | Q96DT6 | Q96E09 | Q96EB6 | Q96ED9-2 | Q96EG3 | Q96EK4 | Q96EZ8 | Q96F86 | Q96FJ2 | Q96FW1 | Q96GG9 | Q96GQ7 | Q96GX1 | Q96H12 | Q96H55 | Q96HA1-2 | Q96HS1 | Q96HS1-2 | Q96I25 | Q96IY1 | Q96IZ5 | Q96IZ7 | Q96JM7 | Q96JM7-2 | Q96JN0 | Q96JY6 | Q96K76 | Q96K78 | Q96KM6 | Q96KP6 | Q96L14 | Q96L91 | Q96LM6 | Q96MH2 | Q96MU7 | Q96MY7 | Q96N21 | Q96PE1 | Q96PQ6 | Q96PU4 | Q96PV4 | Q96Q15 | Q96Q35 | Q96Q45 | Q96Q77 | Q96QF0 | Q96QF0-2 | Q96RD7 | Q96RK4 | Q96RL1 | Q96RR4 | Q96RS6 | Q96S96 | Q96SB4 | Q96SI9 | Q96ST3 | Q96SY0 | Q96T60 | Q96TC7 | Q99081-3 | Q99436 | Q99459 | Q99471 | Q99519 | Q99558 | Q99569 | Q99575 | Q99614 | Q99615 | Q99623 | Q99683 | Q99689 | Q99697 | Q99708 | Q99728 | Q99729-3 | Q99741 | Q99750 | Q99759 | Q997F2 | Q99832 | Q99853 | Q99880 | Q99958 | Q99AU3 | Q9BPU9 | Q9BPX5 | Q9BPZ3 | Q9BPZ7 | Q9BQ15 | Q9BQ24 | Q9BQ39 | Q9BQ95 | Q9BQD3 | Q9BQG0 | Q9BQS8 | Q9BRD0 | Q9BRK4 | Q9BRP4 | Q9BSI4 | Q9BSW7 | Q9BTC0 | Q9BTE3 | Q9BTE3-2 | Q9BTE6-2 | Q9BU70 | Q9BU76 | Q9BUF5 | Q9BUH8 | Q9BUI4 | Q9BUJ2-2 | Q9BUR4 | Q9BUZ4 | Q9BVG8 | Q9BVG8-5 | Q9BVI0 | Q9BVW5 | Q9BW61 | Q9BW85 | Q9BWW9 | Q9BX10-4 | Q9BXC9 | Q9BXF6 | Q9BXK1 | Q9BXL5 | Q9BXL8 | Q9BXP5 | Q9BXY8 | Q9BY76 | Q9BYB0 | Q9BYE7 | Q9BYG5 | Q9BYX2 | Q9BYX4 | Q9BZ95 | Q9BZD4 | Q9BZE4 | Q9BZF3 | Q9BZL4 | Q9C005 | Q9C009 | Q9C0C7 | Q9C0F1 | Q9C0H9 | Q9DGW5 | Q9E7P0 | Q9GZQ8 | Q9GZR2 | Q9GZS1 | Q9GZU8 | Q9GZV5 | Q9H074 | Q9H081 | Q9H0A0 | Q9H0D6 | Q9H0I2 | Q9H0K1 | Q9H0R8 | Q9H0S4 | Q9H160 | Q9H165 | Q9H171 | Q9H1J1 | Q9H1P3 | Q9H1R2-3 | Q9H204 | Q9H3D4 | Q9H3R2 | Q9H3S7 | Q9H444 | Q9H492 | Q9H4B7 | Q9H4L5 | Q9H4X1-2 | Q9H501 | Q9H5I1 | Q9H5V8 | Q9H6F5 | Q9H7D0 | Q9H7D7 | Q9H7E9 | Q9H814 | Q9H832 | Q9H9B4 | Q9H9D4 | Q9H9F9 | Q9H9G7 | Q9H9H4 | Q9H9J4 | Q9H9L3 | Q9HAK2 | Q9HAN9 | Q9HAU0-4 | Q9HAV4 | Q9HAW0 | Q9HBH7 | Q9HBI0 | Q9HC52 | Q9HCJ0 | Q9HCK4 | Q9HCK5 | Q9HCK8 | Q9HCM9-2 | Q9HCU9 | Q9HD26 | Q9HD67 | Q9NP66 | Q9NP77 | Q9NPB6-2 | Q9NPD3 | Q9NPF5 | Q9NPI1 | Q9NPI6 | Q9NPJ1 | Q9NQ76 | Q9NQC3 | Q9NQP4 | Q9NQS7 | Q9NQT8 | Q9NQW6 | Q9NQW7 | Q9NQX5 | Q9NQZ6 | Q9NR30 | Q9NR45 | Q9NRD5 | Q9NRE2 | Q9NRI5-2 | Q9NRJ4 | Q9NRL3 | Q9NS75 | Q9NSK0 | Q9NTI5 | Q9NUC0 | Q9NV70 | Q9NVH1 | Q9NVK5 | Q9NVV9 | Q9NVW2 | Q9NW13 | Q9NWA0 | Q9NWQ8 | Q9NWT6 | Q9NX18 | Q9NX40 | Q9NX58 | Q9NX63 | Q9NX70 | Q9NXB0 | Q9NXC5 | Q9NXV6 | Q9NY61 | Q9NYA1 | Q9NYA3 | Q9NYB0 | Q9NYB5 | Q9NYF8 | Q9NYL2 | Q9NYL9 | Q9NYU2 | Q9NZ43 | Q9NZB2 | Q9NZC7 | Q9NZI8 | Q9NZJ5 | Q9NZM1-6 | Q9NZM4 | Q9NZV5 | Q9NZZ3 | Q9P021 | Q9P035 | Q9P0K1 | Q9P0K1-3 | Q9P0K8 | Q9P0T4 | Q9P0V3 | Q9P127 | Q9P260 | Q9P275 | Q9P289 | Q9P2B4 | Q9P2D1 | Q9P2E9 | Q9P2H0 | Q9P2X3 | Q9QPN3 | Q9UBD5 | Q9UBF8 | Q9UBK2 | Q9UBL3 | Q9UBN6 | Q9UBN7 | Q9UBP5 | Q9UBS4 | Q9UBU8 | Q9UBU9 | Q9UBX2 | Q9UDY2 | Q9UDY4 | Q9UER7 | Q9UEY8 | Q9UGL9 | Q9UH92-3 | Q9UH99 | Q9UHB6 | Q9UHD2 | Q9UHL9 | Q9UHL9-2 | Q9UHR4 | Q9UHR5 | Q9UHV9 | Q9UHX1 | Q9UI10-3 | Q9UIG0 | Q9UII2 | Q9UJ70-2 | Q9UJU2 | Q9UJU6 | Q9UJV9 | Q9UJW9 | Q9UJX4 | Q9UJZ1 | Q9UK58 | Q9UKA9 | Q9UKD1 | Q9UKE5 | Q9UKF6 | Q9UKG1 | Q9UKL0 | Q9UKV0-3 | Q9UKV0-7 | Q9UKV8 | Q9UKX7 | Q9UKY1 | Q9UL03 | Q9UL15 | Q9UL18 | Q9UL51 | Q9ULD2-2 | Q9ULD6 | Q9ULD6-2 | Q9ULG1 | Q9ULH1 | Q9ULM3 | Q9ULR3 | Q9ULV4 | Q9ULV5 | Q9ULW0 | Q9UM11 | Q9UM11-2 | Q9UM54 | Q9UN81 | Q9UNF1 | Q9UPM9 | Q9UPN3 | Q9UPP1 | Q9UPS6 | Q9UPV0 | Q9UPX8 | Q9UPY3 | Q9UPY8 | Q9UPZ3 | Q9UQ16 | Q9UQ88 | Q9UQ90 | Q9UQB3 | Q9UQC1 | Q9UQE7 | Q9UQL6 | Q9WMB5 | Q9WMX2 | Q9WNA9 | Q9WPI5 | Q9Y228 | Q9Y230 | Q9Y247 | Q9Y295 | Q9Y2D5 | Q9Y2D8 | Q9Y2G8 | Q9Y2H6-2 | Q9Y2I6 | Q9Y2J2 | Q9Y2J4 | Q9Y2J4-4 | Q9Y2K6 | Q9Y2K7 | Q9Y2T1 | Q9Y2T2 | Q9Y2T4 | Q9Y2T4-3 | Q9Y2V2 | Q9Y2V7 | Q9Y2W1 | Q9Y2X3 | Q9Y2X7 | Q9Y2X9 | Q9Y3A3 | Q9Y3P9 | Q9Y3Y2 | Q9Y463 | Q9Y468 | Q9Y478 | Q9Y490 | Q9Y4A5 | Q9Y4B4 | Q9Y4K3 | Q9Y4L1 | Q9Y4R8 | Q9Y4W2 | Q9Y4X4 | Q9Y570 | Q9Y572 | Q9Y580 | Q9Y5B9 | Q9Y5J1 | Q9Y5K5 | Q9Y5K8 | Q9Y5N6 | Q9Y5S9 | Q9Y5U2 | Q9Y5V3 | Q9Y5X1 | Q9Y678 | Q9Y6E0 | Q9Y6K9 | Q9Y6Q9 | Q9Y6W6 | U5TQE9
## N 1712
## N non-query with motif: 1601
```

```r
# import into Cytoscape:
## 1. "Import network from file" ending in "_network", choose these types for each column:
### "source" (source)                     "target" (target)                      
### "source_readable" (source attribute)  "target_readable" (target attribute)           
### "human_instances" "top_domain_support4motif_nq" "domain_support4motif_nq" (edge attributes) 
### "type_source" (source attribute)      "type_target" (target attribute)  
## 2. "Import table from files" ending in "_node_table":
### Import data as: Node Table Columns
### Key column for Network: Shared Name
### "node" column in a table is a key column
### choose floating point data type for these columns:
### "p_value" "scale_by" "scale_by2"  (but not for "N_human_instances_bind_domain")
## 3. File -> Import -> style -> motif_domain_network_style.xml (once per cytoscape session)
### In the style pane select "motif_domain_network" style instead of "default" style
## 4. Install yFiles Layout Algorithms. Apps -> App Manager -> find and install
## 5. In the menu Layout select yFiles Hierarchic Layout


# known PDZ motifs
known_PDZ1 = motifSummary(bench_res_cyt, IntAct_net,
             motif = NULL, viral_query = "P03126", human_seed = NULL,
             results_dir = "./results_w_human_dom/")
```

```
## viral query protein with motif: P03126_156_158 = [ST].[LV]$ = Protein E6 (VE6_HPV16)  | P03126_155_158 = [DE]T.[ILV]$ = Protein E6 (VE6_HPV16) 
## 2
## 
## viral targeted human proteins: Q12959 = Disks large homolog 1 | Q14160 = Protein scribble homolog
## 2
## 
## UP support / UP: (P03126_156_158 in Q12959 dataset) 32/46 | (P03126_155_158 in Q14160 dataset) 23/56
## 
## Occurence / SeqNum: (P03126_156_158 in Q12959 dataset) 60/91 | (P03126_155_158 in Q14160 dataset) 32/102
## 
## top-1 domain(s): (P03126_156_158) IPR001478 | (P03126_156_158) IPR008144 | (P03126_156_158) IPR008145 | (P03126_156_158) IPR027417 | (P03126_155_158) IPR001478
## 
## domain_support4motif_nq / protein_with_motif_degree: 18/32 | 6/32 | 6/32 | 8/32 | 11/32
## 
## all with motif: B7Z2Y1 | O60333-3 | O60346 | P03126 | P06427 | P06463 | P0C213 | P16717 | P17386 | P21735 | P22459 | P22460 | P23508 | P24835 | P26555 | P27228 | P33402 | P35222 | P50804 | P53778 | Q01814 | Q13224 | Q14524 | Q15303 | Q15311 | Q6VGS8 | Q7Z628 | Q96A65 | Q96GG9 | Q96PE1 | Q99569 | Q9NS75 | Q9NVW2 | Q9NYB5 | Q9P021 | Q9UDY2 | Q9UQB3
## N 37
## N non-query with motif: 36
```

```r
known_PDZ2 = motifSummary(bench_res_cyt, IntAct_net,
             motif = c("ET..$","[ST].V$", "[DE]T.[ILV]$"), viral_query = "P06463", human_seed = NULL,
             results_dir = "./results_w_human_dom/")
```

```
## viral query protein with motif: P06463_155_158 = [DE]T.[ILV]$ = Protein E6 (VE6_HPV18) 
## 1
## 
## viral targeted human proteins: Q14160 = Protein scribble homolog
## 1
## 
## UP support / UP: (P06463_155_158 in Q14160 dataset) 23/56
## 
## Occurence / SeqNum: (P06463_155_158 in Q14160 dataset) 32/102
## 
## top-1 domain(s): (P06463_155_158) IPR001478
## 
## domain_support4motif_nq / protein_with_motif_degree: 11/81
## 
## all with motif: B7Z2Y1 | O60346 | P03126 | P06427 | P06463 | P0C213 | P21735 | P22460 | P23508 | P24835 | P27228 | P33402 | P35222 | P50804 | P53778 | Q15311 | Q7Z628 | Q9NYB5 | Q9UDY2
## N 19
## N non-query with motif: 18
```

```r
if(FALSE){
    # Candidate PDZ motif
    motifSummary(bench_res_cyt, IntAct_net,
                 motif = NULL, viral_query = "Q2PJP0", human_seed = NULL)
    motifSummary(bench_res_cyt, IntAct_net,
                 motif = c("[DE]..V$", "[ST].[ILV]$", "[ST].V$"), viral_query = "Q2PJP0", human_seed = NULL)
    
    motifSummary(bench_res_cyt, IntAct_net,
                 motif = NULL, viral_query = "P50804", human_seed = NULL)
    
    # candidate SH3 domain motif
    motifSummary(bench_res_cyt, IntAct_net,
                 motif = "P..P.[HKR]", viral_query = NULL, human_seed = NULL)
    
    # candidate WD40 domain motif
    motifSummary(bench_res_cyt, IntAct_net,
                 motif = "DSG", viral_query = NULL, human_seed = NULL)
    motifSummary(bench_res_cyt, IntAct_net,
                 motif = "E.V..G.{0,2}N.{0,1}Q", viral_query = NULL, human_seed = NULL)
    
    # candidate LR.{0,2}G.G.T motif recognised by double-stranded RNA binding domain
    motifSummary(bench_res_cyt, IntAct_net,
                 motif = "LR.{0,2}G.G.T", viral_query = NULL, human_seed = NULL)
    # candidate EI[IV]QQ motif potentially recognised by the EF hand domain
    motifSummary(bench_res_cyt, IntAct_net,
                 motif = "EI[IV]QQ", viral_query = NULL, human_seed = NULL)
    
    # candidate motif L.{0,1}Q.LR potentially recognised by BAG domain
    motifSummary(bench_res_cyt, IntAct_net,
                 motif = "L.{0,1}Q.LR", viral_query = NULL, human_seed = NULL)
    
    ## Check that human proteins in this subnetwork were used to construst "successful" datasets - yes, they were - human protein ID is copied from "interacts_with" column of QSLIMFinder output which contains a seed protein for the dataset
}
```


```r
# lets compare the number of large-scale studies in human vs viral
# load ppi data
IntAct_ppi = loadIntActFTP(dir = "./data_files/IntActRelease_2017Nov13/",
                                                release = "2017Nov13")
# human-viral
all_viral_interaction_ppi = interSpeciesInteractome(taxid1 = 9606, taxid2 = 10239, database = "IntActFTP", format = "tab27",
                                                clean = TRUE, protein_only = TRUE,
                                                MITABdata = IntAct_ppi, directory = "./data_files/",
                                                releaseORdate = "2017Nov13")
# human-human
Full_IntAct_ppi = fullInteractome(taxid = 9606, database = "IntActFTP", format = "tab27",
                              clean = TRUE, protein_only = TRUE,
                              MITABdata = IntAct_ppi, directory = "./data_files/",
                                                releaseORdate = "2017Nov13")


# number of interactions per study
table(all_viral_interaction_ppi$data$Publication_Identifiers)[order(table(all_viral_interaction_ppi$data$Publication_Identifiers), decreasing = T)]
```

```
## 
##       28169297       22810586       22810585       22190034       21911578 
##           3942           3599           1735           1186            902 
##       22761572       25544563       23853584       22458338       18985028 
##            775            579            460            379            354 
##       22028648       20102225       22593156       20064372       21911577 
##            235            232            224            213            209 
##       25616068 unassigned1311       24136289       21994455       17446270 
##            206            206            180            177            171 
##       22014111       24169621       24965446       18457437       24550280 
##            164            156            136            132            118 
##       25158218       26439010       23487458       26808496       21156026 
##            111            108            100             90             84 
##       15140983       15784897       24955142       22522808       20953506 
##             67             66             66             59             55 
##       14993658       22496234       18667692       10487758       12482659 
##             47             46             43             35             32 
##       15491611       22419161       25556234       15607035       24412650 
##             31             31             31             29             28 
##       22124328       24402281       12970734       17696407       12368347 
##             27             27             26             26             24 
##       16112766       24008843       25740999       15567442       20133758 
##             24             24             24             23             23 
##       12692227       16061792       16374509       22072710       12444549 
##             22             22             22             22             21 
##       16601680       17267598       19017798       24284321        8551580 
##             21             21             21             21             21 
##       15846844       16135803       19201886       20088881       21914078 
##             20             20             20             20             20 
##       22056774       22301138       11119584       12970408       15452250 
##             20             20             19             19             19 
##       15983032       21957124       25210169       25599645       10365963 
##             19             19             19             19             18 
##       16227268       18256156       18305892       19403670       23636254 
##             18             18             18             18             18 
##       26771495        9827557       11414803       16140784       16763564 
##             18             18             17             17             17 
##       16917502       20601427       24009510       11907332       16581912 
##             17             17             17             16             16 
##       18922877       24018270       24428437        8627785        9512561 
##             16             16             16             16             16 
##        9651361       14623081       17229681       18946490       23861867 
##             16             15             15             15             15 
##       24478431       25118280       25854864        9683573        9705868 
##             15             15             15             15             15 
##       12427757       12522210        1331501       15175323       16033967 
##             14             14             14             14             14 
##       16237761       16940534       21987769       22776252       23935497 
##             14             14             14             14             14 
##       25944354       11744716       15047801       15084397       15507604 
##             14             13             13             13             13 
##       16983096       17239418       18024891       22685230       22904686 
##             13             13             13             13             13 
##       24418534       24802111       26095031        9223490       10074132 
##             13             13             13             13             12 
##       10722738       12186919       12595532       14557665       15194749 
##             12             12             12             12             12 
##       15364941       16044149       16275648       16534844       18172216 
##             12             12             12             12             12 
##       19008854       19515777       20007272       22640416       23851574 
##             12             12             12             12             12 
##       24263861       24462683       24973451       25253337        9175835 
##             12             12             12             12             12 
##       10381623       10544080       14985081       15824310       16227293 
##             11             11             11             11             11 
##       17412836       17567694       18226242       18498651       22102817 
##             11             11             11             11             11 
##       22416134       22872640       22902366       23576507       26567527 
##             11             11             11             11             11 
##       10228159       10400713       12604805       12743606       12783858 
##             10             10             10             10             10 
##       12813456       14644612       15113910       15669058       16226257 
##             10             10             10             10             10 
##       16249186       16844119       17602943       17728244       18636090 
##             10             10             10             10             10 
##       20133869       20701745       20924359       21297162       21532619 
##             10             10             10             10             10 
##       21667337       21806988       22238231       22301157       23892143 
##             10             10             10             10             10 
##       24335286        7853498        8617242       10675342       11713253 
##             10             10             10              9              9 
##       11733528       12198176       15989969       17024179       17070091 
##              9              9              9              9              9 
##       17686842       18216108       19658097       22045669       22616990 
##              9              9              9              9              9 
##       23638199       24557838       24843023       25122793        9143277 
##              9              9              9              9              9 
##        9188637        9847329       10849009       11086025       11277703 
##              9              9              8              8              8 
##       11571640       12200142       12552007       12649209       12706723 
##              8              8              8              8              8 
##       14524621       14675634       16007207       16103886       16530520 
##              8              8              8              8              8 
##       17023018       17159145       17350572       17932485       18160438 
##              8              8              8              8              8 
##       18688256       19164291       19264780       20639899       20727575 
##              8              8              8              8              8 
##       20862261       21559518       21795333       22334672       22911572 
##              8              8              8              8              8 
##       23096996       23427160       23431397       23552410       23593007 
##              8              8              8              8              8 
##       23986595       24237704       24443581        7624774        9660940 
##              8              8              8              8              8 
##        9687513        9881977       10523853       10644344       10702232 
##              8              8              7              7              7 
##       11842245       11943871       12199906       12730672       15303969 
##              7              7              7              7              7 
##       15567440       15841462       16275649       16322229       16963744 
##              7              7              7              7              7 
##       17360488       17909182       17939992       18198944       18408009 
##              7              7              7              7              7 
##       18682563       18775702       18802093       19322197       19454348 
##              7              7              7              7              7 
##       19605477       19843693       19913487       19966799       20075863 
##              7              7              7              7              7 
##       20562857       20869963       21378963       24643253       24696485 
##              7              7              7              7              7 
##       24706939       25312088       25772236       26559832       26596467 
##              7              7              7              7              7 
##       10082529       10336476       10799599       10873631       11689660 
##              6              6              6              6              6 
##       12359749       12445808       12730686       12740913       14517264 
##              6              6              6              6              6 
##       14528301       15367624       15564746       15611062       15862300 
##              6              6              6              6              6 
##       16439990       16611888       16760386       16799465       17314515 
##              6              6              6              6              6 
##       18266467       18337511       18391203       18488039       18519040 
##              6              6              6              6              6 
##       18579593       18632959       18667510       19153231       19176627 
##              6              6              6              6              6 
##       19264781       19327355       19555689       19850040       20404187 
##              6              6              6              6              6 
##       20554775        2153075       21880768       21915095       21931555 
##              6              6              6              6              6 
##       23015697       23364696       23393263       23414759       23468591 
##              6              6              6              6              6 
##       23486063       23783631       23840900       23890821       24036449 
##              6              6              6              6              6 
##       24970085       25480565       25745166       26638075        8551581 
##              6              6              6              6              6 
##        8553588        8995654        9284049        9326658        9349482 
##              6              6              6              6              6 
##        9472029        9811832        9882311       10329544       11238858 
##              6              6              6              5              5 
##       11915042       12551970       12609975       12730668       12925958 
##              5              5              5              5              5 
##       14638676       14963118       15044019       15073179       15207618 
##              5              5              5              5              5 
##       15492812       15793585       15837194       16051867       16094715 
##              5              5              5              5              5 
##       16257957       16364321       16380082       16405965       16474402 
##              5              5              5              5              5 
##       16641295       16963558       17005671       17022977       17023019 
##              5              5              5              5              5 
##       17074758       17275127       17314511       17485524       18208323 
##              5              5              5              5              5 
##       18551131       18587047       18725644       19197236       20436455 
##              5              5              5              5              5 
##       20631090       21070952       21245344       21270162       21527913 
##              5              5              5              5              5 
##       21808041       21849455       21878328       22048310       22309289 
##              5              5              5              5              5 
##       22916011       23542348       23675303       23857585       25262680 
##              5              5              5              5              5 
##       26789246       26889034        7565739        9192623        9299613 
##              5              5              5              5              5 
##        9405152        9830007        9881696        9990017       10523841 
##              5              5              5              5              4 
##       10873769       10887206       11099378       11145900       11275986 
##              4              4              4              4              4 
##       11278465       11507205       11559788       11573093       11922617 
##              4              4              4              4              4 
##       11934899       11971900       12592401       12620808       14966293 
##              4              4              4              4              4 
##       15140953       15518812       15604276       15767424       15998638 
##              4              4              4              4              4 
##       16227264       16301520       16311516       16376880       16501113 
##              4              4              4              4              4 
##       16690937       16876117       17301141       17348035       17360743 
##              4              4              4              4              4 
##       17666013       17698809       17947517       18045938       18048335 
##              4              4              4              4              4 
##       18059340       18068675       18301737       18562541       18977328 
##              4              4              4              4              4 
##       19150985       19193793       19435845       19590578       19631645 
##              4              4              4              4              4 
##       19651603       20133580       20162731       20534861       20535204 
##              4              4              4              4              4 
##       20638388       21334391       21900157       21994459       22174692 
##              4              4              4              4              4 
##       22641034       23152499       23201263       23328395       23890813 
##              4              4              4              4              4 
##       26789255       26917592        7479821        7565781        7644508 
##              4              4              4              4              4 
##        7815540        8207801        8380896        8411358        8423815 
##              4              4              4              4              4 
##        9223484        9569025        9608663        9927201       10049825 
##              4              4              4              4              3 
##       10074103       10187771       10196247       10318918       10369867 
##              3              3              3              3              3 
##       10374970       10390359       10618715       10779361       11027293 
##              3              3              3              3              3 
##       11238845       11264384       11836380       11878931       12021324 
##              3              3              3              3              3 
##       12034482       12119400       12163591       12356718       12458223 
##              3              3              3              3              3 
##       12600646       14576168       14732683       14737186       14978742 
##              3              3              3              3              3 
##       15016876       15793001       16051843       16103149       16155203 
##              3              3              3              3              3 
##       16477038       16537421       16840345       16932746       16951253 
##              3              3              3              3              3 
##       17123957       17137594       17301785       17310249       17529992 
##              3              3              3              3              3 
##       17616579       18358807       18474220       18539148       19249677 
##              3              3              3              3              3 
##       19369353       19379743       19745842       19901337       20212144 
##              3              3              3              3              3 
##       20462497       21266548       21378965       21835792       22140357 
##              3              3              3              3              3 
##       22623778       22826234       22978549       23001005       23135708 
##              3              3              3              3              3 
##       23213219       23420847       23507139       23523133       23552417 
##              3              3              3              3              3 
##       23831647       23870315       24418540       25036637       25244949 
##              3              3              3              3              3 
##       26098315       26365801       26405230       28046066        7592647 
##              3              3              3              3              3 
##        7896830        8388547        8648648        9018065        9020106 
##              3              3              3              3              3 
##        9203586        9299485        9311810        9374493        9499053 
##              3              3              3              3              3 
##        9973195       10330164       10347196       10498661       10702269 
##              3              2              2              2              2 
##       10702287       10873756       10944537       11160738       11251081 
##              2              2              2              2              2 
##       11298454       11462050       11509578       11773402       11847120 
##              2              2              2              2              2 
##       12604806       12634384       12821782        1321290       15039538 
##              2              2              2              2              2 
##       15144954       15152192       15221897       15507623       15527772 
##              2              2              2              2              2 
##       15554700       15582666       15629770       15825084       15848177 
##              2              2              2              2              2 
##       15893726       15967037       16129692       16381817       16493710 
##              2              2              2              2              2 
##       16611247       16646632       16945160       17006541       17067581 
##              2              2              2              2              2 
##       17097642       17166906       17166913       17170431       17174309 
##              2              2              2              2              2 
##       17267502       17363894       17382937       17428914       17486070 
##              2              2              2              2              2 
##       17717529       17855521       18005740       18042718       18252829 
##              2              2              2              2              2 
##       18267072       18552826       18604270       18616943       18647839 
##              2              2              2              2              2 
##       18701889       18971339       19073181       19104048       19179289 
##              2              2              2              2              2 
##       19293378       19339969       19339972       19424970       19838188 
##              2              2              2              2              2 
##       19966293       20010840       20080564       20080767       20138174 
##              2              2              2              2              2 
##       20142477       20357769       20802522       21068253       21092281 
##              2              2              2              2              2 
##       21098020       21217702        2138977        2162480       21697476 
##              2              2              2              2              2 
##       22004763       22171259       22483108       22491470       22933281 
##              2              2              2              2              2 
##       23576498       23722113       23741449       24034250       24069433 
##              2              2              2              2              2 
##       24335311       24336198       25170080       25193851       26867650 
##              2              2              2              2              2 
##        7483825        7537265        7542742        7680435        7859737 
##              2              2              2              2              2 
##        8386265        9034339        9050995        9168958        9188558 
##              2              2              2              2              2 
##        9311902        9431994        9811611       10226625       10327051 
##              2              2              2              1              1 
##       10488152       10523825       10567268       10570951       10734313 
##              1              1              1              1              1 
##       10816381       10982355       11042264       11226179       11238466 
##              1              1              1              1              1 
##       11251120       11260720       11511370       11576548       11753645 
##              1              1              1              1              1 
##       11807220       11828325       11877480       11907229       12011057 
##              1              1              1              1              1 
##       12086624       12093920       12419245       12610133       12620801 
##              1              1              1              1              1 
##       12859895       12897781        1328865       14527689       15068796 
##              1              1              1              1              1 
##       15109495       15215856       15254196       15292184       15310821 
##              1              1              1              1              1 
##       15380363       15521019       15809419       15890204       16003391 
##              1              1              1              1              1 
##       16126728       16203724       16423054       16501114       16512683 
##              1              1              1              1              1 
##       16595549       16624878       16823848       16878151       16963564 
##              1              1              1              1              1 
##       16996479        1706474       17076584       17254575       17259306 
##              1              1              1              1              1 
##       17293873       17341466       17531282       17540176       17612295 
##              1              1              1              1              1 
##       17662718       17673209       17681535       17713926       17991777 
##              1              1              1              1              1 
##       18325616       18391209        1870208       18836454       19088272 
##              1              1              1              1              1 
##       19581596       19838190       19887642       19920174       20208545 
##              1              1              1              1              1 
##       20529865       20679221       20708612       21094642       21238945 
##              1              1              1              1              1 
##        2189724       21970979       22863293       22904195       22909819 
##              1              1              1              1              1 
##       23063122       23202587       23219818       23320542       23907398 
##              1              1              1              1              1 
##       23954291       24009866       24269682       24618592       24942588 
##              1              1              1              1              1 
##       25686249        2767055        7588762        7609027        7666533 
##              1              1              1              1              1 
##        7671300        8344494        8381428        8607265        8980234 
##              1              1              1              1              1 
##        9067421        9218412        9271389        9343261        9466902 
##              1              1              1              1              1 
##        9467941        9697699        9794763        9807817        9857076 
##              1              1              1              1              1 
##        9874765 
##              1
```

```r
countBaits = function(data){
    viral_ppi_data = copy(data$data)
    viral_prot = viral_ppi_data[, .(prot = c(IDs_interactor_A, IDs_interactor_B), 
                                    bait_prey = c(bait_prey_status_A, bait_prey_status_B),
                                    pmid = c(Publication_Identifiers, Publication_Identifiers))]
    viral_prot = unique(viral_prot)
    viral_prot[, bait_per_study := uniqueN(prot), by = .(pmid, bait_prey)]
    viral_studies = unique(viral_prot[, .(bait_per_study, pmid, bait_prey)])
    viral_studies[order(bait_per_study)][bait_prey == "bait"]
}

viral_studies = countBaits(all_viral_interaction_ppi)
human_studies = countBaits(Full_IntAct_ppi)
```

## Time spend per dataset


```r
motifs = fread("/Users/vk7/Desktop/ebi_projects/viral_project/qslimfinder.Full_IntAct3cloudfixF.FALSE/result/main_result.txt", stringsAsFactors = F)
motifs[, time_in_sec := as.integer(as.POSIXct(strptime(RunTime, "%H:%M:%S")))]
motifs[, time_in_sec := time_in_sec - min(time_in_sec)]
dataset_times = unique(motifs[, .(Dataset, time_in_sec, SeqNum)])
hist(dataset_times$time_in_sec)
plot(dataset_times[SeqNum > 170 & SeqNum < 200]$SeqNum, dataset_times[SeqNum > 170 & SeqNum < 200]$time_in_sec)
mean(dataset_times[SeqNum > 170 & SeqNum < 200]$time_in_sec)
```


```r
session_info. = devtools::session_info()
session_info.
```

```
##  setting  value                       
##  version  R version 3.4.4 (2018-03-15)
##  system   x86_64, linux-gnu           
##  ui       X11                         
##  language (EN)                        
##  collate  en_GB.UTF-8                 
##  tz       Europe/London               
##  date     2018-08-27                  
## 
##  package              * version   date      
##  AnnotationDbi          1.40.0    2018-03-16
##  assertthat             0.2.0     2017-04-11
##  backports              1.1.2     2017-12-13
##  base                 * 3.4.4     2018-03-16
##  bindr                  0.1.1     2018-03-13
##  bindrcpp               0.2.2     2018-03-29
##  Biobase                2.38.0    2018-03-16
##  BiocGenerics         * 0.24.0    2018-03-16
##  BiocParallel           1.12.0    2018-03-26
##  biomaRt              * 2.34.2    2018-03-25
##  Biostrings           * 2.46.0    2018-03-16
##  bit                    1.1-14    2018-05-29
##  bit64                  0.9-7     2017-05-08
##  bitops                 1.0-6     2013-08-17
##  blob                   1.1.1     2018-03-25
##  caTools                1.17.1.1  2018-07-20
##  clustermq            * 0.8.4.1   2018-06-29
##  colorspace             1.3-2     2016-12-14
##  compiler               3.4.4     2018-03-16
##  crayon                 1.3.4     2017-09-16
##  curl                   3.2       2018-03-28
##  data.table           * 1.11.4    2018-05-27
##  datasets             * 3.4.4     2018-03-16
##  DBI                    1.0.0     2018-05-02
##  DelayedArray           0.4.1     2018-03-16
##  devtools               1.13.6    2018-06-27
##  digest                 0.6.15    2018-01-28
##  downloader             0.4       2015-07-09
##  dplyr                  0.7.6     2018-06-29
##  DT                     0.4       2018-01-30
##  evaluate               0.11      2018-07-17
##  gdata                  2.18.0    2017-06-06
##  GenomeInfoDb         * 1.14.0    2018-03-16
##  GenomeInfoDbData       1.0.0     2018-03-16
##  GenomicAlignments      1.14.2    2018-08-01
##  GenomicRanges        * 1.30.3    2018-03-16
##  GGally                 1.4.0     2018-05-17
##  ggplot2              * 3.0.0     2018-07-03
##  glue                   1.3.0     2018-07-17
##  gplots                 3.0.1     2016-03-30
##  graphics             * 3.4.4     2018-03-16
##  grDevices            * 3.4.4     2018-03-16
##  grid                 * 3.4.4     2018-03-16
##  gridExtra            * 2.3       2017-09-09
##  gsubfn                 0.7       2018-03-16
##  gtable                 0.2.0     2016-02-26
##  gtools                 3.8.1     2018-06-26
##  hms                    0.4.2     2018-03-10
##  htmltools              0.3.6     2017-04-28
##  htmlwidgets            1.2       2018-04-19
##  httr                 * 1.3.1     2017-08-20
##  igraph                 1.2.2     2018-07-27
##  IRanges              * 2.12.0    2018-03-16
##  jsonlite               1.5       2017-06-01
##  KernSmooth             2.23-15   2015-06-29
##  knitr                  1.20      2018-02-20
##  labeling               0.3       2014-08-23
##  lattice                0.20-35   2017-03-25
##  lazyeval               0.2.1     2017-10-29
##  magrittr               1.5       2014-11-22
##  Matrix                 1.2-14    2018-04-09
##  matrixStats            0.54.0    2018-07-23
##  memoise                1.1.0     2017-04-21
##  methods              * 3.4.4     2018-03-16
##  MItools              * 0.1.40    2018-08-26
##  munsell                0.5.0     2018-06-12
##  ontologyIndex          2.4       2017-02-06
##  parallel             * 3.4.4     2018-03-16
##  pillar                 1.3.0     2018-07-14
##  pkgconfig              2.0.1     2017-03-21
##  plyr                 * 1.8.4     2016-06-08
##  prettyunits            1.0.2     2015-07-13
##  progress               1.2.0     2018-06-14
##  proto                  1.0.0     2016-10-29
##  PSICQUIC             * 1.16.4    2018-03-26
##  purrr                  0.2.5     2018-05-29
##  qvalue                 2.10.0    2018-03-26
##  R.methodsS3            1.7.1     2016-02-16
##  R.oo                   1.22.0    2018-04-22
##  R.utils                2.6.0     2017-11-05
##  R6                     2.2.2     2017-06-17
##  RColorBrewer         * 1.1-2     2014-12-07
##  Rcpp                   0.12.18   2018-07-23
##  RCurl                  1.95-4.11 2018-07-15
##  reshape                0.8.7     2017-08-06
##  reshape2               1.4.3     2017-12-11
##  rlang                  0.2.1     2018-05-30
##  rmarkdown            * 1.10      2018-06-11
##  ROCR                   1.0-7     2015-03-26
##  rprojroot              1.3-2     2018-01-03
##  Rsamtools              1.30.0    2018-03-26
##  RSQLite                2.1.1     2018-05-06
##  rtracklayer            1.38.3    2018-03-26
##  S4Vectors            * 0.16.0    2018-03-16
##  scales                 0.5.0     2017-08-24
##  splines                3.4.4     2018-03-16
##  stats                * 3.4.4     2018-03-16
##  stats4               * 3.4.4     2018-03-16
##  stringi                1.2.4     2018-07-20
##  stringr                1.3.1     2018-05-10
##  SummarizedExperiment   1.8.1     2018-03-16
##  tibble                 1.4.2     2018-01-22
##  tidyselect             0.2.4     2018-02-26
##  tools                  3.4.4     2018-03-16
##  utils                * 3.4.4     2018-03-16
##  withr                  2.1.2     2018-03-15
##  XML                    3.98-1.13 2018-08-01
##  XVector              * 0.18.0    2018-03-16
##  yaml                   2.2.0     2018-07-25
##  zlibbioc               1.24.0    2018-03-16
##  source                        
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  local                         
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  Bioconductor                  
##  Bioconductor                  
##  Bioconductor                  
##  Bioconductor                  
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  local                         
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  cran (@1.11.4)                
##  local                         
##  CRAN (R 3.4.4)                
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  cran (@0.4)                   
##  CRAN (R 3.4.4)                
##  cran (@0.4)                   
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  Bioconductor                  
##  Bioconductor                  
##  Bioconductor                  
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  cran (@3.0.1)                 
##  local                         
##  local                         
##  local                         
##  CRAN (R 3.4.4)                
##  cran (@0.7)                   
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  cran (@0.4.2)                 
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  local                         
##  Github (vitkl/MItools@87f2b47)
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  local                         
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  cran (@1.2.0)                 
##  cran (@1.0.0)                 
##  Bioconductor                  
##  cran (@0.2.5)                 
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  cran (@0.12.18)               
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  cran (@0.2.1)                 
##  CRAN (R 3.4.4)                
##  cran (@1.0-7)                 
##  CRAN (R 3.4.4)                
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  Bioconductor                  
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  local                         
##  local                         
##  local                         
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  local                         
##  local                         
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  Bioconductor
```

```r
save(list = ls(), file=filename)
```
