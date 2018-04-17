---
title: "Motif_search_strategies"
author: "Vitalii Kleshchevnikov"
date: "29/11/2017"
output: 
  html_document: 
    keep_md: yes
    toc: yes
---



# Introduction

## Strategies

1. Optimal search settings (default):  
- slimlen= 5: Maximum length of SLiMs to return (no. non-wildcard positions)  
- maxwild= 2: Maximum number of consecutive wildcard positions to allow  

2. Comparing interaction datasets (2):  
- Full IntAct  
- Vidal  
- all_viral_interaction  

3. Protein-centered strategy (1)  

4. QSLIMFinder (1)  

Total: 1 \* 1 \* 2 \* 1 \* 1 = 3 conditions  

## load PPI data


```r
# load ppi data
IntAct = loadIntActFTP(dir = "./data_files/IntActRelease_2017Nov13/",
                                                release = "2017Nov13")
```

```
## ... loading local copy ...
```

```
## Read 0.0% of 790620 rowsRead 6.3% of 790620 rowsRead 16.4% of 790620 rowsRead 24.0% of 790620 rowsRead 32.9% of 790620 rowsRead 37.9% of 790620 rowsRead 46.8% of 790620 rowsRead 55.7% of 790620 rowsRead 65.8% of 790620 rowsRead 67.0% of 790620 rowsRead 75.9% of 790620 rowsRead 86.0% of 790620 rowsRead 96.1% of 790620 rowsRead 790620 rows and 42 (of 42) columns from 2.995 GB file in 00:00:24
```

```r
# human-viral
all_viral_interaction3 = interSpeciesInteractome(taxid1 = 9606, taxid2 = 10239, database = "IntActFTP", format = "tab27",
                                                clean = TRUE, protein_only = TRUE,
                                                MITABdata = IntAct, directory = "./data_files/",
                                                releaseORdate = "2017Nov13")
```

```
## Read 5.3% of 189363 rowsRead 189363 rows and 11 (of 11) columns from 0.030 GB file in 00:00:03
```

```r
# human-human
Full_IntAct3 = fullInteractome(taxid = 9606, database = "IntActFTP", format = "tab27",
                              clean = TRUE, protein_only = TRUE,
                              MITABdata = IntAct, directory = "./data_files/",
                                                releaseORdate = "2017Nov13")

Vidal3 = subsetMITABbyPMIDs(MITABdata = Full_IntAct3,
                           PMIDs = c("25416956", "unassigned1304"))
```

## Set standard options


```r
# set standard options
myPPInetwork2SLIMFinder = function(dataset_name = "SLIMFinder", slimlen = 5, maxwild = 2,
                                   interaction_main_set, interaction_query_set,
                                   center_domains = F, analysis_type = "qslimfinder"){
    PPInetwork2SLIMFinder(dataset_name = dataset_name,
                          interaction_main_set = interaction_main_set,
                          interaction_query_set = interaction_query_set,
                          analysis_type = analysis_type,
                          options = paste0("dismask=T consmask=F cloudfix=T probcut=0.3 minwild=0 maxwild=",maxwild," slimlen=",slimlen," alphahelix=F maxseq=800 savespace=0 iuchdir=T extras=2"),
                          domain_res_file = "./processed_data_files/what_we_find_VS_ELM_clust20171019.RData",
                          domain_results_obj = "res_count",
                          center_domains = center_domains,
                          fasta_path = "./data_files/all_human_viral_proteins.fasta",
                          main_set_only = F,
                          domain_pvalue_cutoff = 1, # 0.000005 for test run, 1 for full run
                          SLIMFinder_dir = paste0("./",dataset_name,"/"),
                          LSF_project_path = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/",
                          software_path = "../software/cluster/",
                          length_set1_min = 2, # selected based on benchmarking: minimal set size that returns motifs with Sig < 0.05 = 4+1
                          length_set2_min = 1,
                          write_log = T,
                          N_seq = 200,
                          seed_list = NULL, # NULL for full run, seed_list for test run
                          memory_start = 350,
                          memory_step = 2000)
}

# Download ELM data
elm_filename = paste0("./data_files/", data.table::year(Sys.Date()), 
                      "elms_index.tsv")
if (!file.exists(elm_filename)) {
    download.file("http://elm.eu.org/elms/elms_index.tsv", elm_filename)   
}
```

## Define parameters


```r
analysis_type = c("qslimfinder")
datasets = c("Full_IntAct3", "Vidal3", "all_viral_interaction3")
domain_centered = c(F) # c(F, T)
#slimlen = c(5, 8, 10)
#maxwild = c(2, 3)

dataset_names = data.table(analysis_type = character(), datasets = character(),
                           domain_centered = logical(), dataset_names = character(), 
                           slimlen = integer(), maxwild = integer())
for(a in analysis_type){
    for(d in datasets){
        for(cen in domain_centered){
            dataset_names_temp = data.table(analysis_type = a, datasets = d,
                                            domain_centered = cen, dataset_names = paste(a,d,cen, sep = "."), 
                                            slimlen = 5, maxwild = 2)
            dataset_names = rbind(dataset_names, dataset_names_temp)
        }
    }
}
dataset_names[grepl("BioPlex",datasets), query_dataset := "all_viral_interaction_BioPlex"]
dataset_names[!grepl("BioPlex",datasets), query_dataset := "all_viral_interaction3"]

# remove already finished names
R_data_names = readLines("./RData_from_Motif_search_strategies")
R_data_names = gsub("\\./processed_data_files/QSLIMFinder_instances_h2v_","", R_data_names)
R_data_names = gsub("_clust201712\\.RData\\.zip","", R_data_names)
R_data_names = gsub("_clust201801\\.RData\\.zip","", R_data_names)
R_data_names = gsub("_clust201802\\.RData\\.zip","", R_data_names)
dataset_names_not_done = dataset_names[!dataset_names %in% R_data_names]
```

## Perform motif search


```r
if(nrow(dataset_names_not_done) >= 1){
    # set up parallel processing
    # create cluster
    cores = 3
    cl <- makeCluster(cores)
    # get library support needed to run the code
    clusterEvalQ(cl, {library(MItools); library(rtracklayer)})
    # put objects in place that might be needed for the code
    clusterExport(cl, ls(), envir=environment())
    
    R_data = parSapplyLB(cl = cl, X = 1:nrow(dataset_names_not_done), FUN = function(i){
        R_data_name = myPPInetwork2SLIMFinder(dataset_name = dataset_names_not_done$dataset_names_not_done[i], 
                                              slimlen = dataset_names_not_done$slimlen[i],
                                              maxwild = dataset_names_not_done$maxwild[i],
                                              interaction_main_set = eval(parse(text = dataset_names_not_done$datasets[i])),
                                              interaction_query_set = eval(parse(text = dataset_names_not_done$query_dataset[i])),
                                              center_domains = dataset_names_not_done$domain_centered[i], 
                                              analysis_type = dataset_names_not_done$analysis_type[i])
        
        zip(zipfile = paste0(R_data_name, ".zip"), files = R_data_name)
        unlink(R_data_name)
        
        R_data_names = readLines("./RData_from_Motif_search_strategies")
        R_data_names = c(R_data_names, paste0(R_data_name, ".zip"))
        write(R_data_names, "./RData_from_Motif_search_strategies")
        
        paste0(R_data_name, ".zip")
    }, simplify = TRUE, USE.NAMES = TRUE) # parSapplyLB is parSapply that accounts for different time each dataset might take
    
    stopCluster(cl)
    # paths to RData files containing all objects created by all runs of this pipeline
    print(R_data)
}
```

## Compare discovered motifs to known motifs


```r
resultdirs = paste0("./", dataset_names$dataset_names, "/result/")
software_path = "../software/cluster/"
CompariMotif3_dburl = "http://elm.eu.org/elms/elms_index.tsv"
CompariMotif3_dbpath = "./data_files/"
LSF_project_path = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/"
# compare discovered motifs to ELM
for (resultdir in resultdirs) {
    runCompariMotif3(input_file = paste0(resultdir, "motifs.txt"),
                 slimpath = paste0(software_path, "slimsuite/tools/"),
                 dbpath = CompariMotif3_dbpath,
                 dburl = CompariMotif3_dburl,
                 run = T, with = "db",
                 out_file = paste0(resultdir, "comparimotif.tdt"),
                 LSF_project_path = LSF_project_path)
}
```


```r
Sys.Date. = Sys.Date()
Sys.Date.
```

```
## [1] "2018-04-03"
```

```r
session_info. = devtools::session_info()
session_info.
```

```
## Session info -------------------------------------------------------------
```

```
##  setting  value                       
##  version  R version 3.4.4 (2018-03-15)
##  system   x86_64, linux-gnu           
##  ui       X11                         
##  language (EN)                        
##  collate  en_GB.UTF-8                 
##  tz       Europe/London               
##  date     2018-04-03
```

```
## Packages -----------------------------------------------------------------
```

```
##  package              * version   date      
##  AnnotationDbi          1.40.0    2018-03-16
##  assertthat             0.2.0     2017-04-11
##  backports              1.1.2     2017-12-13
##  base                 * 3.4.4     2018-03-16
##  Biobase                2.38.0    2018-03-16
##  BiocGenerics         * 0.24.0    2018-03-16
##  BiocParallel           1.12.0    2018-03-26
##  biomaRt              * 2.34.2    2018-03-25
##  Biostrings           * 2.46.0    2018-03-16
##  bit                    1.1-12    2014-04-09
##  bit64                  0.9-7     2017-05-08
##  bitops                 1.0-6     2013-08-17
##  blob                   1.1.1     2018-03-25
##  caTools                1.17.1    2014-09-10
##  colorspace             1.3-2     2016-12-14
##  compiler               3.4.4     2018-03-16
##  curl                   3.1       2017-12-12
##  data.table           * 1.10.4-3  2017-10-27
##  datasets             * 3.4.4     2018-03-16
##  DBI                    0.8       2018-03-02
##  DelayedArray           0.4.1     2018-03-16
##  devtools               1.13.5    2018-02-18
##  digest                 0.6.15    2018-01-28
##  downloader             0.4       2015-07-09
##  DT                     0.4       2018-01-30
##  evaluate               0.10.1    2017-06-24
##  gdata                  2.18.0    2017-06-06
##  GenomeInfoDb         * 1.14.0    2018-03-16
##  GenomeInfoDbData       1.0.0     2018-03-16
##  GenomicAlignments      1.14.1    2018-03-26
##  GenomicRanges        * 1.30.3    2018-03-16
##  GGally                 1.3.2     2017-08-02
##  ggplot2                2.2.1     2016-12-30
##  gplots                 3.0.1     2016-03-30
##  graphics             * 3.4.4     2018-03-16
##  grDevices            * 3.4.4     2018-03-16
##  grid                   3.4.4     2018-03-16
##  gsubfn                 0.7       2018-03-16
##  gtable                 0.2.0     2016-02-26
##  gtools                 3.5.0     2015-05-29
##  htmltools              0.3.6     2017-04-28
##  htmlwidgets            1.0       2018-01-20
##  httr                 * 1.3.1     2017-08-20
##  IRanges              * 2.12.0    2018-03-16
##  jsonlite               1.5       2017-06-01
##  KernSmooth             2.23-15   2015-06-29
##  knitr                  1.20      2018-02-20
##  lattice                0.20-35   2017-03-25
##  lazyeval               0.2.1     2017-10-29
##  magrittr               1.5       2014-11-22
##  Matrix                 1.2-12    2017-11-30
##  matrixStats            0.53.1    2018-02-11
##  memoise                1.1.0     2017-04-21
##  methods              * 3.4.4     2018-03-16
##  MItools              * 0.1.36    2018-04-02
##  munsell                0.4.3     2016-02-13
##  ontologyIndex          2.4       2017-02-06
##  parallel             * 3.4.4     2018-03-16
##  pillar                 1.2.1     2018-02-27
##  plyr                 * 1.8.4     2016-06-08
##  prettyunits            1.0.2     2015-07-13
##  progress               1.1.2     2016-12-14
##  proto                  1.0.0     2016-10-29
##  PSICQUIC             * 1.16.4    2018-03-26
##  qvalue                 2.10.0    2018-03-26
##  R.methodsS3            1.7.1     2016-02-16
##  R.oo                   1.21.0    2016-11-01
##  R.utils                2.6.0     2017-11-05
##  R6                     2.2.2     2017-06-17
##  RColorBrewer           1.1-2     2014-12-07
##  Rcpp                   0.12.16   2018-03-13
##  RCurl                  1.95-4.10 2018-01-04
##  reshape                0.8.7     2017-08-06
##  reshape2               1.4.3     2017-12-11
##  rlang                  0.2.0     2018-02-20
##  rmarkdown            * 1.9       2018-03-01
##  ROCR                   1.0-7     2015-03-26
##  rprojroot              1.3-2     2018-01-03
##  Rsamtools              1.30.0    2018-03-26
##  RSQLite                2.0       2017-06-19
##  rtracklayer          * 1.38.3    2018-03-26
##  S4Vectors            * 0.16.0    2018-03-16
##  scales                 0.5.0     2017-08-24
##  splines                3.4.4     2018-03-16
##  stats                * 3.4.4     2018-03-16
##  stats4               * 3.4.4     2018-03-16
##  stringi                1.1.7     2018-03-12
##  stringr                1.3.0     2018-02-19
##  SummarizedExperiment   1.8.1     2018-03-16
##  tibble                 1.4.2     2018-01-22
##  tools                  3.4.4     2018-03-16
##  utils                * 3.4.4     2018-03-16
##  withr                  2.1.2     2018-03-15
##  XML                    3.98-1.10 2018-02-19
##  XVector              * 0.18.0    2018-03-16
##  yaml                   2.1.18    2018-03-08
##  zlibbioc               1.24.0    2018-03-16
##  source                        
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  local                         
##  Bioconductor                  
##  Bioconductor                  
##  Bioconductor                  
##  Bioconductor                  
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  cran (@1.17.1)                
##  CRAN (R 3.4.4)                
##  local                         
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  local                         
##  CRAN (R 3.4.4)                
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  cran (@0.4)                   
##  cran (@0.4)                   
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  Bioconductor                  
##  Bioconductor                  
##  Bioconductor                  
##  Bioconductor                  
##  cran (@1.3.2)                 
##  CRAN (R 3.4.4)                
##  cran (@3.0.1)                 
##  local                         
##  local                         
##  local                         
##  cran (@0.7)                   
##  CRAN (R 3.4.4)                
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
##  local                         
##  Github (vitkl/MItools@d99f98b)
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  local                         
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  cran (@1.0.0)                 
##  Bioconductor                  
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
##  local                         
##  local                         
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  Bioconductor
```

```r
filename = paste0("./processed_data_files/Motif_search_strategies_IntAct_Vidal_viral2",
                  format(Sys.Date(), "%Y%m"),"_",
                  paste0(analysis_type,collapse = "_"),
                  paste0(datasets,collapse = "_"),
                  paste0(domain_centered,collapse = "_"),
                  ".RData")
save(list = ls(), file=filename)
```
