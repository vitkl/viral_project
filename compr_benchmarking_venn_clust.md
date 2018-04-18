---
title: "Comprehensive benchmarking of motif instances (all types) - only IntAct and Vidal"
author: "Vitalii Kleshchevnikov"
date: "18/10/2017"
output: 
  html_document: 
    keep_md: yes
    toc: yes
---


```r
knitr::opts_chunk$set(echo = FALSE, fig.width = 12, fig.height = 8, warning = FALSE, message = FALSE)

packages = c("MItools", "RColorBrewer", "devtools", "grid","gridExtra", "ggplot2", "VennDiagram")
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
normalise = F # if normalised Sig threshold of 0.3 becomes 1 (but this equalises Sig == 0.3 to true 1 - that is when motif was not found)
#ROC:
measure1 = "prec" # "tpr" or "prec"
measure2 = "rec" # "fpr" or "rec"
both_metrics_vs_cutoff = T
single_metric = c("auc","prbe")[1]
single_metric_name = c("Median AUC", "Prec-rec break-even")[1]

datasets = c("qslimfinder.Full_IntAct3.FALSE",
             "qslimfinder.Vidal3.FALSE",
             "qslimfinder.all_viral_interaction3.FALSE")
descriptions = c("human network (full IntAct) searched \nfor motifs present in viral proteins",
                 "human network (Vidal's data only) searched \nfor motifs present in viral proteins",
                 "human network (human-viral data only) searched \nfor motifs present in viral proteins")
motif_setup_month = c("201802",
                      "201802",
                      "201802")
```

## Overview

### 
1. filter_by_domain_data: NULL or "p.value < 0.5"
2. neg_set: "all_instances", "all_proteins"
3. motif_pval_cutoff: 1, precision == recall, precision > 0.5


```
## gTree[venn]
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/not_filt_Sig-1.png)<!-- -->

```
## TableGrob (6 x 1) "arrange": 2 grobs
##   z     cells    name               grob
## 1 1 (1-1,1-1) arrange text[GRID.text.29]
## 2 2 (2-6,1-1) arrange        gTree[venn]
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/not_filt_Sig-2.png)<!-- -->![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/not_filt_Sig-3.png)<!-- -->![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/not_filt_Sig-4.png)<!-- -->![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/not_filt_Sig-5.png)<!-- -->![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/not_filt_Sig-6.png)<!-- -->![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/not_filt_Sig-7.png)<!-- -->![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/not_filt_Sig-8.png)<!-- -->![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/not_filt_Sig-9.png)<!-- -->![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/not_filt_Sig-10.png)<!-- -->![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/not_filt_Sig-11.png)<!-- -->![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/not_filt_Sig-12.png)<!-- -->![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_venn_clust_files/figure-html/not_filt_Sig-13.png)<!-- -->



```
## [1] "2018-04-18"
```

```
##  setting  value                       
##  version  R version 3.4.4 (2018-03-15)
##  system   x86_64, linux-gnu           
##  ui       X11                         
##  language (EN)                        
##  collate  en_GB.UTF-8                 
##  tz       Europe/London               
##  date     2018-04-18                  
## 
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
##  data.table           * 1.10.4-3  2017-10-27
##  datasets             * 3.4.4     2018-03-16
##  DBI                    0.8       2018-03-02
##  DelayedArray           0.4.1     2018-03-16
##  devtools               1.13.5    2018-02-18
##  digest                 0.6.15    2018-01-28
##  downloader             0.4       2015-07-09
##  DT                     0.4       2018-01-30
##  evaluate               0.10.1    2017-06-24
##  futile.logger          1.4.3     2016-07-10
##  futile.options         1.0.0     2010-04-06
##  gdata                  2.18.0    2017-06-06
##  gdtools              * 0.1.7     2018-02-27
##  GenomeInfoDb           1.14.0    2018-03-16
##  GenomeInfoDbData       1.0.0     2018-03-16
##  GenomicAlignments      1.14.1    2018-03-26
##  GenomicRanges          1.30.3    2018-03-16
##  GGally                 1.3.2     2017-08-02
##  ggplot2              * 2.2.1     2016-12-30
##  gplots                 3.0.1     2016-03-30
##  graphics             * 3.4.4     2018-03-16
##  grDevices            * 3.4.4     2018-03-16
##  grid                 * 3.4.4     2018-03-16
##  gridExtra            * 2.3       2017-09-09
##  gsubfn                 0.7       2018-03-16
##  gtable                 0.2.0     2016-02-26
##  gtools                 3.5.0     2015-05-29
##  htmltools              0.3.6     2017-04-28
##  htmlwidgets            1.0       2018-01-20
##  httr                 * 1.3.1     2017-08-20
##  igraph                 1.2.1     2018-03-10
##  IRanges              * 2.12.0    2018-03-16
##  jsonlite               1.5       2017-06-01
##  KernSmooth             2.23-15   2015-06-29
##  knitr                  1.20      2018-02-20
##  lambda.r               1.2       2017-09-16
##  lattice                0.20-35   2017-03-25
##  lazyeval               0.2.1     2017-10-29
##  magrittr               1.5       2014-11-22
##  Matrix                 1.2-12    2017-11-30
##  matrixStats            0.53.1    2018-02-11
##  memoise                1.1.0     2017-04-21
##  methods              * 3.4.4     2018-03-16
##  MItools              * 0.1.36    2018-04-17
##  munsell                0.4.3     2016-02-13
##  ontologyIndex          2.4       2017-02-06
##  parallel             * 3.4.4     2018-03-16
##  pillar                 1.2.1     2018-02-27
##  pkgconfig              2.0.1     2017-03-21
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
##  rtracklayer            1.38.3    2018-03-26
##  S4Vectors            * 0.16.0    2018-03-16
##  scales                 0.5.0     2017-08-24
##  splines                3.4.4     2018-03-16
##  stats                * 3.4.4     2018-03-16
##  stats4               * 3.4.4     2018-03-16
##  stringi                1.1.7     2018-03-12
##  stringr                1.3.0     2018-02-19
##  SummarizedExperiment   1.8.1     2018-03-16
##  svglite                1.2.1     2017-09-11
##  tibble                 1.4.2     2018-01-22
##  tools                  3.4.4     2018-03-16
##  utils                * 3.4.4     2018-03-16
##  VennDiagram            1.6.20    2018-03-28
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
##  local                         
##  CRAN (R 3.4.4)                
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  cran (@0.4)                   
##  cran (@0.4)                   
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
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
##  CRAN (R 3.4.4)                
##  cran (@0.7)                   
##  CRAN (R 3.4.4)                
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
##  CRAN (R 3.4.4)                
##  local                         
##  Github (vitkl/MItools@dbb074c)
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  local                         
##  CRAN (R 3.4.4)                
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
##  CRAN (R 3.4.4)                
##  local                         
##  local                         
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  CRAN (R 3.4.4)                
##  Bioconductor                  
##  CRAN (R 3.4.4)                
##  Bioconductor
```
