# Domain enrichment by permutation
Vitalii Kleshchevnikov  
25/07/2017  

2017-08-01

## Installing/loading packages



## Setting up how many cores to use for analysis and fold-enrichment vs freqency


```r
# cores to use for multiprocessing: is the number of cores is small -> we are not on a cluster -> use all; if the number of cores is large -> we are on cluster -> use only 15 or 7 cores
if(detectCores() <= 4) cores_to_use = detectCores() - 1
if(detectCores() > 4) cores_to_use = 15

# Choose statistic for permutations:
# frequency in a set (of human interacting partners of viral protein) or fold_enrichment? frequency is used if TRUE 
frequency = F
```

## Null hypothesis description

## How often do we observe specific fold enrichment value among interactors of a viral protein

#### Probability of actual fold enrichment under the NULL distribution: How often specific domain fold enrichment among interactors of a specific viral protein (attribute of a viral_protein-human_domain pair) is observed as compared to any domain fold enrichment among interactors of that viral protein (attribute of a viral_protein, NULL distribution), the latter is derived from permutation mimicking as if viral protein was binding to a different set of human proteins   
Using fold enrichent as a statistic accounts for domain frequency in the background distribution. Background distribution of domain frequency in based on all human proteins with known PPI.   

Luckily, the network that we have is a bipartite undirected network (viral proteins - human proteins or human proteins - human domains). That means that we can shuffle interactions (links, edges) on the interface between viral and human proteins which is algorithmically quite fast even for large networks as compared to conventional unipartite undirected networks.  

One of the ways to describe the problem of finding domains that are likely to mediate interaction of human proteins with viral proteins is by drawing a network which has 3 types of elements: viral proteins (V1, V2, V3), human proteins (H1-10) and human domains (D1-10).

```r
knitr::include_graphics("./images/net_start.jpg")
```

<img src="./images/net_start.jpg" width="450px" />


First, we compute the fraction of human interacting partners of V2 that contain domain D6. Then, we compute the fold enrichment by dividing this fraction by the domain frequency among all proteins in the network.

```r
knitr::include_graphics("./images/net_start_calc.jpg")
```

<img src="./images/net_start_calc.jpg" width="450px" />

Next, we permute which human proteins interact with V1, V2 and V3, keeping the number of interaction per both viral and human proteins as well as number of edges (interactions) constant

```r
knitr::include_graphics("./images/net_scram.jpg")
```

<img src="./images/net_scram.jpg" width="450px" />

Finally, we repeat this procedure 1000-5000 times which gives us the probability of specific human domain - specific viral protein co-occurence if the domain and the protein were unrelated.  

We calculate this by comparing real fold enrichment value to the distribution of fold enrichment of randomly associated domains. For each human domain-viral protein value we calculate the fraction of the (random, NULL) distribution that is as high or higher than the real value.  


```r
knitr::include_graphics("./images/net_start_fraq.jpg")
```

<img src="./images/net_start_fraq.jpg" width="200px" />

```r
knitr::include_graphics("./images/net_scram_fraq.jpg")
```

<img src="./images/net_scram_fraq.jpg" width="200px" />


## Load network 

Network is saved by "Map_domains_to_network.Rmd" script which describes how the network and domain data were obtained.


```r
viral_human_net_w_domains_d = fread("./processed_data_files/viral_human_net_w_domains", sep = "\t", stringsAsFactors = F)

# generate minimal information tables
viral_human_net = unique(viral_human_net_w_domains_d[,.(IDs_interactor_viral, IDs_interactor_human, IDs_interactor_viral_degree)])

domains_proteins = unique(viral_human_net_w_domains_d[,.(IDs_interactor_human, IDs_domain_human, domain_frequency)])

viral_human_net_w_domains = unique(viral_human_net_w_domains_d[,.(IDs_interactor_viral, IDs_interactor_human, IDs_domain_human, domain_frequency_per_IDs_interactor_viral, fold_enrichment)])

# remove viral proteins with no corresponding human domains from this version of the data: we can't calculate the probability of viral protein - human domain interaction if there isn't any. However, permutations of viral-human interaction network are performed on all interactions including the one's when human partner doesn't have any domains, viral_human_net contains all interactions.
viral_human_net_w_domains = viral_human_net_w_domains[IDs_domain_human != "",]
```

## Calculate fold enrichment distribution of each viral protein in any domain and plot a few examples


```r
if(!file.exists("./large_processed_data_files/viral_foldEnrichDist.tsv.gz")){
    # calculate fold enrichment distribution
    viral_foldEnrichDist = foldEnrichmentDist(net = viral_human_net, 
                                              protein_annot = domains_proteins, 
                                              N = 1000, cores = cores_to_use, seed = 1, frequency = frequency)
    fwrite(viral_foldEnrichDist, "./large_processed_data_files/viral_foldEnrichDist.tsv", sep = "\t")
    gzip("./large_processed_data_files/viral_foldEnrichDist.tsv", remove = T)
}
gunzip("./large_processed_data_files/viral_foldEnrichDist.tsv.gz", remove = F)
viral_foldEnrichDist = fread("./large_processed_data_files/viral_foldEnrichDist.tsv", sep = "\t", stringsAsFactors = F)
```

```
## 
Read 57.8% of 19255176 rows
Read 19255176 rows and 2 (of 2) columns from 0.382 GB file in 00:00:03
```

```r
unlink("./large_processed_data_files/viral_foldEnrichDist.tsv")
dim(viral_foldEnrichDist)
```

```
## [1] 19255176        2
```

```r
# in how many cases fold enrichment or frequency is 0?
if(frequency) viral_foldEnrichDist[, mean(sampled_domain_frequency_per_set == 0)]
if(!frequency) viral_foldEnrichDist[, mean(sampled_fold_enrichment == 0)]
```

```
## [1] 0.2847763
```

```r
# plot a few random cases
set.seed(1)
plotFoldEnrichmentDist(proteinID = sample(unique(viral_human_net$IDs_interactor_viral), 9), 
                       fold_enrichment_dist = viral_foldEnrichDist, 
                       data = viral_human_net_w_domains, text_lab = F, frequency = frequency)
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](domain_enrichment_permutation_files/figure-html/viral_foldEnrichDist-1.png)<!-- -->

## Calculate P value for each viral protein and each human domain association


```r
if(!file.exists("./processed_data_files/viralProtein_humanDomain_pval.tsv")){
    proctime = proc.time()
    # calculate pvalue
    Pvals = foldEnrichmentPval(fold_enrichment_dist = viral_foldEnrichDist, 
                               data = viral_human_net_w_domains, cores = cores_to_use, frequency = frequency)
    proctime = proc.time() - proctime
    fwrite(Pvals, "./processed_data_files/viralProtein_humanDomain_pval.tsv", sep = "\t")
    proctime
}
Pvals = fread("./processed_data_files/viralProtein_humanDomain_pval.tsv", sep = "\t", stringsAsFactors = F)

# plot pvalue distribution
ggplot(Pvals, aes(x = Pval)) + geom_histogram(bins = 100) + ggtitle("viral protein and human domain association \n pvalue distribution") + theme_light() + xlim(0,1)
```

![](domain_enrichment_permutation_files/figure-html/Pvals-1.png)<!-- -->

```r
# qvalue adjustment procedure fails because of weird p-value distribution
if(frequency) Pvals_nozeros = foldEnrichmentPval(fold_enrichment_dist = viral_foldEnrichDist[sampled_domain_frequency_per_set != 0,], 
                               data = viral_human_net_w_domains, cores = cores_to_use, frequency = frequency)
if(!frequency) Pvals_nozeros = foldEnrichmentPval(fold_enrichment_dist = viral_foldEnrichDist[sampled_fold_enrichment != 0,], 
                               data = viral_human_net_w_domains, cores = cores_to_use, frequency = frequency)
# plot pvalue distribution - no zero domain proteins
ggplot(Pvals_nozeros, aes(x = Pval)) + geom_histogram(bins = 100) + ggtitle("viral protein and human domain association \n pvalue distribution, zero domain proteins discarded") + theme_light() + xlim(0,1)
```

![](domain_enrichment_permutation_files/figure-html/Pvals-2.png)<!-- -->

```r
# discard zero domain proteins for downstream analyses
Pvals = Pvals_nozeros
rm(viral_foldEnrichDist)

# calculate fdr adjusted pvalue
Pvals[, Pval_fdr := p.adjust(Pval, method = "fdr")]
ggplot(Pvals, aes(x = Pval_fdr)) + geom_histogram(bins = 100) + ggtitle("viral protein and human domain association \n FDR adjusted pvalue distribution") + theme_light() + xlim(0,1)
```

![](domain_enrichment_permutation_files/figure-html/Pvals-3.png)<!-- -->

```r
# calculate q-value
Pvals_qvalue = Pvals[, qvalue(Pval)]
Pvals[, Qval := qvalue(Pval)$qvalues]
ggplot((Pvals), aes(x = Qval)) + geom_histogram(bins = 100) + ggtitle("viral protein and human domain association \n qvalue adjusted pvalue distribution") + theme_light() + xlim(0,1)
```

![](domain_enrichment_permutation_files/figure-html/Pvals-4.png)<!-- -->

```r
# plot fold enrichment (or domain_frequency_per_IDs_interactor_viral) vs pvalue
if(frequency) ggplot(Pvals, aes(x = domain_frequency_per_IDs_interactor_viral, y = Pval)) + geom_bin2d() + ggtitle("viral protein and human domain association \n domain frequency per viral protein vs p-value") + ylab("p-value") + xlab("fold enrichment") + scale_x_log10()
if(!frequency) ggplot(Pvals, aes(x = fold_enrichment, y = Pval)) + geom_bin2d() + ggtitle("viral protein and human domain association \n fold enrichment vs p-value") + ylab("p-value") + xlab("fold enrichment") + scale_x_log10()
```

![](domain_enrichment_permutation_files/figure-html/Pvals-5.png)<!-- -->

## Relationships between pvalue statistic and protein/domain properties


```r
# merge results to original data
if(frequency) viral_human_net_w_domains_d2 = viral_human_net_w_domains_d[Pvals, on = c("IDs_interactor_viral", "IDs_domain_human", "domain_frequency_per_IDs_interactor_viral")]
if(!frequency) viral_human_net_w_domains_d2 = viral_human_net_w_domains_d[Pvals, on = c("IDs_interactor_viral", "IDs_domain_human", "fold_enrichment")]

# function to accomodate ggplot2::geom_bin2d in GGally::ggpairs, taken from http://ggobi.github.io/ggally/#custom_functions
d2_bin <- function(data, mapping, ..., low = "#132B43", high = "#56B1F7") {
    ggplot(data = data, mapping = mapping) +
        geom_bin2d(...) +
        scale_fill_gradient(low = low, high = high) +
        scale_y_log10() + scale_x_log10()
}

log10_density = function(data, mapping, ...){
    ggplot(data = data, mapping = mapping) +
        geom_density(...) +
        scale_x_log10()
}

GGally::ggpairs(viral_human_net_w_domains_d2[IDs_domain_human != "",.(domain_count, domain_frequency, 
                                               Taxid_interactor_viral,
                                               IDs_interactor_viral_degree, 
                                               IDs_interactor_human_degree, 
                                               domain_count_per_IDs_interactor_viral,
                                               domain_frequency_per_IDs_interactor_viral,
                                               fold_enrichment,
                                               Pval, Pval_fdr, Qval)], 
                lower = list(continuous = d2_bin), 
                diag = list(continuous = log10_density)) +
    theme_light() +
    theme(strip.text.y = element_text(angle = 0, size = 10),
          strip.text.x = element_text(angle = 90, size = 10))
```

![](domain_enrichment_permutation_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

```r
hist(viral_human_net_w_domains_d2[Pval_fdr < 0.35, domain_count_per_IDs_interactor_viral], main = "number of human proteins with specific domain per viral protein \n the lowest FDR corrected peak (0.35 < pval)", xlab = "how many times domain is repeated")
hist(viral_human_net_w_domains_d2[Pval_fdr < 0.35, domain_count], main = "domain prevalence in the background set \n the lowest FDR corrected peak (0.35 < pval)", xlab = "domain count in background")
```

## Try deleting domains with low background counts


```r
if(!file.exists("./large_processed_data_files/minus_low_background_counts")){
    # recalculate using all domains
    d0 = domainEnrichment(backgr_domain_count = 0, 
                          net = viral_human_net,
                          protein_annot = domains_proteins, 
                          data = viral_human_net_w_domains,
                          N = 500, cores = cores_to_use, seed = 1, 
                          all.data = viral_human_net_w_domains_d, 
                          frequency = frequency, pValPlot = "Pval",
                          remove_zeros = F, calculate_qval = F)
    
    # delete domains with background count <= 1
    d1 = domainEnrichment(backgr_domain_count = 1, 
                          net = viral_human_net,
                          protein_annot = domains_proteins, 
                          data = viral_human_net_w_domains,
                          N = 500, cores = cores_to_use, seed = 1, 
                          all.data = viral_human_net_w_domains_d, 
                          frequency = frequency, pValPlot = "Pval",
                          remove_zeros = F, calculate_qval = F)
    
    # delete domains with background count <= 2
    d2 = domainEnrichment(backgr_domain_count = 2, 
                          net = viral_human_net,
                          protein_annot = domains_proteins, 
                          data = viral_human_net_w_domains,
                          N = 500, cores = cores_to_use, seed = 1, 
                          all.data = viral_human_net_w_domains_d, 
                          frequency = frequency, pValPlot = "Pval",
                          remove_zeros = F, calculate_qval = F)

    # delete domain with background count <= 4
    d4 = domainEnrichment(backgr_domain_count = 4, 
                          net = viral_human_net,
                          protein_annot = domains_proteins, 
                          data = viral_human_net_w_domains,
                          N = 500, cores = cores_to_use, seed = 1, 
                          all.data = viral_human_net_w_domains_d, 
                          frequency = frequency, pValPlot = "Pval",
                          remove_zeros = F, calculate_qval = F)

    # delete domains with background count <= 8
    d8 = domainEnrichment(backgr_domain_count = 8, 
                          net = viral_human_net,
                          protein_annot = domains_proteins, 
                          data = viral_human_net_w_domains,
                          N = 500, cores = cores_to_use, seed = 1, 
                          all.data = viral_human_net_w_domains_d, 
                          frequency = frequency, pValPlot = "Pval",
                          remove_zeros = F, calculate_qval = F)

    # delete domains with background count <= 16
    d16 = domainEnrichment(backgr_domain_count = 16, 
                           net = viral_human_net,
                           protein_annot = domains_proteins, 
                           data = viral_human_net_w_domains,
                           N = 500, cores = cores_to_use, seed = 1, 
                           all.data = viral_human_net_w_domains_d, 
                          frequency = frequency, pValPlot = "Pval",
                          remove_zeros = F, calculate_qval = F)
    
    save(d0, d1, d2, d4, d8, d16, file = "./large_processed_data_files/minus_low_background_counts")
}
```


## Deleting domains with low background frequency counts leads to ...


```r
load(file = "./large_processed_data_files/minus_low_background_counts")
# show p-value distributions and sampled fold_enrichment or frequency distributions
d0$PvalPlot
```

![](domain_enrichment_permutation_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

```r
d0$distPlot
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](domain_enrichment_permutation_files/figure-html/unnamed-chunk-5-2.png)<!-- -->

```r
d1$PvalPlot
```

![](domain_enrichment_permutation_files/figure-html/unnamed-chunk-5-3.png)<!-- -->

```r
d1$distPlot
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](domain_enrichment_permutation_files/figure-html/unnamed-chunk-5-4.png)<!-- -->

```r
d2$PvalPlot
```

![](domain_enrichment_permutation_files/figure-html/unnamed-chunk-5-5.png)<!-- -->

```r
d2$distPlot
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](domain_enrichment_permutation_files/figure-html/unnamed-chunk-5-6.png)<!-- -->

```r
d4$PvalPlot
```

![](domain_enrichment_permutation_files/figure-html/unnamed-chunk-5-7.png)<!-- -->

```r
d4$distPlot
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](domain_enrichment_permutation_files/figure-html/unnamed-chunk-5-8.png)<!-- -->

```r
d8$PvalPlot
```

![](domain_enrichment_permutation_files/figure-html/unnamed-chunk-5-9.png)<!-- -->

```r
d8$distPlot
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](domain_enrichment_permutation_files/figure-html/unnamed-chunk-5-10.png)<!-- -->

```r
d16$PvalPlot
```

![](domain_enrichment_permutation_files/figure-html/unnamed-chunk-5-11.png)<!-- -->

```r
d16$distPlot
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![](domain_enrichment_permutation_files/figure-html/unnamed-chunk-5-12.png)<!-- -->

## R session information


```r
Sys.Date()
```

```
## [1] "2017-08-01"
```

```r
sessionInfo()
```

```
## R version 3.4.0 (2017-04-21)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS Sierra 10.12.6
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] GGally_1.3.1        qvalue_2.8.0        MItools_0.1.7      
##  [4] Biostrings_2.44.1   XVector_0.16.0      PSICQUIC_1.14.0    
##  [7] plyr_1.8.4          httr_1.2.1          biomaRt_2.32.1     
## [10] IRanges_2.10.2      S4Vectors_0.14.3    BiocGenerics_0.22.0
## [13] ggplot2_2.2.1       R.utils_2.5.0       R.oo_1.21.0        
## [16] R.methodsS3_1.7.1   data.table_1.10.4  
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.11               lattice_0.20-35           
##  [3] Rsamtools_1.28.0           rprojroot_1.2             
##  [5] digest_0.6.12              R6_2.2.2                  
##  [7] GenomeInfoDb_1.12.2        backports_1.1.0           
##  [9] RSQLite_2.0                evaluate_0.10             
## [11] zlibbioc_1.22.0            rlang_0.1.1               
## [13] lazyeval_0.2.0             blob_1.1.0                
## [15] Matrix_1.2-10              rmarkdown_1.6             
## [17] gsubfn_0.6-6               labeling_0.3              
## [19] proto_1.0.0                splines_3.4.0             
## [21] BiocParallel_1.10.1        stringr_1.2.0             
## [23] RCurl_1.95-4.8             bit_1.1-12                
## [25] munsell_0.4.3              DelayedArray_0.2.7        
## [27] compiler_3.4.0             rtracklayer_1.36.3        
## [29] htmltools_0.3.6            SummarizedExperiment_1.6.3
## [31] tibble_1.3.3               GenomeInfoDbData_0.99.0   
## [33] matrixStats_0.52.2         XML_3.98-1.9              
## [35] reshape_0.8.6              GenomicAlignments_1.12.1  
## [37] bitops_1.0-6               grid_3.4.0                
## [39] gtable_0.2.0               DBI_0.7                   
## [41] magrittr_1.5               scales_0.4.1              
## [43] stringi_1.1.5              reshape2_1.4.2            
## [45] RColorBrewer_1.1-2         tools_3.4.0               
## [47] bit64_0.9-7                Biobase_2.36.2            
## [49] yaml_2.1.14                AnnotationDbi_1.38.1      
## [51] colorspace_1.3-2           GenomicRanges_1.28.3      
## [53] memoise_1.1.0              knitr_1.16
```
