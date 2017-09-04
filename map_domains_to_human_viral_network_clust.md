# Mapping domains to the human-viral protein interaction network
Vitalii Kleshchevnikov  
29/06/2017  



Date: 2017-09-03 23:08:06

## Map domain information on the protein interaction network

I read interaction data and clean this data to make it more useble. Then, I filter and keep only human-viral interactions.


```r
# load protein-domain mapping
protein_domain_pair = fread(file = "./processed_data_files/protein_domain_pair", stringsAsFactors = F, sep = "\t")
# load PPI data
all_viral_interaction = queryPSICQUICrlib(query = "species:10239",
                format = "tab27",
                database = "imex",
                directory = "./data_files/", 
                releaseORdate = "2017Jul12"
                # ,just_list_releases = T,
                )
```

```
## found local copy of the data from the latest release ... reading into R
```

```r
all_viral_interaction = cleanMITAB(all_viral_interaction)
# pick proteins only
all_viral_interaction$data = all_viral_interaction$data[interactor_IDs_databases_A == "uniprotkb" & interactor_IDs_databases_B == "uniprotkb", ]
# filter only human-viral interactions
all_viral_interaction$data = all_viral_interaction$data[Taxid_interactor_A == "9606" | Taxid_interactor_B == "9606",]
all_viral_interaction$data = unique(all_viral_interaction$data)
```

Both the network and the domain data contain more information than necessary for identifying domains likely to mediate interaction. Selecting only what's necessary: viral_protein_UniprotID - human_protein_UniprotID - human_domain_InterProID.  

First, we need to rearrange interactions so that one column contains viral proteins and the other contains human proteins (in the database some interactions are stored viral-human order and some in human-viral order).   


```r
all_viral_interaction_ab = all_viral_interaction$data[Taxid_interactor_A == "9606" & Taxid_interactor_B != "9606",
                                                 .(IDs_interactor_human = IDs_interactor_A, 
                                                   IDs_interactor_viral = IDs_interactor_B, 
                                                   Taxid_interactor_human = Taxid_interactor_A, 
                                                   Taxid_interactor_viral = Taxid_interactor_B)]
all_viral_interaction_ba = all_viral_interaction$data[Taxid_interactor_A != "9606" & Taxid_interactor_B == "9606",
                                                 .(IDs_interactor_human = IDs_interactor_B,
                                                   IDs_interactor_viral = IDs_interactor_A,
                                                   Taxid_interactor_human = Taxid_interactor_B,
                                                   Taxid_interactor_viral = Taxid_interactor_A)]
all_viral_interaction_simp = unique(rbind(all_viral_interaction_ab, all_viral_interaction_ba))
```

Next, the frequency of each domain is calculated and simplified domain information (IDs_protein, IDs_domain, domain_type) is integrated with the simplified network data described above.  


```r
protein_domain_pair_temp = copy(protein_domain_pair)[, IDs_interactor_human := IDs_protein][, IDs_protein := NULL][, IDs_domain_human := IDs_domain][, IDs_domain := NULL]

# calculate total number of human proteins
protein_domain_pair_temp[, N_prot_w_interactors := length(unique(IDs_interactor_human))]
# calculate domain count and frequency
protein_domain_pair_temp[, domain_count := length(unique(IDs_interactor_human)), by = IDs_domain_human]
protein_domain_pair_temp[, domain_frequency := domain_count / N_prot_w_interactors]

# calculate network descriptive stats
# viral protein degree
all_viral_interaction_simp[, IDs_interactor_viral_degree := length(unique(IDs_interactor_human)), by = IDs_interactor_viral]
# human protein degree
all_viral_interaction_simp[, IDs_interactor_human_degree := length(unique(IDs_interactor_viral)), by = IDs_interactor_human]

# I keep all interactions even if a human protein has no known domain (all.x = F, all.y = T)
viral_human_w_domains = merge(protein_domain_pair_temp, all_viral_interaction_simp, all.x = F, all.y = T, by = "IDs_interactor_human", allow.cartesian = T)

# human domains per viral protein
viral_human_w_domains[, IDs_domain_human_per_IDs_interactor_viral := length(unique(IDs_domain_human)), by = IDs_interactor_viral]
# viral protein per human domain
viral_human_w_domains[, IDs_interactor_viral_per_IDs_domain_human := length(unique(IDs_interactor_viral)), by = IDs_domain_human]

# domain count but per viral protein human domain instances (how many proteins the domain is located in) per viral protein (ID) and human domain (ID)
viral_human_w_domains[, domain_count_per_IDs_interactor_viral := length(unique(IDs_interactor_human)), by = .(IDs_interactor_viral, IDs_domain_human)]
viral_human_w_domains[is.na(IDs_domain_human), domain_count_per_IDs_interactor_viral := 0]
# domain frequency but per viral protein
viral_human_w_domains[, domain_frequency_per_IDs_interactor_viral := domain_count_per_IDs_interactor_viral / IDs_interactor_viral_degree, by = IDs_interactor_viral]
viral_human_w_domains[is.na(IDs_domain_human), domain_frequency_per_IDs_interactor_viral := 0]
# fold enrichment
viral_human_w_domains[, fold_enrichment := domain_frequency_per_IDs_interactor_viral / domain_frequency]
viral_human_w_domains[is.na(IDs_domain_human), fold_enrichment := 0]


# save resulting network
fwrite(viral_human_w_domains, file = "./processed_data_files/viral_human_net_w_domains", sep = "\t")
```

## Summary of the network

The plot below show the relationships (2D histogram), distribution density (on the diagonal) and pearson correlation for a number of parameters characterising human domains, viral proteins or it's interactions:

1. domain_frequency is the number of human proteins with a particular domain divided by the total number of human proteins (the attribute of a human domain)
2. IDs_interactor_viral_degree is the number of interactions each viral protein has (the attribute of a viral protein)
3. IDs_interactor_human_degree is the number of interactions each human protein has (the attribute of a human protein)
4. IDs_domain_human_per_IDs_interactor_viral is the number of domain types that are present in proteins which a particular viral protein interacts with (the attribute of a viral protein)


```r
# function to accomodate ggplot2::geom_bin2d in GGally::ggpairs, taken from http://ggobi.github.io/ggally/#custom_functions
d2_bin <- function(data, mapping, ..., low = "#132B43", high = "#56B1F7") {
    ggplot(data = data, mapping = mapping) +
        geom_bin2d(...) +
        scale_fill_gradient(low = low, high = high) +
        scale_y_log10() + scale_x_log10() + annotation_logticks()
}

log10_density = function(data, mapping, ...){
    ggplot(data = data, mapping = mapping) +
        geom_density(...) +
        scale_x_log10() + annotation_logticks()
}

d2_bin_plot = GGally::ggpairs(viral_human_w_domains[!is.na(IDs_domain_human),.(domain_frequency, 
                                                                 IDs_interactor_viral_degree, 
                                                                 IDs_interactor_human_degree, 
                                                                 IDs_domain_human_per_IDs_interactor_viral, 
                                                                 IDs_interactor_viral_per_IDs_domain_human,
                                                                 domain_count_per_IDs_interactor_viral,
                                                                 domain_frequency_per_IDs_interactor_viral,
                                                                 fold_enrichment)], 
                lower = list(continuous = d2_bin), 
                diag = list(continuous = log10_density)) +
    theme_light() +
    theme(strip.text.y = element_text(angle = 0, size = 10),
          strip.text.x = element_text(angle = 90, size = 10))
d2_bin_plot
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/map_domains_to_human_viral_network_clust_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

## R session information


```r
save(list = ls(), file="./processed_data_files/map_domains_to_human_viral_network_clust.RData")
R.utils::gzip(filename = "./processed_data_files/map_domains_to_human_viral_network_clust.RData",
              destname = "./processed_data_files/map_domains_to_human_viral_network_clust.RData.gz",
              remove = T, overwrite = T)

Sys.Date()
```

```
## [1] "2017-09-03"
```

```r
sessionInfo()
```

```
## R version 3.4.1 (2017-06-30)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Red Hat Enterprise Linux Server 7.3 (Maipo)
## 
## Matrix products: default
## BLAS: /hps/nobackup/research/petsalaki/users/vitalii/R/lib64/R/lib/libRblas.so
## LAPACK: /hps/nobackup/research/petsalaki/users/vitalii/R/lib64/R/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    parallel  methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] GGally_1.3.2         MItools_0.1.15       Biostrings_2.44.1   
##  [4] XVector_0.16.0       rtracklayer_1.36.3   GenomicRanges_1.28.3
##  [7] GenomeInfoDb_1.12.2  ggplot2_2.2.1        PSICQUIC_1.14.0     
## [10] plyr_1.8.4           httr_1.3.1           biomaRt_2.32.1      
## [13] IRanges_2.10.2       S4Vectors_0.14.3     UniProt.ws_2.16.0   
## [16] BiocGenerics_0.22.0  RCurl_1.95-4.8       bitops_1.0-6        
## [19] RSQLite_2.0          R.utils_2.5.0        R.oo_1.21.0         
## [22] R.methodsS3_1.7.1    downloader_0.4       data.table_1.10.4   
## [25] rmarkdown_1.6       
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.12               lattice_0.20-35           
##  [3] Rsamtools_1.28.0           rprojroot_1.2             
##  [5] digest_0.6.12              R6_2.2.2                  
##  [7] backports_1.1.0            evaluate_0.10.1           
##  [9] zlibbioc_1.22.0            rlang_0.1.2               
## [11] lazyeval_0.2.0             blob_1.1.0                
## [13] Matrix_1.2-11              qvalue_2.8.0              
## [15] gsubfn_0.6-6               labeling_0.3              
## [17] proto_1.0.0                splines_3.4.1             
## [19] BiocParallel_1.10.1        stringr_1.2.0             
## [21] bit_1.1-12                 munsell_0.4.3             
## [23] DelayedArray_0.2.7         compiler_3.4.1            
## [25] htmltools_0.3.6            SummarizedExperiment_1.6.3
## [27] tibble_1.3.4               GenomeInfoDbData_0.99.0   
## [29] matrixStats_0.52.2         XML_3.98-1.9              
## [31] reshape_0.8.7              GenomicAlignments_1.12.1  
## [33] grid_3.4.1                 gtable_0.2.0              
## [35] DBI_0.7                    magrittr_1.5              
## [37] scales_0.5.0               stringi_1.1.5             
## [39] reshape2_1.4.2             RColorBrewer_1.1-2        
## [41] tools_3.4.1                bit64_0.9-7               
## [43] Biobase_2.36.2             yaml_2.1.14               
## [45] AnnotationDbi_1.38.1       colorspace_1.3-2          
## [47] memoise_1.1.0              knitr_1.17
```
