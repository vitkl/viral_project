# remove_redundant_domains
Vitalii Kleshchevnikov  
3/09/2017  



Date: 2017-10-22 16:53:39

## Read InterProScan result and filter for "Domain", "Active_site", "Binding_site", "Conserved_site", "PTM" signatures

I read InterProScan result and the InterPro_entry_types file. InterProScan output format description: https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats  


```r
# read InterProScan result, download and add InterPro Entry Types information, extract from relevant columns and add names metadata and sequence length information
InterProScan_result = readInterProGFF3("./processed_data_files/all_human_viral_protein_domains102017.gff3.gz", processed = F)
InterProScan_result = addInterProEntryTypes(InterProScan_result, "./data_files/entry.list")
InterProScan_domains = SubsetByInterProEntryType(InterProScan_result, c("Domain"))
```

## Remove redundancy in the identified domains

Signatures of the same domain from different InterPro member databases.  


```r
#### collapse domains by family (not only by single InterProID removing multiple the databases), In domain families: assign all domains to the top domain.
#InterProTree = loadInterProTree(filename = "./data_files/interpro_ParentChildTreeFile.txt")
# find all InterProIDs that have children InterProIDs
#InterProTree_level1 = getLevelXchildren(tree = InterProTree, level = 1) 
#has_children = InterProScan_domains$Dbxref %in% InterProTree_level1$allchildren

# map top InterProIDs in the hierarchy when InterProIDs have have children InterProIDs
#InterProScan_domains$topInterProID = InterProScan_domains$Dbxref
#InterProScan_domains$topInterProID[has_children] = InterProTree_level1$level1[match(InterProScan_domains$Dbxref, #InterProTree_level1$allchildren)][has_children]
#### 
InterProScan_domains = unique(InterProScan_domains)
# collapse on InterProID contained in Dbxref column
InterProScan_domains_nonred = collapseByInterProID(InterProScan_features = InterProScan_domains, id_col = "Dbxref")
```

We used to have 7221 domains and 86143 protein-domain pairs.   
After combining the signatures of the same domain from different InterPro member databases we have 4761 domains and 47348 protein-domain pairs.   
Finally, after considering only one domain per domain family we are left with 0 domains and 15940 protein-domain pairs.  


```r
# save results
export(InterProScan_domains_nonred, con = "./processed_data_files/InterProScan_domains_nonredundant.gff3", format = "gff3")
gzip(filename = "./processed_data_files/InterProScan_domains_nonredundant.gff3", destname = "./processed_data_files/InterProScan_domains_nonredundant.gff3.gz", overwrite = T, remove = T)

#InterProScan_domains_nonred2 = readInterProGFF3("./processed_data_files/InterProScan_domains_nonredundant.gff3.gz", processed = T)

# take protein-domain pair discarding range information
protein_domain_pair = unique(data.table(IDs_protein = as.character(seqnames(InterProScan_domains_nonred)),
                                        IDs_domain = as.character(InterProScan_domains_nonred$Dbxref),
                                        all_IDs_domain = as.character(InterProScan_domains_nonred$Name),
                                        domain_type = InterProScan_domains_nonred$ENTRY_TYPE))
# save simplified table
fwrite(protein_domain_pair, file = "./processed_data_files/protein_domain_pair", sep = "\t")
```

## R session information


```r
save(list = ls(), file="./processed_data_files/remove_redundant_domains_clust.RData")
R.utils::gzip(filename = "./processed_data_files/remove_redundant_domains_clust.RData",
              destname = "./processed_data_files/remove_redundant_domains_clust.RData.gz",
              remove = T, overwrite = T)

Sys.Date()
```

```
## [1] "2017-10-22"
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
##  [1] GGally_1.3.2         MItools_0.1.27       Biostrings_2.44.1   
##  [4] XVector_0.16.0       rtracklayer_1.36.3   GenomicRanges_1.28.3
##  [7] GenomeInfoDb_1.12.2  ggplot2_2.2.1        PSICQUIC_1.14.0     
## [10] plyr_1.8.4           httr_1.3.1           biomaRt_2.32.1      
## [13] IRanges_2.10.2       S4Vectors_0.14.3     UniProt.ws_2.16.0   
## [16] BiocGenerics_0.22.0  RCurl_1.95-4.8       bitops_1.0-6        
## [19] RSQLite_2.0          R.utils_2.5.0        R.oo_1.21.0         
## [22] R.methodsS3_1.7.1    downloader_0.4       data.table_1.10.4-2 
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
## [15] gsubfn_0.6-6               proto_1.0.0               
## [17] splines_3.4.1              BiocParallel_1.10.1       
## [19] stringr_1.2.0              bit_1.1-12                
## [21] munsell_0.4.3              DelayedArray_0.2.7        
## [23] compiler_3.4.1             htmltools_0.3.6           
## [25] SummarizedExperiment_1.6.3 tibble_1.3.4              
## [27] GenomeInfoDbData_0.99.0    matrixStats_0.52.2        
## [29] XML_3.98-1.9               reshape_0.8.7             
## [31] GenomicAlignments_1.12.1   grid_3.4.1                
## [33] ontologyIndex_2.4          jsonlite_1.5              
## [35] gtable_0.2.0               DBI_0.7                   
## [37] magrittr_1.5               scales_0.5.0              
## [39] stringi_1.1.5              reshape2_1.4.2            
## [41] RColorBrewer_1.1-2         tools_3.4.1               
## [43] bit64_0.9-7                Biobase_2.36.2            
## [45] yaml_2.1.14                AnnotationDbi_1.38.1      
## [47] colorspace_1.3-2           memoise_1.1.0             
## [49] knitr_1.17
```
