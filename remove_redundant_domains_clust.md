# remove_redundant_domains
Vitalii Kleshchevnikov  
3/09/2017  



Date: 2017-09-03 19:55:05

## Read InterProScan result and filter for "Domain", "Active_site", "Binding_site", "Conserved_site", "PTM" signatures

I read InterProScan result and the InterPro_entry_types file. InterProScan output format description: https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats  


```r
# read InterProScan result, download and add InterPro Entry Types information, extract from relevant columns and add names metadata and sequence length information
InterProScan_result = readInterProGFF3("./processed_data_files/all_human_viral_protein_domains.gff3.gz")
InterProScan_result = addInterProEntryTypes(InterProScan_result, "./data_files/entry.list")
# create a subset that contains "Domain", "Active_site", "Binding_site", "Conserved_site", "PTM" signatures
InterProScan_domains = SubsetByInterProEntryType(InterProScan_result)
```

## Remove redundancy in the identified domains (signatures of the same domain from different InterPro member databases)


```r
InterProScan_domains_list = split(InterProScan_domains, paste0(seqnames(InterProScan_domains),"|", InterProScan_domains$Dbxref))
InterProScan_domains_nonred = GRanges()

for (i in 1:length(InterProScan_domains_list)) {
    Nseq = length(InterProScan_domains_list[[i]])
    reduced = reduce(InterProScan_domains_list[[i]], with.revmap=T)
    metadata_orig = mcols(InterProScan_domains_list[[i]])
    
    metadata_new = metadata_orig[-which(colnames(metadata_orig) %in% c("status", "signature_desc", "md5", "ID", "type", "phase"))]
    metadata_new$source = (as.character(metadata_new$source))
    
    mapping = character(Nseq)
    for (ind in 1:length(reduced$revmap)) {
        number = reduced$revmap[[ind]]
        mapping[number] = ind
    }
    metadata_new_list = split(metadata_new, mapping)
    new_metadata_new_list = SplitDataFrameList()
    for (indic in 1:length(metadata_new_list)) {
        temp_metadata = metadata_new_list[[indic]]
        temp_metadata$source = paste0(temp_metadata$source, collapse="|")
        temp_metadata$score = paste0(temp_metadata$score, collapse="|")
        temp_metadata$Target = paste0(temp_metadata$Target, collapse="|")
        temp_metadata$Name = paste0(temp_metadata$Name, collapse="|")
        # convert CharacterList to character first 
        temp_metadata$Ontology_term = sapply(temp_metadata$Ontology_term, paste0, collapse=",")
        temp_metadata$Ontology_term = paste0(temp_metadata$Ontology_term, collapse="|")
        new_metadata_new_list[[indic]] = unique(temp_metadata)
    }
    metadata_new = do.call("rbind", new_metadata_new_list)
    mcols(reduced) = metadata_new
    InterProScan_domains_nonred = c(InterProScan_domains_nonred, reduced)
}

InterProScan_domains_nonred = unique(InterProScan_domains_nonred) 
```


```r
# save results
export(InterProScan_domains_nonred, con = "./processed_data_files/InterProScan_domains_nonredundant.gff3", format = "gff3")

# take protein-domain pair discarding range information
protein_domain_pair = unique(data.table(IDs_protein = as.character(seqnames(InterProScan_domains_nonred)), IDs_domain = as.character(InterProScan_domains_nonred$Dbxref), domain_type = InterProScan_domains_nonred$ENTRY_TYPE))
# save simplified table
fwrite(protein_domain_pair, file = "./processed_data_files/protein_domain_pair", sep = "\t")
```



```r
positionalDistance = function(InterProScan_domains, maxgap = 100){
    # generate Granges which contain feature start
    domain_start = resize(InterProScan_domains, width = 1, fix="start", use.names=TRUE)
    # generate Granges which contain feature end
    domain_end = resize(InterProScan_domains, width = 1, fix="end", use.names=TRUE)
    
    # find features with overlapping start (including feature with itself)
    overlap_start = findOverlaps(query = domain_start, maxgap = maxgap)
    # find features with overlapping end (including feature with itself)
    overlap_end = findOverlaps(query = domain_end, maxgap = maxgap)
    
    # remove overlapping start of feature with itself
    overlap_start = overlap_start[(queryHits(overlap_start) != subjectHits(overlap_start))]
    # remove overlapping end of feature with itself
    overlap_end = overlap_end[(queryHits(overlap_end) != subjectHits(overlap_end))]
    
    # find features in which both the start and the end overlap
    overlap_start_n_end = intersect(overlap_start, overlap_end)
    
    # calculate the distance between positions of overlapping start features
    overlap_both_start_dist = distance(domain_start[queryHits(overlap_start_n_end)], domain_start[subjectHits(overlap_start_n_end)])
    # calculate the distance between positions of overlapping end features
    overlap_both_end_dist = distance(domain_end[queryHits(overlap_start_n_end)], domain_end[subjectHits(overlap_start_n_end)])
    
    # sum these distances
    overlap_total_dist = overlap_both_start_dist + overlap_both_end_dist
    
    # map distances to domains
    
}
# generate Granges which contain feature start
domain_start = resize(InterProScan_domains, width = 1, fix="start", use.names=TRUE)
# generate Granges which contain feature end
domain_end = resize(InterProScan_domains, width = 1, fix="end", use.names=TRUE)

# find features with overlapping start (including feature with itself)
overlap_start = findOverlaps(query = domain_start, maxgap = 100)
# find features with overlapping end (including feature with itself)
overlap_end = findOverlaps(query = domain_end, maxgap = 100)

# remove overlapping start of feature with itself
overlap_start = overlap_start[(queryHits(overlap_start) != subjectHits(overlap_start))]
# remove overlapping end of feature with itself
overlap_end = overlap_end[(queryHits(overlap_end) != subjectHits(overlap_end))]

# find features in which both the start and the end overlap
overlap_start_n_end = intersect(overlap_start, overlap_end)

# calculate the distance between positions of overlapping start features
overlap_both_start_dist = distance(domain_start[queryHits(overlap_start_n_end)], domain_start[subjectHits(overlap_start_n_end)])
hist(overlap_both_start_dist, breaks = seq(0,100,1))
# calculate the distance between positions of overlapping end features
overlap_both_end_dist = distance(domain_end[queryHits(overlap_start_n_end)], domain_end[subjectHits(overlap_start_n_end)])
hist(overlap_both_end_dist, breaks = seq(0,100,1))

# sum these distances
overlap_total_dist = overlap_both_start_dist + overlap_both_end_dist
overlap_total_dist = data.table(overlap_total_dist = overlap_total_dist,
                                queryHits = InterProScan_domains[queryHits(overlap_start_n_end)]$Dbxref,
                                subjectHits = InterProScan_domains[subjectHits(overlap_start_n_end)]$Dbxref)
hist(overlap_total_dist, breaks = seq(0,200,1))
abline(v=10)

# select distance difference cutoff
overlap_start_n_end = overlap_start_n_end[overlap_total_dist < 10]
# generate non-redundant domain annotatations by keeping only the first domain signature among overlapping signatures
overlap_start_n_end_nonred = overlap_start_n_end[queryHits(overlap_start_n_end) < subjectHits(overlap_start_n_end)]
InterProScan_domains_nonred = InterProScan_domains[unique(queryHits(overlap_start_n_end_nonred))]
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
## [15] gsubfn_0.6-6               proto_1.0.0               
## [17] splines_3.4.1              BiocParallel_1.10.1       
## [19] stringr_1.2.0              bit_1.1-12                
## [21] munsell_0.4.3              DelayedArray_0.2.7        
## [23] compiler_3.4.1             htmltools_0.3.6           
## [25] SummarizedExperiment_1.6.3 tibble_1.3.4              
## [27] GenomeInfoDbData_0.99.0    matrixStats_0.52.2        
## [29] XML_3.98-1.9               reshape_0.8.7             
## [31] GenomicAlignments_1.12.1   grid_3.4.1                
## [33] gtable_0.2.0               DBI_0.7                   
## [35] magrittr_1.5               scales_0.5.0              
## [37] stringi_1.1.5              reshape2_1.4.2            
## [39] RColorBrewer_1.1-2         tools_3.4.1               
## [41] bit64_0.9-7                Biobase_2.36.2            
## [43] yaml_2.1.14                AnnotationDbi_1.38.1      
## [45] colorspace_1.3-2           memoise_1.1.0             
## [47] knitr_1.17
```
