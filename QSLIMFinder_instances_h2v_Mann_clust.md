# QSLIMFinder: domain instances, human to search (Mann data only), viral to filter
Vitalii Kleshchevnikov  
04/09/2017  



## use all the interacting partners of human proteins with significant domains to seach for motifs and their viral interactors to filter

QSLIMfinder options: 
- seqin - all the interacting partners of human proteins with significant domains (including viral proteins); FASTA sequences
- query - the viral interacting partners of human proteins with significant domains; names of the FASTA sequences in seqin
- resdir - directory where to store all results
- resfile - file where to write
- dismask - mask disordered regions using IUPRED
- consmask - relative conservation masking
- probcut - FDR-corrected probability cutoff
- qregion - Mask all but the region of the query from (and including) residue X to residue Y [0,-1] (useful for regions sufficient/necessary to interact)

## load PPI data


```r
# load ppi data
# human-human
all_human_interaction = fullInteractome(taxid = 9606, database = "IntActFTP", format = "tab27",
                                        clean = TRUE, protein_only = TRUE,
                                        directory = "./data_files/")
```

```
## ... looking for the date of the latest IntAct release ...
```

```
## Warning in download.file(input, tt, method = method, mode = "wb", quiet = !
## showProgress): downloaded length 591 != reported length 0
```

```
## ... looking for the date of the latest IntAct release ...
```

```
## Warning in download.file(input, tt, method = method, mode = "wb", quiet = !
## showProgress): downloaded length 591 != reported length 0
```

```
## ... loading local copy ...
```

```
## 
Read 0.0% of 785947 rows
Read 1.3% of 785947 rows
Read 10.2% of 785947 rows
Read 19.1% of 785947 rows
Read 28.0% of 785947 rows
Read 35.6% of 785947 rows
Read 43.3% of 785947 rows
Read 50.9% of 785947 rows
Read 58.5% of 785947 rows
Read 67.4% of 785947 rows
Read 76.3% of 785947 rows
Read 85.2% of 785947 rows
Read 87.8% of 785947 rows
Read 98.0% of 785947 rows
Read 785947 rows and 42 (of 42) columns from 2.940 GB file in 00:00:23
```

```r
all_human_interaction = subsetMITABbyPMIDs(MITABdata = all_human_interaction,
                                           PMIDs = c("26496610"))
# human-viral
all_viral_interaction = interSpeciesInteractome(taxid1 = 9606, taxid2 = 10239, database = "IntActFTP", format = "tab27",
                                                clean = TRUE, protein_only = TRUE, 
                                                directory = "./data_files/")
```

```
## ... looking for the date of the latest IntAct release ...
```

```
## Warning in download.file(input, tt, method = method, mode = "wb", quiet = !
## showProgress): downloaded length 591 != reported length 0
```

```
## ... looking for the date of the latest IntAct release ...
```

```
## Warning in download.file(input, tt, method = method, mode = "wb", quiet = !
## showProgress): downloaded length 591 != reported length 0
```

```
## ... loading local copy ...
```

```
## 
Read 0.0% of 785947 rows
Read 8.9% of 785947 rows
Read 19.1% of 785947 rows
Read 25.4% of 785947 rows
Read 34.4% of 785947 rows
Read 43.3% of 785947 rows
Read 50.9% of 785947 rows
Read 58.5% of 785947 rows
Read 67.4% of 785947 rows
Read 76.3% of 785947 rows
Read 86.5% of 785947 rows
Read 96.7% of 785947 rows
Read 785947 rows and 42 (of 42) columns from 2.940 GB file in 00:00:18
```

## load results of which domains are likely to mediate interaction

```r
# load the domain analysis results
load("./processed_data_files/what_we_find_VS_ELM_clust20171019.RData")
domain_res = res_count
rm(list = ls()[!ls() %in% c("domain_res", "all_human_interaction", "all_viral_interaction")])

# choose pvalue cutoff:
plot(domain_res)
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/QSLIMFinder_instances_h2v_Mann_clust_files/figure-html/load_results-1.png)<!-- -->

```r
proteins_w_signif_domains = unique(domain_res$data_with_pval[p.value <= 1, IDs_interactor_human])
```

### prepare sets of protein sequences and instructions for SLIMFinder


```r
# load FASTA
all.fasta = readAAStringSet(filepath = "./data_files/all_human_viral_proteins.fasta", format = "fasta")

# remove interactions if one of the interactors does not have a FASTA sequence
all_human_interaction = removeInteractionNoFASTA(all_human_interaction, all.fasta)
all_viral_interaction = removeInteractionNoFASTA(all_viral_interaction, all.fasta)

#all_viral_interaction$data = all_viral_interaction$data[IDs_interactor_A == "O60503",]
#forSLIMFinder = listInteractionSubsetFASTA(interaction_set1 = all_human_interaction, 
#                                           interaction_set2 = all_viral_interaction, 
#                                           seed_id_vect = "O60503",
#                                           fasta = all.fasta,
#                                           single_interact_from_set2 = T, set1_only = T)

forSLIMFinder = listInteractionSubsetFASTA(interaction_set1 = all_human_interaction, 
                                           interaction_set2 = all_viral_interaction, 
                                           seed_id_vect = proteins_w_signif_domains,
                                           fasta = all.fasta,
                                           single_interact_from_set2 = T, set1_only = F)

forSLIMFinder_Ready = filterInteractionSubsetFASTA_list(forSLIMFinder,  length_set1_min = 2, length_set2_min = 1)
# how many viral-human pairs in which viral proteins interact with human protein (with relevant domain) but are not associated with the domain
sum(!domainProteinPairMatch(forSLIMFinder_Ready, domain_res, remove = F))
```

```
## [1] 0
```

```r
# how many viral-human pairs in which viral proteins interact with human protein (with relevant domain) and are associated with the domain
sum(domainProteinPairMatch(forSLIMFinder_Ready, domain_res, remove = F))
```

```
## [1] 10152
```

```r
# total
forSLIMFinder_Ready$length
```

```
## [1] 10152
```

```r
#filter
domain_filt = domain_res
domain_filt$data_with_pval = domain_filt$data_with_pval[p.value <= 1,]
forSLIMFinder_Ready = domainProteinPairMatch(forSLIMFinder_Ready, domain_filt, remove = T)

### sanity check
obs = unique(domain_filt$data_with_pval[, .(IDs_interactor_human, IDs_interactor_viral)])
setkey(obs)

res = tstrsplit(gsub("interactors_of\\.|\\.","", names(forSLIMFinder_Ready$fasta_subset_list)), ":")
res = data.table(idA = res[[1]], idB = res[[2]])
res_unkeyed = copy(res)
res = unique(res[complete.cases(res),])
setkey(res)

mean(paste0(res$idA, res$idB) %in% paste0(obs$IDs_interactor_human, obs$IDs_interactor_viral))
```

```
## [1] 1
```

```r
query_in_FASTA = sapply(1:length(forSLIMFinder_Ready$fasta_subset_list), function(i, forSLIMFinder){
    res = tstrsplit(gsub("interactors_of\\.|\\.","", names(forSLIMFinder$fasta_subset_list)), ":")
    res[[2]][i] %in% names(forSLIMFinder$fasta_subset_list[[i]])
}, forSLIMFinder_Ready)
mean(query_in_FASTA)
```

```
## [1] 1
```

```r
### sanity check end
```

## Generate bash commands to run QSLIMFinder


```r
SLIMFinder_dir = "./SLIMFinder_Mann/"
if(!dir.exists(SLIMFinder_dir)) dir.create(SLIMFinder_dir)

forSLIMFinder_file_list = writeInteractionSubsetFASTA_list(interactionFASTA_list = forSLIMFinder_Ready,
                                                           dir = SLIMFinder_dir)

QSLIMFinderCommand(file_list = forSLIMFinder_file_list)
```

```
## [1] "bsub -n 1 -q research-rh7 -M 100 -R \"rusage[mem=100]\" python /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/../software/slimsuite/tools/qslimfinder.py blast+path=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/../software/ncbi_blast_2.6.0/bin/ iupath=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/../software/iupred/iupred dismask=T consmask=F cloudfix=F probcut=0.1 resdir=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./SLIMFinder_Mann/output/interactors_of.A0FGR8.P0DOE9./ resfile=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./SLIMFinder_Mann/output/interactors_of.A0FGR8.P0DOE9./main_result seqin=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./SLIMFinder_Mann/input/fasta/interactors_of.A0FGR8.P0DOE9.fas query=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./SLIMFinder_Mann/input/query/interactors_of.A0FGR8.P0DOE9.fas "
```

```r
all_commands = mQSLIMFinderCommand(file_list = forSLIMFinder_file_list, 
                                   slimpath = "../software/cluster/slimsuite/tools/qslimfinder.py", 
                                   blast = "../software/cluster/ncbi_blast_2.6.0/bin/", 
                                   iupred = "../software/cluster/iupred/iupred",
                                   options = "dismask=T consmask=F cloudfix=T probcut=0.3 iuchdir=T",
                                   LSF_cluster_mode = T,
                                   LSF_project_path = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/",
                                   LSF_cluster_par = "bsub -n 1 -q research-rh7 -M 100 -R \"rusage[mem=100]\"",
                                   write_log = F)
```

## Run QSLIMFinder


```r
#if(tryCatch({system("bjobs", intern = T)}, error = function(e) e)$message != "error in running command") onLSF = T else onLSF = F
runQSLIMFinder(commands = all_commands, file_list = forSLIMFinder_file_list, onLSF = T)
```

```
## character(0)
```

## Read results


```r
if(!dir.exists("./SLIMFinder_Mann/result/")) dir.create("./SLIMFinder_Mann/result/")
QSLIMFinder_main_result = readQSLIMFinderMain(outputfile = forSLIMFinder_file_list$outputfile)
QSLIMFinder_occurence = readQSLIMFinderOccurence(outputdir = forSLIMFinder_file_list$outputdir)
fwrite(QSLIMFinder_main_result, "./SLIMFinder_Mann/result/main_result.txt", sep = "\t")
fwrite(QSLIMFinder_occurence, "./SLIMFinder_Mann/result/occurence.txt", sep = "\t")
writePatternList(QSLIMFinder_main_result, filename = "./SLIMFinder_Mann/result/motifs.txt")
```

## Compare motifs to ELM


```r
runCompariMotif3(slimpath = "../software/cluster/slimsuite/tools/",
                 run = T, with = "db",
                 out_file = "./SLIMFinder_Mann/result/comparimotif.tdt")
```

```
## [1] "python /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/../software/cluster/slimsuite/tools/comparimotif_V3.py motifs=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./SLIMFinder/result/motifs.txt searchdb=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project//hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./data_files/2017-10-29elms_index.tsv unmatched=T motific=T resfile=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./SLIMFinder_Mann/result/comparimotif.tdt"
```

```r
runCompariMotif3(slimpath = "../software/cluster/slimsuite/tools/",
                 run = T, with = "self",
                 out_file = "./SLIMFinder_Mann/result/comparimotif_with_self.tdt")
```

```
## [1] "python /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/../software/cluster/slimsuite/tools/comparimotif_V3.py motifs=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./SLIMFinder/result/motifs.txt  unmatched=T motific=T resfile=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./SLIMFinder_Mann/result/comparimotif_with_self.tdt"
```


```r
tar("./SLIMFinder_Mann/input.gz","./SLIMFinder_Mann/input/", compression='gzip')
#tar("./SLIMFinder_Vidal/log_dir.gz","./SLIMFinder_Vidal/log_dir/", compression='gzip')
#tar("./SLIMFinder_Vidal/output.gz","./SLIMFinder_Vidal/output", compression='gzip')
unlink("./SLIMFinder_Mann/input/")
#unlink("./SLIMFinder_Vidal/log_dir/")
#unlink("./SLIMFinder_Vidal/output/")
filename = paste0("./processed_data_files/QSLIMFinder_instances_h2v_Mann_clust",gsub("-","",Sys.Date()),".RData")
save(list = ls(), file=filename)
Sys.Date()
```

```
## [1] "2017-10-29"
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
##  [1] ontologyIndex_2.4    scales_0.5.0         RColorBrewer_1.1-2  
##  [4] GGally_1.3.2         ggplot2_2.2.1        rtracklayer_1.36.3  
##  [7] GenomicRanges_1.28.3 GenomeInfoDb_1.12.2  MItools_0.1.29      
## [10] Biostrings_2.44.1    XVector_0.16.0       data.table_1.10.4-3 
## [13] PSICQUIC_1.14.0      plyr_1.8.4           httr_1.3.1          
## [16] biomaRt_2.32.1       IRanges_2.10.2       S4Vectors_0.14.3    
## [19] BiocGenerics_0.22.0  rmarkdown_1.6       
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.13               lattice_0.20-35           
##  [3] Rsamtools_1.28.0           rprojroot_1.2             
##  [5] digest_0.6.12              R6_2.2.2                  
##  [7] backports_1.1.1            RSQLite_2.0               
##  [9] evaluate_0.10.1            zlibbioc_1.22.0           
## [11] rlang_0.1.2                lazyeval_0.2.0            
## [13] curl_3.0                   blob_1.1.0                
## [15] R.utils_2.5.0              R.oo_1.21.0               
## [17] Matrix_1.2-11              qvalue_2.8.0              
## [19] gsubfn_0.6-6               proto_1.0.0               
## [21] splines_3.4.1              BiocParallel_1.10.1       
## [23] downloader_0.4             stringr_1.2.0             
## [25] RCurl_1.95-4.8             bit_1.1-12                
## [27] munsell_0.4.3              DelayedArray_0.2.7        
## [29] compiler_3.4.1             htmltools_0.3.6           
## [31] SummarizedExperiment_1.6.3 tibble_1.3.4              
## [33] GenomeInfoDbData_0.99.0    matrixStats_0.52.2        
## [35] XML_3.98-1.9               reshape_0.8.7             
## [37] GenomicAlignments_1.12.1   bitops_1.0-6              
## [39] R.methodsS3_1.7.1          grid_3.4.1                
## [41] jsonlite_1.5               gtable_0.2.0              
## [43] DBI_0.7                    magrittr_1.5              
## [45] stringi_1.1.5              reshape2_1.4.2            
## [47] tools_3.4.1                bit64_0.9-7               
## [49] Biobase_2.36.2             yaml_2.1.14               
## [51] AnnotationDbi_1.38.1       colorspace_1.3-2          
## [53] memoise_1.1.0              knitr_1.17
```
