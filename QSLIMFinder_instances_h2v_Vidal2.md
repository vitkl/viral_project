# QSLIMFinder: domain instances, human to search (Vidal data only), viral to filter
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
Read 0.0% of 794782 rows
Read 6.3% of 794782 rows
Read 15.1% of 794782 rows
Read 16.4% of 794782 rows
Read 22.6% of 794782 rows
Read 30.2% of 794782 rows
Read 36.5% of 794782 rows
Read 42.8% of 794782 rows
Read 49.1% of 794782 rows
Read 55.4% of 794782 rows
Read 59.1% of 794782 rows
Read 66.7% of 794782 rows
Read 73.0% of 794782 rows
Read 80.5% of 794782 rows
Read 88.1% of 794782 rows
Read 96.9% of 794782 rows
Read 794782 rows and 42 (of 42) columns from 3.017 GB file in 00:00:31
```

```r
all_human_interaction = subsetMITABbyPMIDs(MITABdata = all_human_interaction,
                                           PMIDs = c("25416956", "unassigned1304"))
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
Read 0.0% of 794782 rows
Read 7.5% of 794782 rows
Read 15.1% of 794782 rows
Read 21.4% of 794782 rows
Read 28.9% of 794782 rows
Read 35.2% of 794782 rows
Read 41.5% of 794782 rows
Read 46.6% of 794782 rows
Read 51.6% of 794782 rows
Read 56.6% of 794782 rows
Read 64.2% of 794782 rows
Read 70.5% of 794782 rows
Read 76.8% of 794782 rows
Read 83.0% of 794782 rows
Read 91.8% of 794782 rows
Read 98.1% of 794782 rows
Read 794782 rows and 42 (of 42) columns from 3.017 GB file in 00:00:26
## 
Read 37.5% of 186853 rows
Read 186853 rows and 11 (of 11) columns from 0.030 GB file in 00:00:04
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

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/QSLIMFinder_instances_h2v_Vidal2_files/figure-html/load_results-1.png)<!-- -->

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
                                           seed_id_vect = proteins_w_signif_domains, # c("O00459", "O00571") proteins_w_signif_domains
                                           fasta = all.fasta,
                                           single_interact_from_set2 = T, set1_only = F)
```

```
## Warning in listSingleInteractFromSet2(subset1 = subset,
## single_interact_from_set2, : set 2 is empty, name:P27449
```

```r
forSLIMFinder_Ready = filterInteractionSubsetFASTA_list(forSLIMFinder,  length_set1_min = 2, length_set2_min = 1)
# how many viral-human pairs in which viral proteins interact with human protein (with relevant domain) but are not associated with the domain
sum(!domainProteinPairMatch(forSLIMFinder_Ready, domain_res, remove = F))
```

```
## [1] 3
```

```r
# how many viral-human pairs in which viral proteins interact with human protein (with relevant domain) and are associated with the domain
sum(domainProteinPairMatch(forSLIMFinder_Ready, domain_res, remove = F))
```

```
## [1] 10269
```

```r
# total
forSLIMFinder_Ready$length
```

```
## [1] 10272
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
res_ = tstrsplit(gsub("interactors_of\\.|\\.","", names(forSLIMFinder_Ready$fasta_subset_list)), ":")
query_in_FASTA = sapply(1:length(forSLIMFinder_Ready$fasta_subset_list), function(i, forSLIMFinder_, res_){
    res_[[2]][i] %in% names(forSLIMFinder_$fasta_subset_list[[i]]) # query from name in FASTA
    #res_[[2]][i] %in% forSLIMFinder_$interaction_subset[[i]]$ids_set2 # query from name is the query from interaction set 
}, forSLIMFinder_Ready, res_)
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
SLIMFinder_dir = "./Vidal2/"
if(!dir.exists(SLIMFinder_dir)) dir.create(SLIMFinder_dir)

forSLIMFinder_file_list = writeInteractionSubsetFASTA_list(interactionFASTA_list = forSLIMFinder_Ready,
                                                           dir = SLIMFinder_dir)

QSLIMFinderCommand(file_list = forSLIMFinder_file_list)
```

```
## [1] "bsub -n 1 -q research-rh7 -M 100 -R \"rusage[mem=100]\" python /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/../software/slimsuite/tools/qslimfinder.py blast+path=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/../software/ncbi_blast_2.6.0/bin/ iupath=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/../software/iupred/iupred dismask=T consmask=F cloudfix=F probcut=0.1 resdir=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Vidal2/output/interactors_of.A0FGR8.P0DOE9./ resfile=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Vidal2/output/interactors_of.A0FGR8.P0DOE9./main_result seqin=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Vidal2/input/fasta/interactors_of.A0FGR8.P0DOE9.fas query=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Vidal2/input/query/interactors_of.A0FGR8.P0DOE9.fas "
```

```r
all_commands = mQSLIMFinderCommand(file_list = forSLIMFinder_file_list, 
                                   slimpath = "../software/cluster/slimsuite/tools/", 
                                   blast = "../software/cluster/ncbi_blast_2.6.0/bin/", 
                                   iupred = "../software/cluster/iupred/iupred",
                                   options = "dismask=T consmask=F cloudfix=T probcut=0.3 iuchdir=T",
                                   LSF_cluster_mode = T,
                                   LSF_project_path = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/",
                                   LSF_cluster_par = "bsub -n 1 -q research-rh7 -M 100 -R \"rusage[mem=100]\"",
                                   write_log = F, log_dir = "./Vidal2/log_dir/")
```

## Run QSLIMFinder


```r
#if(tryCatch({system("bjobs", intern = T)}, error = function(e) e)$message != "error in running command") onLSF = T else onLSF = F
runQSLIMFinder(commands_list = all_commands, file_list = forSLIMFinder_file_list, onLSF = T)
```

```
## character(0)
```

## Read results


```r
if(!dir.exists("./Vidal2/result/")) dir.create("./Vidal2/result/")
QSLIMFinder_main_result = readQSLIMFinderMain(outputfile = forSLIMFinder_file_list$outputfile)
QSLIMFinder_occurence = readQSLIMFinderOccurence(outputdir = forSLIMFinder_file_list$outputdir)
fwrite(QSLIMFinder_main_result, "./Vidal2/result/main_result.txt", sep = "\t")
fwrite(QSLIMFinder_occurence, "./Vidal2/result/occurence.txt", sep = "\t")
writePatternList(QSLIMFinder_main_result, filename = "./Vidal2/result/motifs.txt")
```

## Compare motifs to ELM


```r
runCompariMotif3(input_file = "./Vidal2/result/motifs.txt",
                 slimpath = "../software/cluster/slimsuite/tools/",
                 run = T, with = "db",
                 out_file = "./Vidal2/result/comparimotif.tdt")
```

```
## [1] "python /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/../software/cluster/slimsuite/tools/comparimotif_V3.py motifs=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Vidal2/result/motifs.txt searchdb=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project//hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./data_files/2018-02-19elms_index.tsv unmatched=T motific=T resfile=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Vidal2/result/comparimotif.tdt"
```

```r
runCompariMotif3(input_file = "./Vidal2/result/motifs.txt",
                 slimpath = "../software/cluster/slimsuite/tools/",
                 run = T, with = "self",
                 out_file = "./Vidal2/result/comparimotif_with_self.tdt")
```

```
## [1] "python /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/../software/cluster/slimsuite/tools/comparimotif_V3.py motifs=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Vidal2/result/motifs.txt  unmatched=T motific=T resfile=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Vidal2/result/comparimotif_with_self.tdt"
```

```r
#R.utils::gzip("./Vidal2/result/comparimotif_with_self.compare.tdt", "./Vidal2/result/comparimotif_with_self.compare.tdt.gz")
#R.utils::gzip("./Vidal2/result/comparimotif_with_self.compare.xgmml", "./Vidal2/result/comparimotif_with_self.compare.xgmml.gz")
```


```r
tar("./Vidal2/input.gz","./Vidal2/input/", compression='gzip')
#tar("./Vidal2/log_dir.gz","./Vidal2/log_dir/", compression='gzip')
#tar("./Vidal2/output.gz","./Vidal2/output", compression='gzip')
unlink("./Vidal2/input/")
#unlink("./Vidal2/log_dir/")
#unlink("./Vidal2/output/")
filename = paste0("./processed_data_files/QSLIMFinder_instances_h2v_Vidal2_clust",format(Sys.Date(), "%Y%m"),".RData")
save(list = ls(), file=filename)
Sys.Date()
```

```
## [1] "2018-02-19"
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
##  [7] GenomicRanges_1.28.3 GenomeInfoDb_1.12.2  MItools_0.1.36      
## [10] Biostrings_2.44.1    XVector_0.16.0       data.table_1.10.4-3 
## [13] PSICQUIC_1.14.0      plyr_1.8.4           httr_1.3.1          
## [16] biomaRt_2.32.1       IRanges_2.10.2       S4Vectors_0.14.3    
## [19] BiocGenerics_0.22.0  rmarkdown_1.6       
## 
## loaded via a namespace (and not attached):
##  [1] Biobase_2.36.2             bit64_0.9-7               
##  [3] jsonlite_1.5               splines_3.4.1             
##  [5] gsubfn_0.6-6               R.utils_2.6.0             
##  [7] gtools_3.5.0               blob_1.1.0                
##  [9] GenomeInfoDbData_0.99.0    Rsamtools_1.28.0          
## [11] yaml_2.1.16                RSQLite_2.0               
## [13] backports_1.1.1            lattice_0.20-35           
## [15] downloader_0.4             digest_0.6.12             
## [17] qvalue_2.8.0               colorspace_1.3-2          
## [19] htmltools_0.3.6            Matrix_1.2-11             
## [21] R.oo_1.21.0                XML_3.98-1.9              
## [23] zlibbioc_1.22.0            gdata_2.18.0              
## [25] BiocParallel_1.10.1        tibble_1.3.4              
## [27] DT_0.4                     SummarizedExperiment_1.6.3
## [29] ROCR_1.0-7                 lazyeval_0.2.0            
## [31] proto_1.0.0                magrittr_1.5              
## [33] memoise_1.1.0              evaluate_0.10.1           
## [35] R.methodsS3_1.7.1          gplots_3.0.1              
## [37] tools_3.4.1                matrixStats_0.52.2        
## [39] stringr_1.2.0              munsell_0.4.3             
## [41] DelayedArray_0.2.7         AnnotationDbi_1.38.1      
## [43] compiler_3.4.1             caTools_1.17.1            
## [45] rlang_0.1.2                grid_3.4.1                
## [47] RCurl_1.95-4.10            htmlwidgets_1.0           
## [49] bitops_1.0-6               gtable_0.2.0              
## [51] curl_3.0                   DBI_0.7                   
## [53] reshape_0.8.7              reshape2_1.4.2            
## [55] R6_2.2.2                   GenomicAlignments_1.12.1  
## [57] knitr_1.17                 bit_1.1-12                
## [59] rprojroot_1.2              KernSmooth_2.23-15        
## [61] stringi_1.1.6              Rcpp_0.12.13
```
