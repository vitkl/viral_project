# QSLIMFinder: domain instances, human to search, viral to filter
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
IntAct = loadIntActFTP(dir = "./data_files/IntActRelease_2017Nov13/")
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
Read 44.0% of 794782 rows
Read 50.3% of 794782 rows
Read 56.6% of 794782 rows
Read 59.1% of 794782 rows
Read 66.7% of 794782 rows
Read 73.0% of 794782 rows
Read 79.3% of 794782 rows
Read 86.8% of 794782 rows
Read 94.4% of 794782 rows
Read 99.4% of 794782 rows
Read 794782 rows and 42 (of 42) columns from 3.017 GB file in 00:00:36
```

```r
# human-human
all_human_interaction = fullInteractome(taxid = 9606, database = "IntActFTP", format = "tab27",
                                        clean = TRUE, protein_only = TRUE,
                                        MITABdata = IntAct, directory = "./data_files/")
```

```
## ... looking for the date of the latest IntAct release ...
```

```
## Warning in download.file(input, tt, method = method, mode = "wb", quiet = !
## showProgress): downloaded length 591 != reported length 0
```

```r
# human-viral
all_viral_interaction = interSpeciesInteractome(taxid1 = 9606, taxid2 = 10239, database = "IntActFTP", format = "tab27",
                                                clean = TRUE, protein_only = TRUE,
                                                MITABdata = IntAct, directory = "./data_files/")
```

```
## ... looking for the date of the latest IntAct release ...
```

```
## Warning in download.file(input, tt, method = method, mode = "wb", quiet = !
## showProgress): downloaded length 591 != reported length 0
```

```r
dataset_name = "Full_IntAct2" # refer to MItools::mBenchmarkMotifs
SLIMFinder_dir = paste0("./",dataset_name,"/")
LSF_project_path = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/"
analysis_type = "qslimfinder" # or "slimfinder"
options = "dismask=T consmask=F cloudfix=T probcut=0.3 iuchdir=T"
software_path = "../software/cluster/" # "../software/"
write_log = T
```

## load results of which domains are likely to mediate interaction

```r
# load the domain analysis results
load("./processed_data_files/what_we_find_VS_ELM_clust20171019.RData")
domain_res = res_count
rm(list = ls()[!ls() %in% c("domain_res", "all_human_interaction", "all_viral_interaction", "dataset_name", "SLIMFinder_dir", "LSF_project_path", "analysis_type", "options", "software_path", "write_log")])

# choose pvalue cutoff:
plot(domain_res)
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/QSLIMFinder_instances_h2v_FullIntAct2_files/figure-html/load_results-1.png)<!-- -->

```r
proteins_w_signif_domains = unique(domain_res$data_with_pval[p.value <= 1, IDs_interactor_human])
#proteins_w_signif_domains = extractInteractors(all_viral_interaction, taxid = 9606, inverse_filter = T)
```

### prepare sets of protein sequences and instructions for SLIMFinder


```r
# load FASTA
all.fasta = readAAStringSet(filepath = "./data_files/all_human_viral_proteins.fasta", format = "fasta")

# remove interactions if one of the interactors does not have a FASTA sequence
all_human_interaction = removeInteractionNoFASTA(all_human_interaction, all.fasta)
all_viral_interaction = removeInteractionNoFASTA(all_viral_interaction, all.fasta)

forSLIMFinder = listInteractionSubsetFASTA(interaction_set1 = all_human_interaction, 
                                           interaction_set2 = all_viral_interaction, 
                                           seed_id_vect = proteins_w_signif_domains,
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
## [1] 4
```

```r
# how many viral-human pairs in which viral proteins interact with human protein (with relevant domain) and are associated with the domain
sum(domainProteinPairMatch(forSLIMFinder_Ready, domain_res, remove = F))
```

```
## [1] 11830
```

```r
# total
forSLIMFinder_Ready$length
```

```
## [1] 11834
```

```r
#filter
domain_filt = domain_res
domain_filt$data_with_pval = domain_filt$data_with_pval[p.value <= 1,]
forSLIMFinder_Ready = domainProteinPairMatch(forSLIMFinder_Ready, domain_filt, remove = T)
```

## Generate bash commands to run QSLIMFinder


```r
if(!dir.exists(SLIMFinder_dir)) dir.create(SLIMFinder_dir)

forSLIMFinder_file_list = writeInteractionSubsetFASTA_list(interactionFASTA_list = forSLIMFinder_Ready,
                                                           dir = SLIMFinder_dir, analysis_type = analysis_type)

QSLIMFinderCommand(file_list = forSLIMFinder_file_list, 
                                   slimpath = paste0(software_path, "slimsuite/tools/"), 
                                   blast = paste0(software_path, "ncbi_blast_2.6.0/bin/"), 
                                   iupred = paste0(software_path, "iupred/iupred"),
                                   options = options,
                                   LSF_project_path = LSF_project_path,
                                   LSF_cluster_par = "bsub -n 1 -q research-rh7 -M 100 -R \"rusage[mem=100]\"",
                                   analysis_type = analysis_type)
```

```
## [1] "bsub -n 1 -q research-rh7 -M 100 -R \"rusage[mem=100]\" python /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/../software/cluster/slimsuite/tools/qslimfinder.py blast+path=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/../software/cluster/ncbi_blast_2.6.0/bin/ iupath=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/../software/cluster/iupred/iupred dismask=T consmask=F cloudfix=T probcut=0.3 iuchdir=T resdir=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Full_IntAct2/output/interactors_of.A0FGR8.P0DOE9./ resfile=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Full_IntAct2/output/interactors_of.A0FGR8.P0DOE9./main_result seqin=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Full_IntAct2/input/fasta/interactors_of.A0FGR8.P0DOE9.fas query=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Full_IntAct2/input/query/interactors_of.A0FGR8.P0DOE9.fas "
```

```r
all_commands = mQSLIMFinderCommand(file_list = forSLIMFinder_file_list, 
                                   slimpath = paste0(software_path, "slimsuite/tools/"), 
                                   blast = paste0(software_path, "ncbi_blast_2.6.0/bin/"), 
                                   iupred = paste0(software_path, "iupred/iupred"),
                                   options = options,
                                   LSF_cluster_mode = T,
                                   LSF_project_path = LSF_project_path,
                                   LSF_cluster_par = "bsub -n 1 -q research-rh7 -M 100 -R \"rusage[mem=100]\"",
                                   analysis_type = analysis_type,
                                   write_log = write_log,
                                   log_dir = paste0(SLIMFinder_dir, "log_dir/"))
#all_commands = groupQSLIMFinderCommand(commands = all_commands,
#                                      InteractionSubsetFASTA_list = forSLIMFinder_Ready,
#                                       sh_dir = "/sh_dir/",
#                                       LSF_project_path = LSF_project_path,
#                                       dataset_name = dataset_name, N_seq = 200, write_log = write_log)
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
resultdir = paste0(SLIMFinder_dir, "result/")
if(!dir.exists(resultdir)) dir.create(resultdir)
QSLIMFinder_main_result = readQSLIMFinderMain(outputfile = forSLIMFinder_file_list$outputfile)
QSLIMFinder_occurence = readQSLIMFinderOccurence(outputdir = forSLIMFinder_file_list$outputdir)
fwrite(QSLIMFinder_main_result, paste0(resultdir, "main_result.txt"), sep = "\t")

fwrite(QSLIMFinder_occurence, paste0(resultdir, "occurence.txt"), sep = "\t")
writePatternList(QSLIMFinder_main_result, filename = paste0(resultdir, "motifs.txt"))
```

## Compare motifs to ELM


```r
runCompariMotif3(input_file = paste0(resultdir, "motifs.txt"),
                 slimpath = paste0(software_path, "slimsuite/tools/"),
                 run = T, with = "db",
                 out_file = paste0(resultdir, "comparimotif.tdt"))
```

```
## [1] "python /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/../software/cluster/slimsuite/tools/comparimotif_V3.py motifs=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Full_IntAct2/result/motifs.txt searchdb=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project//hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./data_files/2018-02-18elms_index.tsv unmatched=T motific=T resfile=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Full_IntAct2/result/comparimotif.tdt"
```

```r
runCompariMotif3(input_file = paste0(resultdir, "motifs.txt"),
                 slimpath = paste0(software_path, "slimsuite/tools/"),
                 run = T, with = "self",
                 out_file = paste0(resultdir, "comparimotif_with_self.tdt"))
```

```
## [1] "python /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/../software/cluster/slimsuite/tools/comparimotif_V3.py motifs=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Full_IntAct2/result/motifs.txt  unmatched=T motific=T resfile=/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/./Full_IntAct2/result/comparimotif_with_self.tdt"
```

```r
#R.utils::gzip(paste0(resultdir, "comparimotif_with_self.compare.tdt"), paste0(resultdir, "comparimotif_with_self.compare.tdt.gz"))
#R.utils::gzip(paste0(resultdir, "comparimotif_with_self.compare.xgmml"), paste0(resultdir, "comparimotif_with_self.compare.xgmml.gz"))
```


```r
tar(paste0(SLIMFinder_dir, "input.gz"),paste0(SLIMFinder_dir, "input/"), compression='gzip')
#tar(paste0(SLIMFinder_dir, "log_dir.gz"),paste0(SLIMFinder_dir, "log_dir/"), compression='gzip')
#tar(paste0(SLIMFinder_dir, "output.gz"),paste0(SLIMFinder_dir, "output"), compression='gzip')
unlink(paste0(SLIMFinder_dir, "input/"), recursive = T)
#unlink(paste0(SLIMFinder_dir, "log_dir/"))
#unlink(paste0(SLIMFinder_dir, "output/"))
filename = paste0("./processed_data_files/QSLIMFinder_instances_h2v_",dataset_name,"_clust",format(Sys.Date(), "%Y%m"),".RData")
save(list = ls(), file=filename)
Sys.Date()
```

```
## [1] "2018-02-18"
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
##  [1] rtracklayer_1.36.3   GenomicRanges_1.28.3 GenomeInfoDb_1.12.2 
##  [4] MItools_0.1.36       Biostrings_2.44.1    XVector_0.16.0      
##  [7] data.table_1.10.4-3  PSICQUIC_1.14.0      plyr_1.8.4          
## [10] httr_1.3.1           biomaRt_2.32.1       IRanges_2.10.2      
## [13] S4Vectors_0.14.3     BiocGenerics_0.22.0  rmarkdown_1.6       
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
## [17] RColorBrewer_1.1-2         qvalue_2.8.0              
## [19] colorspace_1.3-2           htmltools_0.3.6           
## [21] Matrix_1.2-11              R.oo_1.21.0               
## [23] XML_3.98-1.9               zlibbioc_1.22.0           
## [25] scales_0.5.0               gdata_2.18.0              
## [27] BiocParallel_1.10.1        tibble_1.3.4              
## [29] ggplot2_2.2.1              DT_0.4                    
## [31] SummarizedExperiment_1.6.3 ROCR_1.0-7                
## [33] lazyeval_0.2.0             proto_1.0.0               
## [35] magrittr_1.5               memoise_1.1.0             
## [37] evaluate_0.10.1            GGally_1.3.2              
## [39] R.methodsS3_1.7.1          gplots_3.0.1              
## [41] tools_3.4.1                matrixStats_0.52.2        
## [43] stringr_1.2.0              munsell_0.4.3             
## [45] DelayedArray_0.2.7         AnnotationDbi_1.38.1      
## [47] compiler_3.4.1             caTools_1.17.1            
## [49] rlang_0.1.2                grid_3.4.1                
## [51] RCurl_1.95-4.10            htmlwidgets_1.0           
## [53] bitops_1.0-6               gtable_0.2.0              
## [55] curl_3.0                   DBI_0.7                   
## [57] reshape_0.8.7              reshape2_1.4.2            
## [59] R6_2.2.2                   GenomicAlignments_1.12.1  
## [61] knitr_1.17                 bit_1.1-12                
## [63] rprojroot_1.2              ontologyIndex_2.4         
## [65] KernSmooth_2.23-15         stringi_1.1.6             
## [67] Rcpp_0.12.13
```
