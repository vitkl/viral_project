# Motif_search_strategies
Vitalii Kleshchevnikov  
29/11/2017  



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

Total: 1 \* 1 \* 2 \* 1 \* 1 = 2 conditions  

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
Read 16.4% of 794782 rows
Read 22.6% of 794782 rows
Read 32.7% of 794782 rows
Read 36.5% of 794782 rows
Read 46.6% of 794782 rows
Read 55.4% of 794782 rows
Read 59.1% of 794782 rows
Read 69.2% of 794782 rows
Read 79.3% of 794782 rows
Read 89.3% of 794782 rows
Read 99.4% of 794782 rows
Read 794782 rows and 42 (of 42) columns from 3.017 GB file in 00:00:24
```

```r
# human-viral
all_viral_interaction2 = interSpeciesInteractome(taxid1 = 9606, taxid2 = 10239, database = "IntActFTP", format = "tab27",
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
# human-human
Full_IntAct2 = fullInteractome(taxid = 9606, database = "IntActFTP", format = "tab27",
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
Vidal2 = subsetMITABbyPMIDs(MITABdata = Full_IntAct2,
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
                          length_set1_min = 4, # selected based on benchmarking: minimal set size that returns motifs with Sig < 0.05 = 4+1
                          length_set2_min = 1,
                          write_log = T,
                          N_seq = 200,
                          seed_list = NULL, # NULL for full run, seed_list for test run
                          memory_start = 350,
                          memory_step = 2000)
}

# Download ELM data
elm_filename = paste0("/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/", "./data_files/", Sys.Date(), 
                      "elms_index.tsv")
if (!file.exists(elm_filename)) 
    download.file("http://elm.eu.org/elms/elms_index.tsv", elm_filename)
```

## Define parameters


```r
analysis_type = c("qslimfinder")
datasets = c("Full_IntAct2", "Vidal2", "all_viral_interaction2")
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
dataset_names[!grepl("BioPlex",datasets), query_dataset := "all_viral_interaction2"]

# remove already finished names
R_data_names = readLines("./RData_from_Motif_search_strategies")
R_data_names = gsub("\\./processed_data_files/QSLIMFinder_instances_h2v_","", R_data_names)
R_data_names = gsub("_clust201712\\.RData\\.zip","", R_data_names)
R_data_names = gsub("_clust201801\\.RData\\.zip","", R_data_names)
R_data_names = gsub("_clust201802\\.RData\\.zip","", R_data_names)
dataset_names = dataset_names[!dataset_names %in% R_data_names]
```

## Perform motif search


```r
# set up parallel processing
# create cluster
cores = 3
cl <- makeCluster(cores)
# get library support needed to run the code
clusterEvalQ(cl, {library(MItools); library(rtracklayer)})
```

```
## [[1]]
##  [1] "rtracklayer"   "GenomicRanges" "GenomeInfoDb"  "MItools"      
##  [5] "Biostrings"    "XVector"       "data.table"    "PSICQUIC"     
##  [9] "plyr"          "httr"          "biomaRt"       "IRanges"      
## [13] "S4Vectors"     "stats4"        "BiocGenerics"  "parallel"     
## [17] "stats"         "graphics"      "grDevices"     "utils"        
## [21] "datasets"      "methods"       "base"         
## 
## [[2]]
##  [1] "rtracklayer"   "GenomicRanges" "GenomeInfoDb"  "MItools"      
##  [5] "Biostrings"    "XVector"       "data.table"    "PSICQUIC"     
##  [9] "plyr"          "httr"          "biomaRt"       "IRanges"      
## [13] "S4Vectors"     "stats4"        "BiocGenerics"  "parallel"     
## [17] "stats"         "graphics"      "grDevices"     "utils"        
## [21] "datasets"      "methods"       "base"         
## 
## [[3]]
##  [1] "rtracklayer"   "GenomicRanges" "GenomeInfoDb"  "MItools"      
##  [5] "Biostrings"    "XVector"       "data.table"    "PSICQUIC"     
##  [9] "plyr"          "httr"          "biomaRt"       "IRanges"      
## [13] "S4Vectors"     "stats4"        "BiocGenerics"  "parallel"     
## [17] "stats"         "graphics"      "grDevices"     "utils"        
## [21] "datasets"      "methods"       "base"
```

```r
# put objects in place that might be needed for the code
clusterExport(cl, ls(), envir=environment())

R_data = parSapplyLB(cl = cl, X = 1:nrow(dataset_names), FUN = function(i){
    R_data_name = myPPInetwork2SLIMFinder(dataset_name = dataset_names$dataset_names[i], 
                            slimlen = dataset_names$slimlen[i], maxwild = dataset_names$maxwild[i],
                            interaction_main_set = eval(parse(text = dataset_names$datasets[i])),
                            interaction_query_set = eval(parse(text = dataset_names$query_dataset[i])),
                            center_domains = dataset_names$domain_centered[i], 
                            analysis_type = dataset_names$analysis_type[i])
    
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
```

```
## [1] "./processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.all_viral_interaction2.FALSE_clust201802.RData.zip"
```
