# Motif_search_strategies
Vitalii Kleshchevnikov  
29/11/2017  



# Introduction

## Strategies

1. Optimal search settings (3*2):  
- slimlen= 5,8,10 : Maximum length of SLiMs to return (no. non-wildcard positions)  
- maxwild= 2,3 : Maximum number of consecutive wildcard positions to allow  

2. Comparing interaction datasets (4):  
- Full IntAct  
- BioPlex  
- Vidal  
- randomised BioPlex (keeping the degree and the number of interactions intact)  

3. Domain-centered vs protein-centered strategy (2)  

4. QSLIMFinder vs SLIMFinder (2)  

Total: 3 \* 2 \* 4 \* 2 \* 2 = 96 conditions  

This is not feasible, so I will first compare interaction datasets, domain- vs protein- centered strategy and QSLIMFinder vs SLIMFinder (16 conditions).  
I will then test multiple motif settings on a selected dataset (different run of the pipeline, 6 conditions).  

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
Read 0.0% of 790620 rows
Read 6.3% of 790620 rows
Read 15.2% of 790620 rows
Read 16.4% of 790620 rows
Read 22.8% of 790620 rows
Read 31.6% of 790620 rows
Read 36.7% of 790620 rows
Read 44.3% of 790620 rows
Read 51.9% of 790620 rows
Read 59.4% of 790620 rows
Read 68.3% of 790620 rows
Read 77.2% of 790620 rows
Read 86.0% of 790620 rows
Read 96.1% of 790620 rows
Read 790620 rows and 42 (of 42) columns from 2.995 GB file in 00:00:29
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

```
## 
Read 5.4% of 184566 rows
Read 184566 rows and 11 (of 11) columns from 0.029 GB file in 00:00:04
```

```r
# need to simplify isoform and post-processed chain IDs for human proteins
all_viral_interaction_BioPlex = copy(all_viral_interaction)
all_viral_interaction_BioPlex$data$IDs_interactor_A = gsub("-[[:digit:]]+$","",all_viral_interaction_BioPlex$data$IDs_interactor_A)
all_viral_interaction_BioPlex$data$pair_id = paste0(all_viral_interaction_BioPlex$data$IDs_interactor_A, "|", all_viral_interaction_BioPlex$data$IDs_interactor_B)
all_viral_interaction_BioPlex$data = unique(all_viral_interaction_BioPlex$data)

# human-human
Full_IntAct = fullInteractome(taxid = 9606, database = "IntActFTP", format = "tab27",
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
BioPlex = loadBioplex(dir = "./data_files/",
                      url = "http://bioplex.hms.harvard.edu/data/BioPlex_2.3_interactionList.tsv")
```

```
## ... loading local copy ...
```

```r
Vidal = subsetMITABbyPMIDs(MITABdata = Full_IntAct,
                           PMIDs = c("25416956", "unassigned1304"))

randomised_BioPlex = copy(BioPlex)
randomised_BioPlex$data$Publication_Identifiers = NA
randomised_BioPlex$data$Confidence_values = NA
randomised_BioPlex$data = unique(randomised_BioPlex$data)
set.seed(1)
randomised_BioPlex$data$IDs_interactor_B = randomised_BioPlex$data$IDs_interactor_B[sample.int(length(randomised_BioPlex$data$IDs_interactor_B))]
randomised_BioPlex$data$pair_id = paste0(randomised_BioPlex$data$IDs_interactor_A, "|", randomised_BioPlex$data$IDs_interactor_B)
```

## Set standard options


```r
seed_list = c("O00264","O00264-PRO_0000038428","O00469","O14545","O14979","O15031","O15269","O15270","O43149","O43164","O43294","O60518","O60825","O75477","O95639","P00403","P01892","P02795","P05023","P05546-PRO_0000037940","P05546-PRO_0000037946","P06576","P06703","P08195","P08962","P09543","P11166","P13667","P16615","P16615-PRO_0000045596","P16989","P21912","P21964","P27824","P35613","P35658","P38435","P49257","P49755","P50395","P50395-PRO_0000037966","P51570","P51648-PRO_0000045596","P51648","P51659","P52630","P52756","P54136","P61201","P61619","P63151","P82979","Q00587","Q02539","Q03188-PRO_0000037551","Q03518","Q08999","Q13105","Q13586","Q13724","Q13724-PRO_0000045596","Q14677","Q14789","Q15369","Q49A26","Q53H12","Q5C9Z4","Q5H9R7","Q5HYI8","Q5SXM2","Q5T3F8-PRO_0000037548","Q5VTE0","Q86YA3","Q8IZY2","Q8N448","Q8N7W2","Q8NBS9-PRO_0000045602","Q8NI60","Q8TED0","Q8WVX9","Q92538","Q92805","Q969X5","Q96ST2","Q96T23","Q96TA2","Q99442","Q9BQ70","Q9BX40","Q9H6S0-PRO_0000037576","Q9H9L4-PRO_0000037566","Q9NP80-PRO_0000037548","Q9NSY0","Q9NXH8","Q9NZ01","Q9NZM1","Q9UJZ1","Q9UPN7","Q9UPU5","Q9Y3T9-PRO_0000037713","Q9Y4F1","Q9Y4W2","Q9Y5M8","Q9Y6N5")
# set standard options
myPPInetwork2SLIMFinder = function(dataset_name = "SLIMFinder", slimlen = 5, maxwild = 2,
                                   interaction_main_set = Full_IntAct, interaction_query_set = all_viral_interaction,
                                   center_domains = F, analysis_type = "qslimfinder"){
    PPInetwork2SLIMFinder(dataset_name = dataset_name,
                          interaction_main_set = interaction_main_set,
                          interaction_query_set = all_viral_interaction,
                          analysis_type = analysis_type,
                          options = paste0("dismask=T consmask=F cloudfix=T probcut=0.3 minwild=0 maxwild=",maxwild," slimlen=",slimlen," alphahelix=F maxseq=1500 savespace=0 iuchdir=T extras=2"),
                          path2domain_enrich = "./processed_data_files/what_we_find_VS_ELM_clust20171019.RData",
                          domain_enrich_object = "res_count",
                          center_domains = center_domains,
                          fasta_path = "./data_files/all_human_viral_proteins.fasta",
                          main_set_only = F,
                          domain_pvalue_cutoff = 0.000005, # 0.000005 for test run, 1 for full run
                          SLIMFinder_dir = paste0("./",dataset_name,"/"),
                          LSF_project_path = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/",
                          software_path = "../software/cluster/",
                          length_set1_min = 4, # selected based on benchmarking: minimal set size that returns motifs with Sig < 0.05
                          length_set2_min = 1,
                          write_log = T,
                          N_seq = 200,
                          seed_list = seed_list, # NULL for full run
                          memory_start = 350,
                          memory_step = 1000)
}

# Download ELM data
elm_filename = paste0("/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/", "./data_files/", Sys.Date(), 
                      "elms_index.tsv")
if (!file.exists(elm_filename)) 
    download.file("http://elm.eu.org/elms/elms_index.tsv", elm_filename)
```

## Define parameters


```r
analysis_type = c("qslimfinder", "slimfinder")
datasets = c("Full_IntAct", "BioPlex", "Vidal", "randomised_BioPlex", "all_viral_interaction")
domain_centered = c(F, T) # c(F, T)
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
dataset_names[!grepl("BioPlex",datasets), query_dataset := "all_viral_interaction"]
```

## Perform motif search


```r
# set up parallel processing
# create cluster
cores = 8
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
## 
## [[4]]
##  [1] "rtracklayer"   "GenomicRanges" "GenomeInfoDb"  "MItools"      
##  [5] "Biostrings"    "XVector"       "data.table"    "PSICQUIC"     
##  [9] "plyr"          "httr"          "biomaRt"       "IRanges"      
## [13] "S4Vectors"     "stats4"        "BiocGenerics"  "parallel"     
## [17] "stats"         "graphics"      "grDevices"     "utils"        
## [21] "datasets"      "methods"       "base"         
## 
## [[5]]
##  [1] "rtracklayer"   "GenomicRanges" "GenomeInfoDb"  "MItools"      
##  [5] "Biostrings"    "XVector"       "data.table"    "PSICQUIC"     
##  [9] "plyr"          "httr"          "biomaRt"       "IRanges"      
## [13] "S4Vectors"     "stats4"        "BiocGenerics"  "parallel"     
## [17] "stats"         "graphics"      "grDevices"     "utils"        
## [21] "datasets"      "methods"       "base"         
## 
## [[6]]
##  [1] "rtracklayer"   "GenomicRanges" "GenomeInfoDb"  "MItools"      
##  [5] "Biostrings"    "XVector"       "data.table"    "PSICQUIC"     
##  [9] "plyr"          "httr"          "biomaRt"       "IRanges"      
## [13] "S4Vectors"     "stats4"        "BiocGenerics"  "parallel"     
## [17] "stats"         "graphics"      "grDevices"     "utils"        
## [21] "datasets"      "methods"       "base"         
## 
## [[7]]
##  [1] "rtracklayer"   "GenomicRanges" "GenomeInfoDb"  "MItools"      
##  [5] "Biostrings"    "XVector"       "data.table"    "PSICQUIC"     
##  [9] "plyr"          "httr"          "biomaRt"       "IRanges"      
## [13] "S4Vectors"     "stats4"        "BiocGenerics"  "parallel"     
## [17] "stats"         "graphics"      "grDevices"     "utils"        
## [21] "datasets"      "methods"       "base"         
## 
## [[8]]
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
    myPPInetwork2SLIMFinder(dataset_name = dataset_names$dataset_names[i], slimlen = dataset_names$slimlen[i], maxwild = dataset_names$maxwild[i],
                            interaction_main_set = eval(parse(text = dataset_names$datasets[i])),
                            interaction_query_set = eval(parse(text = dataset_names$query_dataset[i])),
                            center_domains = dataset_names$domain_centered[i], analysis_type = dataset_names$analysis_type[i])
}, simplify = TRUE, USE.NAMES = TRUE) # parSapplyLB is parSapply that accounts for different time each dataset might take

stopCluster(cl)

# paths to RData files containing all objects created by all runs of this pipeline
print(R_data)
```

```
##  [1] "./processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.Full_IntAct.FALSE_clust201712.RData"          
##  [2] "./processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.Full_IntAct.TRUE_clust201712.RData"           
##  [3] "./processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.BioPlex.FALSE_clust201712.RData"              
##  [4] "./processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.BioPlex.TRUE_clust201712.RData"               
##  [5] "./processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.Vidal.FALSE_clust201712.RData"                
##  [6] "./processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.Vidal.TRUE_clust201712.RData"                 
##  [7] "./processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.randomised_BioPlex.FALSE_clust201712.RData"   
##  [8] "./processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.randomised_BioPlex.TRUE_clust201712.RData"    
##  [9] "./processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.all_viral_interaction.FALSE_clust201712.RData"
## [10] "./processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.all_viral_interaction.TRUE_clust201712.RData" 
## [11] "./processed_data_files/QSLIMFinder_instances_h2v_slimfinder.Full_IntAct.FALSE_clust201712.RData"           
## [12] "./processed_data_files/QSLIMFinder_instances_h2v_slimfinder.Full_IntAct.TRUE_clust201712.RData"            
## [13] "./processed_data_files/QSLIMFinder_instances_h2v_slimfinder.BioPlex.FALSE_clust201712.RData"               
## [14] "./processed_data_files/QSLIMFinder_instances_h2v_slimfinder.BioPlex.TRUE_clust201712.RData"                
## [15] "./processed_data_files/QSLIMFinder_instances_h2v_slimfinder.Vidal.FALSE_clust201712.RData"                 
## [16] "./processed_data_files/QSLIMFinder_instances_h2v_slimfinder.Vidal.TRUE_clust201712.RData"                  
## [17] "./processed_data_files/QSLIMFinder_instances_h2v_slimfinder.randomised_BioPlex.FALSE_clust201712.RData"    
## [18] "./processed_data_files/QSLIMFinder_instances_h2v_slimfinder.randomised_BioPlex.TRUE_clust201712.RData"     
## [19] "./processed_data_files/QSLIMFinder_instances_h2v_slimfinder.all_viral_interaction.FALSE_clust201712.RData" 
## [20] "./processed_data_files/QSLIMFinder_instances_h2v_slimfinder.all_viral_interaction.TRUE_clust201712.RData"
```
