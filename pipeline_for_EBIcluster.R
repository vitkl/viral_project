# Pipeline for de-novo discovery of SLIMs using viral-human protein interaction network 

## cluster command line commands to run code in this RScript (pipeline_for_EBIcluster.R)
# bsub -n 1 -q research-rh7 -M 8000 -R "rusage[mem=8000]" -o /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.log Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R
# bsub -n 3 -q research-rh7 -M 25000 -R "rusage[mem=25000]" -o /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.log Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R
# bsub -n 8 -q research-rh7 -M 16000 -R "rusage[mem=16000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R
## for domain enrichment:
# bsub -n 32 -q research-rh7 -M 64000 -R "rusage[mem=64000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R

## Pipeline
######## uncomment relevant parts of the code
library(rmarkdown)

### Download interaction data, retrieve protein sequences, run InterProScan to predict domains
#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/interactions_and_sequences.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/interactions_and_sequences_clust.html")

### Remove redundant domain predictions
#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/remove_redundant_domains.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/remove_redundant_domains_clust.html")

### Combine domain and protein interaction data
#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/map_domains_to_human_viral_network.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/map_domains_to_human_viral_network_clust.html")

### Estimate which domains are likely to mediate interaction
#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_count_justFisher.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_count_justFisher_clust.html")

### Analyse degree distribution in human and human-viral protein interaction network
#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/Degree_distribution_in_the_network.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/Degree_distribution_in_the_network_clust.html")

### Motif search using multiple strategies
#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/Motif_search_strategies.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/Motif_search_strategies_clust.html")

# In progress: check which datasets were already processed by motif search pipeline 
# cat /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/RData_from_Motif_search_strategies

### Motif prediction benchmarking using known motif instances (from ELM database)
# rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_strateg.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/compr_benchmarking_strateg.html")



#bsub -n 16 -q research-rh7 -M 16000 -R "rusage[mem=16000]" -Is $SHELL
#library(rmarkdown);
#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_copy.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_copy_clust.html")

### TEST for reproducibility

### Motif search using multiple strategies: IntAct Vidal only_viral
rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/Motif_search_strategies_IntAct_Vidal_viral.Rmd", 
                  output_format = "html_document", 
                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/Motif_search_strategies_IntAct_Vidal_viral.html")

#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/Motif_search_strategies_IntAct_Vidal_viral2.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/Motif_search_strategies_IntAct_Vidal_viral_set_len2.html")

### re-run the same pipeline as an Rmd file
#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/QSLIMFinder_instances_h2v.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/QSLIMFinder_instances_h2v_FullIntAct2.html")

#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/QSLIMFinder_instances_h2v_Vidal.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/QSLIMFinder_instances_h2v_Vidal2.html")
