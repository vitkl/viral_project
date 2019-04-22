# Pipeline for de-novo discovery of SLIMs using viral-human protein interaction network 

## cluster command line commands to run code in this RScript (pipeline_for_EBIcluster.R)
# bsub -n 1 -q research-rh7 -M 16000 -R "rusage[mem=16000]" -o /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.log Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R
# bsub -n 4 -q research-rh7 -M 64000 -R "rusage[mem=64000]" -o /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.log Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R
# bsub -n 8 -q research-rh7 -M 16000 -R "rusage[mem=16000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R
## for domain enrichment:
# bsub -n 32 -q research-rh7 -M 64000 -R "rusage[mem=64000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R

## Pipeline
######## uncomment relevant parts of the code
library(rmarkdown)
proj_path = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/"
proj_path = "/nfs/research1/petsalaki/users/vitalii/vitalii/viral_project/"

### 1 Download interaction data, retrieve protein sequences, run InterProScan to predict domains
# bsub -n 1 -q research-rh7 -M 8000 -R "rusage[mem=8000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R
#rmarkdown::render(input = paste0(proj_path, "interactions_and_sequences.Rmd"), 
#                  output_format = "html_document", 
#                  output_file=paste0(proj_path, "interactions_and_sequences_clust.html"))

### 2 Run InterProScan to identify domains
# bsub -n 32 -q research-rh7 -M 64000 -R "rusage[mem=64000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R
#rmarkdown::render(input = paste0(proj_path, "interactions_and_sequences.Rmd"), 
#                  output_format = "html_document", 
#                  output_file=paste0(proj_path, "interactions_and_sequences_clust.html"))

### 3 Remove redundant domain predictions
# bsub -n 1 -q research-rh7 -M 8000 -R "rusage[mem=8000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R
#rmarkdown::render(input = paste0(proj_path, "remove_redundant_domains.Rmd"), 
#                  output_format = "html_document", 
#                  output_file=paste0(proj_path, "remove_redundant_domains_clust.html"))

### 4 Combine domain and protein interaction data
# bsub -n 1 -q research-rh7 -M 8000 -R "rusage[mem=8000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R
#rmarkdown::render(input = paste0(proj_path, "map_domains_to_human_viral_network.Rmd"), 
#                  output_format = "html_document", 
#                  output_file=paste0(proj_path, "map_domains_to_human_viral_network_clust.html"))

### 4 Estimate which domains are likely to mediate interactions of viral with human proteins
# bsub -n 32 -q research-rh7 -M 64000 -R "rusage[mem=64000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R
#rmarkdown::render(input = paste0(proj_path, "predict_domain_viral_clust_count.Rmd"), 
#                  output_format = "html_document", 
#                  output_file=paste0(proj_path, "predict_domain_viral_clust_count.html"))

### 4.1 (optional) Estimate which domains are likely to mediate interactions between human proteins
#bsub -n 1 -q research-rh7 -M 32000 -R "rusage[mem=32000] select[mem>32000]" Rscript /nfs/research1/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R
#/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R
#rmarkdown::render(input = paste0(proj_path, "predict_domain_human_clust_count.Rmd"), 
#                  output_format = "html_document", 
#                  output_file=paste0(proj_path, "predict_domain_human_clust_count.html"))

### 6 Motif search using multiple datasets: IntAct Vidal only_viral
#rmarkdown::render(input = paste0(proj_path, "Motif_search_strategies_IntAct_Vidal_viral.Rmd"), 
#                  output_format = "html_document", 
#                  output_file=paste0(proj_path, "Motif_search_strategies_IntAct_Vidal_viral.html"))

### 7 Motif prediction benchmarking using known motif instances (from ELM database) + other output for the paper
# bsub -n 1 -q research-rh7 -M 16000 -R "rusage[mem=16000] select[mem>16000]" Rscript /nfs/research1/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R
# /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/pipeline_for_EBIcluster.R
rmarkdown::render(input = paste0(proj_path, "compr_benchmarking_venn.Rmd"), 
                  output_format = "html_document", 
                  output_file=paste0(proj_path, "compr_benchmarking_venn_clust.html"))

### 8. Extra
### 8.1 Analyse degree distribution in human and human-viral protein interaction network
#rmarkdown::render(input = paste0(proj_path, "Degree_distribution_in_the_network.Rmd"), 
#                  output_format = "html_document", 
#                  output_file=paste0(proj_path, "Degree_distribution_in_the_network_clust.html"))

####### IGNORE everything below ############


# ???????????????
# In progress: check which datasets were already processed by motif search pipeline 
# cat /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/RData_from_Motif_search_4_known_motifs

### Motif search using multiple strategies
#rmarkdown::render(input = paste0(proj_path, "Motif_search_strategies.Rmd"), 
#                  output_format = "html_document", 
#                  output_file=paste0(proj_path, "Motif_search_strategies_clust.html"))

# In progress: check which datasets were already processed by motif search pipeline 
# cat /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/RData_from_Motif_search_strategies
# end ????????????????????????????



#bsub -n 16 -q research-rh7 -M 16000 -R "rusage[mem=16000]" -Is $SHELL
#library(rmarkdown);
#rmarkdown::render(input = paste0(proj_path, "what_we_find_VS_ELM_copy.Rmd"), 
#                  output_format = "html_document", 
#                  output_file=paste0(proj_path, "what_we_find_VS_ELM_copy_clust.html"))

### TEST for reproducibility

### Motif search using multiple strategies: IntAct Vidal only_viral
#rmarkdown::render(input = paste0(proj_path, "Motif_search_strategies_IntAct_Vidal_viral.Rmd"), 
#                  output_format = "html_document", 
#                  output_file=paste0(proj_path, "Motif_search_strategies_IntAct_Vidal_viral.html"))

#rmarkdown::render(input = paste0(proj_path, "Motif_search_strategies_IntAct_Vidal_viral2.Rmd"), 
#                  output_format = "html_document", 
#                  output_file=paste0(proj_path, "Motif_search_strategies_IntAct_Vidal_viral_clust.html"))

### re-run the same pipeline as an Rmd file
#rmarkdown::render(input = paste0(proj_path, "QSLIMFinder_instances_h2v.Rmd"), 
#                  output_format = "html_document", 
#                  output_file=paste0(proj_path, "QSLIMFinder_instances_h2v_FullIntAct2.html"))

#rmarkdown::render(input = paste0(proj_path, "QSLIMFinder_instances_h2v_Vidal.Rmd"), 
#                  output_format = "html_document", 
#                  output_file=paste0(proj_path, "QSLIMFinder_instances_h2v_Vidal2.html"))
