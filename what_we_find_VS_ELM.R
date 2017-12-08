# knit RMD on cluster
library(rmarkdown)
#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/interactions_and_sequences.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/interactions_and_sequences_clust.html")

#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/remove_redundant_domains.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/remove_redundant_domains_clust.html")

#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/map_domains_to_human_viral_network.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/map_domains_to_human_viral_network_clust.html")

#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_count_justFisher.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_count_justFisher_clust.html")

#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/Degree_distribution_in_the_network.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/Degree_distribution_in_the_network_clust.html")

rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/Motif_search_strategies.Rmd", 
                  output_format = "html_document", 
                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/Motif_search_strategies_clust.html")


# bsub -n 1 -q research-rh7 -M 8000 -R "rusage[mem=8000]" -o /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/SLIMFinder/QSLIMFinder_instances_h2v.log Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM.R
# bsub -n 8 -q research-rh7 -M 32000 -R "rusage[mem=32000]" -o /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM.log Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM.R
# bsub -n 8 -q research-rh7 -M 16000 -R "rusage[mem=16000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM.R
# bsub -n 32 -q research-rh7 -M 64000 -R "rusage[mem=64000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM.R

#bsub -n 16 -q research-rh7 -M 16000 -R "rusage[mem=16000]" -Is $SHELL
#library(rmarkdown);
#rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_copy.Rmd", 
#                  output_format = "html_document", 
#                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_copy_clust.html")
