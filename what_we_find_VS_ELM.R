# knit RMD on cluster
library(rmarkdown)
rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/map_domains_to_network.Rmd", 
                  output_format = "html_document", 
                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/map_domains_to_network_clust.html")
rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM.Rmd", 
                  output_format = "html_document", 
                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust.html")

# bsub -n 16 -q research-rh7 -M 16000 -R "rusage[mem=16000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM.R
# bsub -n 32 -q research-rh7 -M 32000 -R "rusage[mem=32000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM.R