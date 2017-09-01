# knit RMD on cluster
library(rmarkdown)
rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM.Rmd", 
                  output_format = "html_document", 
                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust.html")

# bsub -n 16 -q research-rh7 -M 4000 -R "rusage[mem=4000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM.R
# bsub -n 16 -q research-rh7 -M 12288 -R "rusage[mem=12288]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM.R