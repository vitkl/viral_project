# knit RMD on cluster
library(rmarkdown)
rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/interactions_and_sequences.Rmd", 
                  output_format = "html_document", 
                  output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/interactions_and_sequences_clust.html")

# bsub -n 16 -q research-rh7 -M 16000 -R "rusage[mem=16000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/knit_interactions_and_sequences_on_cluster.R
# bsub -n 16 -q research-rh7 -M 12288 -R "rusage[mem=12288]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/knit_interactions_and_sequences_on_cluster.R