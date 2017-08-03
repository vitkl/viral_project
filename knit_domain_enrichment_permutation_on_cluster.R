# knit RMD on cluster
if(parallel::detectCores() > 4){
    library(rmarkdown)
    rmarkdown::render(input = "/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/domain_enrichment_permutation.Rmd", 
                      output_format = "html_document", 
                      output_file="/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/domain_enrichment_permutation_clust.html")
}

# bsub -n 16 -q research-rh7 -M 16000 -R "rusage[mem=16000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/knit_domain_enrichment_permutation_on_cluster.R
# bsub -n 16 -q research-rh7 -M 12288 -R "rusage[mem=12288]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/knit_domain_enrichment_permutation_on_cluster.R
# bsub -n 16 -q research-rh7 -M 20000 -R "rusage[mem=20000]" Rscript /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/knit_domain_enrichment_permutation_on_cluster.R

# knit locally
if(parallel::detectCores() <= 4){

#########################################################################################################
directory = "./"
# knit locally foldEnrichment AND remove_zeros FALSE

library(rmarkdown)
RMD = readLines("./domain_enrichment_permutation.Rmd")
if(grep("frequency", RMD[51])) RMD[51] = "frequency = FALSE"
if(grep("remove_zeros", RMD[53])) RMD[53] = "remove_zeros = FALSE"

if(grep("statistic", RMD[64])) RMD[64] = "- **statistic:** fold enrichment of a domain among interacting partners of viral protein  "
if(grep("proteins with missing domains", RMD[68])) RMD[68] = "- **proteins with missing domains:** inlcude for permutations, inlcude in permutation-derived distribution before calculating p-value  "

write(RMD, paste0(directory,"domain_enrichment_permutation_foldEnr_remove_zerosFALSE.Rmd"))

rmarkdown::render(input = paste0(directory,"domain_enrichment_permutation_foldEnr_remove_zerosFALSE.Rmd"), 
                  output_format = "html_document", 
                  output_file= paste0(directory,"domain_enrichment_permutation_foldEnr_remove_zerosFALSE.html"))

#rm(list = ls()[ls() != "directory"])
#########################################################################################################
directory = "./"
# knit locally foldEnrichment AND remove_zeros = TRUE

library(rmarkdown)
RMD = readLines("./domain_enrichment_permutation.Rmd")
if(grep("frequency", RMD[51])) RMD[51] = "frequency = FALSE"
if(grep("remove_zeros", RMD[53])) RMD[53] = "remove_zeros = TRUE"

if(grep("statistic", RMD[64])) RMD[64] = "- **statistic:** fold enrichment of a domain among interacting partners of viral protein  "
if(grep("proteins with missing domains", RMD[68])) RMD[68] = "- **proteins with missing domains:** inlcude for permutations, remove from permutation-derived distribution before calculating p-value  "

write(RMD, paste0(directory,"domain_enrichment_permutation_foldEnr_remove_zerosTRUE.Rmd"))

rmarkdown::render(input = paste0(directory,"domain_enrichment_permutation_foldEnr_remove_zerosTRUE.Rmd"), 
                  output_format = "html_document", 
                  output_file= paste0(directory,"domain_enrichment_permutation_foldEnr_remove_zerosTRUE.html"))

#rm(list = ls()[ls() != "directory"])
#########################################################################################################
directory = "./"
# knit locally frequency AND remove_zeros = FALSE

RMD = readLines("./domain_enrichment_permutation.Rmd")
if(grep("frequency", RMD[51])) RMD[51] = "frequency = TRUE"
if(grep("remove_zeros", RMD[53])) RMD[53] = "remove_zeros = FALSE"

if(grep("statistic", RMD[64])) RMD[64] = "- **statistic:** frequency of a domain among interacting partners of viral protein  "
if(grep("proteins with missing domains", RMD[68])) RMD[68] = "- **proteins with missing domains:** inlcude for permutations, inlcude in permutation-derived distribution before calculating p-value  "

write(RMD, paste0(directory,"domain_enrichment_permutation_frequency_remove_zerosFALSE.Rmd"))

library(rmarkdown)
rmarkdown::render(input = paste0(directory,"domain_enrichment_permutation_frequency_remove_zerosFALSE.Rmd"), 
                  output_format = "html_document", 
                  output_file=paste0(directory,"domain_enrichment_permutation_frequency_remove_zerosFALSE.html"))

#rm(list = ls()[ls() != "directory"])
#########################################################################################################
directory = "./"
# knit locally frequency AND remove_zeros = TRUE

RMD = readLines("./domain_enrichment_permutation.Rmd")
if(grep("frequency", RMD[51])) RMD[51] = "frequency = TRUE"
if(grep("remove_zeros", RMD[53])) RMD[53] = "remove_zeros = TRUE"

if(grep("statistic", RMD[64])) RMD[64] = "- **statistic:** frequency of a domain among interacting partners of viral protein  "
if(grep("proteins with missing domains", RMD[68])) RMD[68] = "- **proteins with missing domains:** inlcude for permutations, remove from permutation-derived distribution before calculating p-value  "

write(RMD, paste0(directory,"domain_enrichment_permutation_frequency_remove_zerosTRUE.Rmd"))

library(rmarkdown)
rmarkdown::render(input = paste0(directory,"domain_enrichment_permutation_frequency_remove_zerosTRUE.Rmd"), 
                  output_format = "html_document", 
                  output_file=paste0(directory,"domain_enrichment_permutation_frequency_remove_zerosTRUE.html"))

#rm(list = ls()[ls() != "directory"])
#########################################################################################################
}
