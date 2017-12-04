# testing r packages
# clusterProfiler
# create documentation and namespace from ##' comments in the .R
library(devtools)
library(roxygen2)
setwd("~/Desktop/vitalii/other/clusterProfiler")
document()
# install
setwd("..")
install("clusterProfiler")

# test
library(devtools)
install_github(repo = "vitkl/clusterProfiler")
library(clusterProfiler)
library(data.table)
# create proteins annotated with domains
set.seed(1)
x = data.table(protein = paste0("P",c(83:100,sample(1:100, 162, replace = T),81:100)), domain = c(rep("A",20),rep("B",20),rep("C",30),rep("D",40),rep("E",10),rep("F",20),rep("G",20),rep("H",40)))
# test categ_dist without cluster mode
x_simili_noclust = categ_dist(mapping_table = x, parallel = F)
# test categ_dist in cluster mode
x_simili_clust = categ_dist(mapping_table = x, parallel = T)
all.equal(x_simili_noclust, x_simili_clust)
# create proteins annotated with domains: large
set.seed(1)
x = data.table(protein = paste0("P",c(83:100,sample(1:9962, 9962, replace = T),81:100)), domain = c(rep("D1",18), paste0("D", sample(1:200, 9962, replace = T)),rep("D200",20)))
# test categ_dist without cluster mode
time1 = proc.time()
x_simili_noclust = categ_dist(mapping_table = x, parallel = F)
(time1end = proc.time() - time1)
# test categ_dist in cluster mode
time2 = proc.time()
x_simili_clust = categ_dist(mapping_table = x, parallel = T)
(time2end = proc.time() - time2)
# test categ_dist in cluster mode
time2 = proc.time()
x_simili_clust = categ_dist(mapping_table = x, parallel = T, cores_to_use = 3)
(time2end = proc.time() - time2)
# test categ_dist in cluster mode
time3 = proc.time()
x_simili_clust = categ_dist(mapping_table = x, parallel = T, cores_to_use = 2)
(time2end = proc.time() - time3)

all.equal(x_simili_noclust, x_simili_clust)

# testing on cluster
bsub -n 8 -q research-rh7 -M 8000 -R "rusage[mem=8000]" 
bsub -n 4 -q research-rh7 -M 8000 -R "rusage[mem=8000]" 