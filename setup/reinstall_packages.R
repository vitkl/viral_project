# reinstall packages
# write(rownames(installed.packages()), paste0("installed_packages_",Sys.Date(),".txt"))

packages = read.delim(paste0("/Users/vitalii/Desktop/toggle_switches/installed_packages_","2017-04-30",".txt"), header = F, stringsAsFactors = F)$V1
# install CRAN packages
if(mean(packages %in% rownames(installed.packages())) != 1) {
      install.packages(packages[!packages %in% rownames(installed.packages())])
}
# install Bioconductor core
if(mean(packages %in% rownames(installed.packages())) != 1) {
      source("https://bioconductor.org/biocLite.R")
      biocLite()
}
# install Bioconductor packages
if(mean(packages %in% rownames(installed.packages())) != 1) {
      source("https://bioconductor.org/biocLite.R")
      biocLite(packages[!packages %in% rownames(installed.packages())])
}

library(devtools)
install_github(paste0("genomicsclass/",c("dagdata", "ERBS", "GSE5859", "GSE5859Subset", "maPooling", "ph525x", "tissuesGeneExpression")))

install_github(repo = "vitkl/viral_project", subdir = "queryPSICQUIC")
install_github(repo = "vitkl/clusterProfiler")
install_github(repo = "vitkl/DOSE")
install_github(repo = "vitkl/GOSemSim")
