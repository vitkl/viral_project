source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite(c("PSICQUIC", "qvalue"))
devtools::install_github("vitkl/MItools", dependencies = T)
library(MItools)
# load ppi data
IntAct = loadIntActFTP(dir = "./data_files/IntActRelease_2017Nov13/")

# human-human
Full_IntAct = fullInteractome(taxid = 9606, database = "IntActFTP", format = "tab27",
                              clean = TRUE, protein_only = TRUE,
                              MITABdata = IntAct, directory = "./data_files/")
Full_IntAct_degree = edgelist2degree(Full_IntAct$data)
BioPlex = loadBioplex(dir = "./data_files/",
                      url = "http://bioplex.hms.harvard.edu/data/BioPlex_2.3_interactionList.tsv", uniprot_id = F)
BioPlex_degree = edgelist2degree(BioPlex$data)
# subset both published and unpublished Vidal group data
Vidal = subsetMITABbyPMIDs(MITABdata = Full_IntAct,
                           PMIDs = c("25416956", "unassigned1304"))
Vidal_degree = edgelist2degree(Vidal$data)
# subset Mattias Mann 2015 paper data
Mann = subsetMITABbyPMIDs(MITABdata = full,
                          PMIDs = "26496610")
Mann_degree = edgelist2degree(Mann$data)