# Identify viral-human interactions and download protein sequences
Vitalii Kleshchevnikov  
22/06/2017  




## Seaching for viral-human PPI  

I download all protein-protein interactions from IntAct (IMEx databases) involving viral proteins using PSICQUIC service in MITAB2.7 format (https://psicquic.github.io/MiqlReference27.html).


```r
query_PSICQUIC(query = "species:10239",
                format = "tab27",
                database = "imex",
                file = "./data_files/human_viral_interactions.txt")
```

```
## Warning in if (N_interactions > 0) {: the condition has length > 1 and only
## the first element will be used
```

```
## Warning in N_SPECIES_ID_interactome[indices] = N_interactions: number of
## items to replace is not a multiple of replacement length
```

```
##             query                                      file format
##  1: species:10239 ./data_files/human_viral_interactions.txt  tab27
##  2: species:10239 ./data_files/human_viral_interactions.txt  tab27
##  3: species:10239 ./data_files/human_viral_interactions.txt  tab27
##  4: species:10239 ./data_files/human_viral_interactions.txt  tab27
##  5: species:10239 ./data_files/human_viral_interactions.txt  tab27
##  6: species:10239 ./data_files/human_viral_interactions.txt  tab27
##  7: species:10239 ./data_files/human_viral_interactions.txt  tab27
##  8: species:10239 ./data_files/human_viral_interactions.txt  tab27
##  9: species:10239 ./data_files/human_viral_interactions.txt  tab27
## 10: species:10239 ./data_files/human_viral_interactions.txt  tab27
## 11: species:10239 ./data_files/human_viral_interactions.txt  tab27
##     all.databases n.interactions.in.database database.not.active
##  1:        IntAct                      10329                    
##  2:          MINT                       8770                    
##  3:       bhf-ucl                          0                    
##  4:         MPIDB                          8                    
##  5:      MatrixDB                     <html>                    
##  6:         HPIDb                       3452                    
##  7:      I2D-IMEx                         25                    
##  8: InnateDB-IMEx                          0                    
##  9:        MolCon                          2                    
## 10:       UniProt                        696                    
## 11:        MBInfo                          0
```

```r
all_viral_interaction = fread("./data_files/human_viral_interactions.txt", stringsAsFactors = F)
```

I clean the data in the table to make it more useble. Then, I filter and keep only human-viral interactions.


```r
colnames(all_viral_interaction) = unlist(strsplit(readLines("ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.txt",  n = 1), "\t"))

# changing column names to data.table-compatible format
{
colnames(all_viral_interaction) = gsub(" ","_",colnames(all_viral_interaction))
colnames(all_viral_interaction) = gsub("\\(|\\)","",colnames(all_viral_interaction))
colnames(all_viral_interaction) = gsub("#","",colnames(all_viral_interaction))
}
# cleaning Taxid "taxid:9606(human)|taxid:9606(Homo sapiens)" to 9606
{
all_viral_interaction[, Taxid_interactor_A := gsub("taxid:|\\(.*$","",Taxid_interactor_A)]
all_viral_interaction[, Taxid_interactor_B := gsub("taxid:|\\(.*$","",Taxid_interactor_B)]
all_viral_interaction[, Host_organisms := gsub("taxid:|\\(.*$","",Host_organisms)]
# saving identifier types and cleaning interactor ids
all_viral_interaction[, interactor_IDs_databases_A := gsub(":.*$","",IDs_interactor_A)]
all_viral_interaction[, interactor_IDs_databases_B := gsub(":.*$","",IDs_interactor_B)]
all_viral_interaction[, IDs_interactor_A := gsub("^.*:","",IDs_interactor_A)]
all_viral_interaction[, IDs_interactor_B := gsub("^.*:","",IDs_interactor_B)]
# isoform "-1" is a canonical sequence, IntAct uses isoform "-1" when it's clear that the isoform is "-1" and a canonical identifier if it's not clear which isoform was used in the experiment. Removing isoform sign "-1":
all_viral_interaction[, IDs_interactor_A := gsub("-1$", "", IDs_interactor_A)]
all_viral_interaction[, IDs_interactor_B := gsub("-1$", "", IDs_interactor_B)]
# removing interactions if at least one interactor has non-uniprot id
all_viral_interaction = all_viral_interaction[interactor_IDs_databases_A == "uniprotkb",][interactor_IDs_databases_B == "uniprotkb",]
# cleaning other information
all_viral_interaction[, bait_prey_status_A := gsub("^.*\\(|\\)","",Experimental_roles_interactor_A)]
all_viral_interaction[, bait_prey_status_B := gsub("^.*\\(|\\)","",Experimental_roles_interactor_B)]
all_viral_interaction[, Publication_Identifiers := gsub("^.*pubmed:|\\|.*$","",Publication_Identifiers)]
all_viral_interaction[, Confidence_values := gsub("^intact-miscore:","",Confidence_values)]
all_viral_interaction[, Confidence_values := gsub("-","NA",Confidence_values)]
all_viral_interaction[, Confidence_values := as.numeric(Confidence_values)]
#all_viral_interaction[, Interaction_identifiers := unlist(gsubfn::strapplyc(Interaction_identifiers,"EBI-[[:digit:]]+",simplify = T)), by =Interaction_identifiers]
# generating unique identifier for interacting pairs
all_viral_interaction[, pair_id := apply(data.table(IDs_interactor_A,IDs_interactor_B,stringsAsFactors = F), 1,
                                               function(a) { z = sort(a)
                                               paste0(z[1],"_",z[2]) })]
all_viral_interaction[, pair_species_id := apply(data.table(Taxid_interactor_A,Taxid_interactor_B,stringsAsFactors = F), 1,
                                               function(a) { z = sort(a)
                                               paste0(z[1],"_",z[2]) })]
}
```

```
## Warning in eval(jsub, SDenv, parent.frame()): NAs introduced by coercion
```

```r
# filter only human-viral interactions
all_viral_interaction = all_viral_interaction[Taxid_interactor_A == "9606" | Taxid_interactor_B == "9606",]
all_viral_interaction = unique(all_viral_interaction)

# select proteins
viral_proteins = unique(c(all_viral_interaction[Taxid_interactor_A != "9606", IDs_interactor_A], all_viral_interaction[Taxid_interactor_B != "9606", IDs_interactor_B]))
human_proteins = unique(c(all_viral_interaction[Taxid_interactor_A == "9606", IDs_interactor_A], all_viral_interaction[Taxid_interactor_B == "9606", IDs_interactor_B]))
```

There are 787 viral_proteins (including isoforms and postprocessed chains).  
There are 4423 human_proteins (including isoforms and postprocessed chains).  

## Finding protein sequences

Interacting viral and human proteins belong to three categories:  
- proteins identified by canonical UniProt identifier (P04637): no isoform exist or impossible to distinguish isoforms from the published result  
- proteins identified to an isoform (P04637-1)  
- proteins that are cleaved into the functional fragments from a precursor protein (P04591-PRO_0000261216)  

These categories require different approaches of retrieving sequences.   
1. Retrieving sequences of proteins with canonical UniProt identifier is possible with R package for UniProt webservices or UniProt REST API. These proteins also may not require InterProScan to identify domains.  
2. Isoform sequences are accessible using UniProt REST API: http://www.uniprot.org/uniprot/P04637-2.fasta.
3. Post-processed chains are not straightforward to map as they are defined by the interval of the canonical UniProt sequence (http://www.uniprot.org/uniprot/P04591.fasta). First, we need to retrieve post-processed chain position from Uniprot in gff format: http://www.uniprot.org/uniprot/?query=PRO_0000261216&format=gff. The search return all sequence features from a given protein, we select only post-processed chain we are interested in and then use position specified to subset the sequence.

### Retrieving canonical and isoform sequences


```r
# Using UniProt.ws package to retrieve sequences
# uniprot = UniProt.ws(taxId = 9606)
# select(uniprot, keys = "P04637", columns = c("UNIPROTKB","SEQUENCE"), keytype = "UNIPROTKB")

# Filtering sequence names by group
canonical_human = human_proteins[-c(grep("-[[:digit:]]+", human_proteins, value = F),grep("-PRO_[[:digit:]]+$", human_proteins, value = F))]
canonical_viral = viral_proteins[-c(grep("-[[:digit:]]+", viral_proteins, value = F),grep("-PRO_[[:digit:]]+$", viral_proteins, value = F))]
isoform_human = grep("-[[:digit:]]+", human_proteins, value = T)
isoform_viral = grep("-[[:digit:]]+", viral_proteins, value = T)
postproc_human = grep("-PRO_[[:digit:]]+$", human_proteins, value = T)
postproc_viral = grep("-PRO_[[:digit:]]+$", viral_proteins, value = T)

canonical_and_isoform = c(canonical_human, canonical_viral, isoform_human, isoform_viral)

# dowloading FASTA for canonical_and_isoform sequences if file with this data doesn't exist or doesn't contain all sequences
# download_fasta(uniprot_ac = canonical_and_isoform, file_name = "./data_files/canonical_and_isoform.fasta")
```

Some protein ids have been UniParc-ed. I will comment the code necessary to download the sequences instead of dealing with the problem as the problem is being generated by de-sync between UniProt and IntAct releases and is, therefore, temporary.  

### Retrieving post-processed chain sequences


```r
# Filtering sequence names by group
postproc_human = grep("-PRO_[[:digit:]]+$", human_proteins, value = T)
postproc_viral = grep("-PRO_[[:digit:]]+$", viral_proteins, value = T)
postproc = c(gsub("^[[:alnum:]]+-","",postproc_human), gsub("^[[:alnum:]]+-","",postproc_viral))

# dowloading FASTA for post-processed sequences if file with this data doesn't exist or doesn't contain all sequences
download_fasta_postproc(postproc_id = postproc, file_name = "./data_files/postproc.fasta")
```

```
## all sequences for given post-processed chain ID are alredy downloaded
```

### Combining sequences


```r
all_human_viral_proteins.fasta = append(readAAStringSet("./data_files/canonical_and_isoform.fasta"), readAAStringSet("./data_files/postproc.fasta"))
writeXStringSet(all_human_viral_proteins.fasta,
                file = "./data_files/all_human_viral_proteins.fasta",
                format="fasta")
```

5218 FASTA sequences total.
