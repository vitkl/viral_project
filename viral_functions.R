## viral project functions
##' function to download from Uniprot and write fasta sequence files (renaming the UniProt-given name of the sequence to UniprotAC only) given generic (P04637) or isoform (P04637-2) UniprotAC. If the file already contains the sequences for given UniprotAC you get a message, if some sequences are missing all sequences will be reloaded.
download_fasta = function(uniprot_ac, file_name){
    if(if(file.exists(file_name)) mean(uniprot_ac %in% fasta.index(file_name)$desc) != 1 else T){
        
        all_fasta = AAStringSet()
        sequences_per_id = numeric(length = length(uniprot_ac))
        names(sequences_per_id) = uniprot_ac
        for(protein in uniprot_ac){
            new_fasta = readAAStringSet(paste0("http://www.uniprot.org/uniprot/", protein,".fasta"))
            # if the search has yielded fasta file with 1 sequence then rename to UniprotAC, otherwise save original name and deal with the problem manually
            if(length(new_fasta) == 1) names(new_fasta) = protein
            all_fasta = append(all_fasta, new_fasta)
            sequences_per_id[protein] = length(new_fasta)
            
            if(which(uniprot_ac == protein) %in% seq(1,length(uniprot_ac),10)){
                # write data at each 10th sequence in case R crashes
                writeXStringSet(all_fasta, file = file_name, format="fasta")
                # limit the number of requests per second to 10
                Sys.sleep(1)
            }
        }
        # write final result
        writeXStringSet(all_fasta, file = file_name, format="fasta")
    }
    if(if(file.exists(file_name)) mean(uniprot_ac %in% fasta.index(file_name)$desc) == 1 else F){
        message("all sequences for given uniprot_ac are alredy downloaded")
    }
}

##' function to download from Uniprot and write fasta sequence files (renaming the UniProt-given name of the sequence to proteinID-FeatureID: P04591-PRO_0000261216) given post-processed chain ID, such as PRO_0000261216. If the file already contains the sequences for given post-processed chain ID you get a message, if some sequences are missing all sequences will be reloaded.
download_fasta_postproc = function(postproc_id, file_name){
    if(if(file.exists(file_name)) mean(postproc_id %in% gsub("^[[:alnum:]]+-","",fasta.index(file_name)$desc)) != 1 else T){
        postproc_fasta = AAStringSet()
        for(feature in postproc_id){
            # Read gff from UniProt
            new_feature = import(paste0("http://www.uniprot.org/uniprot/?query=",feature,"&format=gff"), format="gff")
            # subset features with ID
            new_feature = new_feature[!is.na(new_feature$ID)]
            # select the feature we have searched for
            new_feature = new_feature[new_feature$ID == feature]
            # read FASTA for the protein
            protein = as.character(seqnames(new_feature))
            proteinFASTA = readAAStringSet(paste0("http://www.uniprot.org/uniprot/", protein,".fasta"))
            # select the first and the only sequence in a set and subset it by feature range
            postproc_fasta_new = AAStringSet(proteinFASTA[[1]][ranges(new_feature)])
            names(postproc_fasta_new) = paste0(protein, "-", feature)
            postproc_fasta = append(postproc_fasta, postproc_fasta_new)
            
            if(which(postproc_id == feature) %in% seq(1,length(postproc_id),10)){
                # write data at each 10th sequence in case R crashes
                writeXStringSet(postproc_fasta, file = file_name, format="fasta")
                # limit the number of requests per second to 10
                Sys.sleep(1)
            }
        }
        # write final result
        writeXStringSet(postproc_fasta, file = file_name, format="fasta")
    }
    if(if(file.exists(file_name)) mean(postproc_id %in% gsub("^[[:alnum:]]+-","",fasta.index(file_name)$desc)) == 1 else F){
        message("all sequences for given post-processed chain ID are alredy downloaded")
    }
}