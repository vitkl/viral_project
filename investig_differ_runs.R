


# load ppi data
IntAct = loadIntActFTP(dir = "../viral_project/data_files/IntActRelease_2017Nov13/", release = "2017Nov13")
#... loading local copy ...
#Read 790620 rows and 42 (of 42) columns from 2.995 GB file in 00:00:19
# human-human
all_human_interaction = fullInteractome(taxid = 9606, database = "IntActFTP", format = "tab27",
                                        clean = TRUE, protein_only = TRUE,
                                        MITABdata = IntAct, directory = "../viral_project/data_files/",
                                        releaseORdate = "2017Nov13")
# human-viral
all_viral_interaction = interSpeciesInteractome(taxid1 = 9606, taxid2 = 10239, database = "IntActFTP", format = "tab27",
                                                clean = TRUE, protein_only = TRUE,
                                                MITABdata = IntAct, directory = "../viral_project/data_files/",
                                                releaseORdate = "2017Nov13")
mod_new_datasets = lapply(mod_new_datasets, function(mod_new_dataset) {
    (sum(all_viral_interaction$data$IDs_interactor_A %in% mod_new_dataset[1] & all_viral_interaction$data$IDs_interactor_B %in% mod_new_dataset[2]) >= 1)
})
mod_new_datasets = Reduce(c, mod_new_datasets)
mean(mod_new_datasets)

############################ replicating old pipeline INPUT
Full_IntAct2 = R.utils::env(load("/Users/vk7/Desktop/ebi_projects/viral_project/processed_data_files/QSLIMFinder_instances_h2v_Full_IntAct2_clust201802.RData"))
Vidal2 = R.utils::env(load("/Users/vk7/Desktop/ebi_projects/viral_project/processed_data_files/QSLIMFinder_instances_h2v_Vidal2_clust201802.RData"))
Full_IntAct1 = R.utils::env(load("/Users/vk7/Desktop/ebi_projects/viral_project/processed_data_files/QSLIMFinder_instances_h2v_Full_IntAct_clust201710.Rdata"))
Vidal1 = R.utils::env(load("/Users/vk7/Desktop/ebi_projects/viral_project/processed_data_files/QSLIMFinder_instances_h2v_Vidal_clust201710.RData"))
Full_IntAct3 = R.utils::env(load("/Users/vk7/Desktop/ebi_projects/viral_project/processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.Full_IntAct3.FALSE_clust201802.RData"))
Vidal3 = R.utils::env(load("/Users/vk7/Desktop/ebi_projects/viral_project/processed_data_files/QSLIMFinder_instances_h2v_qslimfinder.Vidal3.FALSE_clust201802.RData"))

# sequences in file list checker
seqFileCheck = function(dataset, old = T) {
    # reorder sequence names
    datasets1 = dataset$forSLIMFinder_file_list[,.(interactors_of, QSLIMFinder_query, sequences)]
    datasets1[, sequences := {
        seq_names = strsplit(datasets1[,sequences], "\\.")
        seq_names = lapply(seq_names, sort)
        seq_names = lapply(seq_names, paste0, collapse = ".")
        unlist(seq_names)
    }]
    in_inter_data1 = sapply(1:nrow(datasets1), function(ind, datasets1, dataset){
        if(old){
            sum(dataset$all_human_interaction$data[IDs_interactor_A %in% datasets1$interactors_of[ind],
                                                   mean(unlist(strsplit(datasets1[ind,sequences], "\\.")) %in% IDs_interactor_B)],
                dataset$all_human_interaction$data[IDs_interactor_B %in% datasets1$interactors_of[ind],
                                                   mean(unlist(strsplit(datasets1[ind,sequences], "\\.")) %in% IDs_interactor_A)],
                dataset$all_viral_interaction$data[IDs_interactor_A %in% datasets1$interactors_of[ind],
                                                   mean(unlist(strsplit(datasets1[ind,sequences], "\\.")) %in% IDs_interactor_B)])
        } else {
            sum(dataset$interaction_main_set$data[IDs_interactor_A %in% datasets1$interactors_of[ind],
                                                  mean(unlist(strsplit(datasets1[ind,sequences], "\\.")) %in% IDs_interactor_B)],
                dataset$interaction_main_set$data[IDs_interactor_B %in% datasets1$interactors_of[ind],
                                                  mean(unlist(strsplit(datasets1[ind,sequences], "\\.")) %in% IDs_interactor_A)],
                dataset$interaction_query_set$data[IDs_interactor_A %in% datasets1$interactors_of[ind],
                                                   mean(unlist(strsplit(datasets1[ind,sequences], "\\.")) %in% IDs_interactor_B)])
        }
        
    }, datasets1, dataset)
    list(datasets = datasets1, in_inter_data = in_inter_data, frac_in_inter_data = mean(in_inter_data1))
}

# all datasets are constructed from original interaction data
Vidal11 = seqFileCheck(Vidal1)
Vidal22 = seqFileCheck(Vidal2)
Vidal33 = seqFileCheck(Vidal3, old = F)
Vidal11$frac_in_inter_data
Vidal22$frac_in_inter_data
Vidal33$frac_in_inter_data

Full_IntAct11 = seqFileCheck(Full_IntAct1)
Full_IntAct22 = seqFileCheck(Full_IntAct2)
Full_IntAct33 = seqFileCheck(Full_IntAct3, old = F)
Full_IntAct11$frac_in_inter_data
Full_IntAct22$frac_in_inter_data
Full_IntAct33$frac_in_inter_data

# comparing old and new interactions datasets
# human-human
mean(Vidal1$all_human_interaction$data$pair_id %in% Vidal2$all_human_interaction$data$pair_id)
mean(Full_IntAct1$all_human_interaction$data$pair_id %in% Full_IntAct2$all_human_interaction$data$pair_id)
# huma-viral
mean(Full_IntAct1$all_viral_interaction$data$pair_id %in% Full_IntAct2$all_viral_interaction$data$pair_id)

# QSLIMFinder datasets are the same, except ~30
mean(Full_IntAct11$datasets[,paste0(interactors_of, QSLIMFinder_query)] %in%
         Full_IntAct22$datasets[,paste0(interactors_of, QSLIMFinder_query)])
mean(Full_IntAct22$datasets[,paste0(interactors_of, QSLIMFinder_query)] %in%
         Full_IntAct11$datasets[,paste0(interactors_of, QSLIMFinder_query)])
# but 20% contain one or more different sequence
mean(Full_IntAct11$datasets[,paste0(interactors_of, QSLIMFinder_query, sequences)] %in%
         Full_IntAct22$datasets[,paste0(interactors_of, QSLIMFinder_query, sequences)])
mean(Full_IntAct22$datasets[,paste0(interactors_of, QSLIMFinder_query, sequences)] %in%
         Full_IntAct11$datasets[,paste0(interactors_of, QSLIMFinder_query, sequences)])

# old pipeline produces the same datasets as new pipeline
mean(Full_IntAct22$datasets[,paste0(interactors_of, QSLIMFinder_query, sequences)] %in%
         Full_IntAct33$datasets[,paste0(interactors_of, QSLIMFinder_query, sequences)])
mean(Full_IntAct33$datasets[,paste0(interactors_of, QSLIMFinder_query, sequences)] %in%
         Full_IntAct22$datasets[,paste0(interactors_of, QSLIMFinder_query, sequences)])

# what if original pipeline did no filtering (how many sequences in each set) at all? no - filtering was done (https://github.com/vitkl/viral_project/blob/12b87236f01aaa4056eab6e7faca23d4d3fbb5ab/QSLIMFinder_instances_h2v.Rmd)
# Are the qslimfinder commands the same?
mean(Full_IntAct1$all_commands$run %in% Full_IntAct2$all_commands$run)
# nope, first need to delete dir path that different
run_tmp1 = gsub("/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/\\./SLIMFinder", "", Full_IntAct1$all_commands$run)
run_tmp2 = gsub("/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/\\./Full_IntAct2", "", Full_IntAct2$all_commands$run)
# almost identical commands
mean(run_tmp1 %in% run_tmp2)
mean(run_tmp2 %in% run_tmp1)

# So, did the listInteractionSubsetFASTA screw up packaging sequences into datasets
match_names = match(names(Full_IntAct1$forSLIMFinder_Ready$interaction_subset), names(Full_IntAct2$forSLIMFinder_Ready$interaction_subset))
# names that did match - the same number - 0.997133
mean(!is.na(match_names))
matching_datasets = sapply(1:length(match_names), function(ind){
    #mean(Full_IntAct1$forSLIMFinder_Ready$interaction_subset[[ind]]$ids_set1 %in%
    #          Full_IntAct2$forSLIMFinder_Ready$interaction_subset[[match_names[ind]]]$ids_set1)
    #mean(Full_IntAct1$forSLIMFinder_Ready$interaction_subset[[ind]]$ids_set2 %in%
    #       Full_IntAct2$forSLIMFinder_Ready$interaction_subset[[match_names[ind]]]$ids_set2)
    #mean(names(Full_IntAct1$forSLIMFinder_Ready$fasta_subset_list[[ind]]) %in%
    #       names(Full_IntAct2$forSLIMFinder_Ready$fasta_subset_list[[match_names[ind]]])) == 1
    # invert comparison
    mean(names(Full_IntAct2$forSLIMFinder_Ready$fasta_subset_list[[match_names[ind]]]) %in%
             names(Full_IntAct1$forSLIMFinder_Ready$fasta_subset_list[[ind]])) == 1
})
mean(matching_datasets, na.rm = T)
matching_datasets_TRUE = sapply(matching_datasets, isTRUE)
mean(matching_datasets_TRUE)
# have a look at a couple of examples:
matching_datasets[!matching_datasets_TRUE][[10]]
Full_IntAct1$forSLIMFinder$interaction_subset[!matching_datasets_TRUE][[10]]$ids_set1
Full_IntAct2$forSLIMFinder$interaction_subset[match_names[which(!matching_datasets_TRUE)]][[10]]$ids_set1

all.equal(Full_IntAct1$forSLIMFinder$interaction_subset[!matching_datasets_TRUE][10]$`interactors_of.A0FGR8:P0C739.`$combined_MITAB$Interaction_identifiers,
          Full_IntAct2$forSLIMFinder$interaction_subset[match_names[which(!matching_datasets_TRUE)]][10]$`interactors_of.A0FGR8:P0C739.`$combined_MITAB$Interaction_identifiers)

############################  replicating old pipeline  RESULTS
intact3 = fread("/Users/vk7/Desktop/ebi_projects/viral_project/qslimfinder.Full_IntAct3.FALSE/result/occurence.txt", stringsAsFactors = F)
intact_old = fread("/Users/vk7/Desktop/ebi_projects/viral_project/Full_IntAct2/result/occurence.txt", stringsAsFactors = F)
intact_very_old = fread("/Users/vk7/Desktop/ebi_projects/viral_project/Full_IntAct/result/occurence.txt", stringsAsFactors = F)

old_datasets = unique(intact_old$Dataset)
very_old_datasets = unique(intact_very_old$Dataset)
new_datasets = unique(intact3$Dataset)

mean(very_old_datasets %in% old_datasets)
mean(old_datasets %in% very_old_datasets)
mean(old_datasets %in% new_datasets)
mean(new_datasets %in% old_datasets)
mean(new_datasets %in% very_old_datasets)
### conclusion
# 1. all datasets are constructed from the original interaction data
# 2. QSLIMFinder datasets (seed + query) are the same, except ~30, but 20% contain at least one different sequence per (seed + query). Most of these are due to new run containing more interacting partners for 16.55814% for seeds (the default is to use download new release if available). 5.7% though contained interactions unique to old run. Why? no idea. Only 0.15% of human-human and 0.48% of human-viral interactions are unique to old PPI dataset
# 3. bash commands are almost identical (0.997133) driven by difference in 30 datasets
# 4. Results are very different between the old and the new run: 2820 datasets returned results in February but only 955 in October (new pipeline is 2890)


mod_new_datasets = gsub("interactors_of\\.", "", new_datasets)
mod_new_datasets = strsplit(mod_new_datasets, "\\.")

############################ old pipeline vs new pipeline
mean(Full_IntAct2$forSLIMFinder_file_list[,paste0(interactors_of, QSLIMFinder_query, sequences)] %in%
         Full_IntAct3$forSLIMFinder_file_list[,paste0(interactors_of, QSLIMFinder_query, sequences)])

all.equal(Full_IntAct2$forSLIMFinder$fasta_subset_list[order(names(Full_IntAct2$forSLIMFinder$fasta_subset_list))],
          Full_IntAct3$forSLIMFinder$fasta_subset_list[order(names(Full_IntAct3$forSLIMFinder$fasta_subset_list))])

mean(names(Full_IntAct2$forSLIMFinder$fasta_subset_list) %in% names(Full_IntAct3$forSLIMFinder$fasta_subset_list))
mean(names(Full_IntAct3$forSLIMFinder$fasta_subset_list) %in% names(Full_IntAct2$forSLIMFinder$fasta_subset_list))

mean(names(Full_IntAct2$forSLIMFinder_Ready$fasta_subset_list) %in% names(Full_IntAct3$forSLIMFinder_Ready$fasta_subset_list))
mean(names(Full_IntAct3$forSLIMFinder_Ready$fasta_subset_list) %in% names(Full_IntAct2$forSLIMFinder_Ready$fasta_subset_list))

all.equal(Full_IntAct2$forSLIMFinder_Ready$fasta_subset_list[order(names(Full_IntAct2$forSLIMFinder_Ready$fasta_subset_list))],
          Full_IntAct3$forSLIMFinder_Ready$fasta_subset_list[order(names(Full_IntAct3$forSLIMFinder_Ready$fasta_subset_list))])

all.equal(Full_IntAct2$forSLIMFinder_file_list[,.(interactors_of, QSLIMFinder_query, sequences)], Full_IntAct3$forSLIMFinder_file_list[,.(interactors_of, QSLIMFinder_query, sequences)])

Full_IntAct2$options
Full_IntAct3$options

all.equal(Full_IntAct2$QSLIMFinder_main_result, Full_IntAct3$QSLIMFinder_main_result)


# I started with identical datasets
all.equal(Full_IntAct2$forSLIMFinder_file_list[,.(interactors_of, QSLIMFinder_query, sequences)],
          Full_IntAct3$forSLIMFinder_file_list[,.(interactors_of, QSLIMFinder_query, sequences)])
# run qslimfinder using these options:
Full_IntAct2$options
Full_IntAct3$options
# but ended up with different results (including jobs that started but found nothing)
mean(unique(Full_IntAct2$QSLIMFinder_main_result$Dataset) %in%
         Full_IntAct3$QSLIMFinder_main_result$Dataset)
mean(unique(Full_IntAct3$QSLIMFinder_main_result$Dataset) %in%
         Full_IntAct2$QSLIMFinder_main_result$Dataset)

# but ended up with different results (including only jobs found something)
mean(unique(Full_IntAct2$QSLIMFinder_main_result[!is.na(IC),Dataset]) %in%
         Full_IntAct3$QSLIMFinder_main_result[!is.na(IC),Dataset])
mean(unique(Full_IntAct3$QSLIMFinder_main_result[!is.na(IC),Dataset]) %in%
         Full_IntAct2$QSLIMFinder_main_result[!is.na(IC),Dataset])
# jobs that found nothing
mean(unique(Full_IntAct2$QSLIMFinder_main_result[is.na(IC),Dataset]) %in%
         Full_IntAct3$QSLIMFinder_main_result[is.na(IC),Dataset])
mean(unique(Full_IntAct3$QSLIMFinder_main_result[is.na(IC),Dataset]) %in%
         Full_IntAct2$QSLIMFinder_main_result[is.na(IC),Dataset])

# sequences in fasta subset checker

