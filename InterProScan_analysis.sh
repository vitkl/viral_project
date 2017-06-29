ssh vitalii@ebi-cli
bsub rsync vitalii@petsalakigroup02-ml.windows.ebi.ac.uk:~/Desktop/vitalii /nfs/research/petsalaki/users/vitalii -r -av
## with match lookup 
bsub rsync /nfs/research/petsalaki/users/vitalii /hps/nobackup/research/petsalaki/users/ -r -av
bsub -n 16 -q research-rh7 -M 16384 -R "rusage[mem=16384]" /hps/nobackup/research/petsalaki/users/vitalii/my_interproscan/interproscan-5.24-63.0/interproscan.sh -i /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/data_files/all_human_viral_proteins.fasta -f tsv -iprlookup -goterms -b /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/processed_data_files/all_human_viral_protein_domains.tsv
bsub rsync /hps/nobackup/research/petsalaki/users/vitalii /nfs/research/petsalaki/users/ -r -av

bsub -Is $SHELL 
rsync /nfs/research/petsalaki/users/vitalii/vitalii vitalii@petsalakigroup02-ml.windows.ebi.ac.uk:~/Desktop/ -r -av
exit

## without match lookup 
bsub -n 16 -q research-rh7 -M 16384 -R "rusage[mem=16384]" /hps/nobackup/research/petsalaki/users/vitalii/my_interproscan/interproscan-5.24-63.0/interproscan.sh -i /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/data_files/all_human_viral_proteins.fasta -dp -f tsv -iprlookup -goterms -b /hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/processed_data_files/all_human_viral_protein_domains_nolookup.tsv
bsub rsync /hps/nobackup/research/petsalaki/users/vitalii /nfs/research/petsalaki/users/ -r -av