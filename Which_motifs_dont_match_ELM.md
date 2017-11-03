# DomainResultsTable: viral protein - human domain pairs
Vitalii Kleshchevnikov  
17/08/2017  



# Load results


```r
name = "./DomainResultsTable_Vidal20171019.RData"
if(!file.exists(name)){
    load("./processed_data_files/what_we_find_VS_ELM_clust20171019.RData")
    occurence_QSLIMFinder = fread("./SLIMFinder_Vidal/result/occurence.txt")
    patterns_QSLIMFinder = fread("./SLIMFinder_Vidal/result/main_result.txt")
    comparimotif_wdb = fread("./SLIMFinder_Vidal/result/comparimotif.compare.tdt")
    rm(list = ls()[!ls() %in% c("printTable","XYZ.p.adjust","res_count", "sequential_filter", "name", "occurence_QSLIMFinder", "comparimotif_wdb", "patterns_QSLIMFinder")])
    save.image(file = name)
} else {
    load(name)
}
doman_viral_pairs = T # if false - add human proteins containing domains
motifs = T # based on Vidal's data
```

## Empirical p-value for seing a domain a number of times
### What is the chance of randomly seeing any domain the observed number of times among all proteins that interact with a specific viral protein


```r
fdr_motifs = 1
fdr_pval_thresh = 1
X = printMotifDomainTable(res_count, doman_viral_pairs = doman_viral_pairs, motifs = motifs, destfile = "./domains_fdr_corrected_only_with_significant_motifs_empirical_p_value.tsv",
               only_with_motifs = T, fdr_motifs = fdr_motifs, fdr_pval_thresh = fdr_pval_thresh,
               occurence_QSLIMFinder = occurence_QSLIMFinder,
               comparimotif_wdb = comparimotif_wdb,
               patterns_QSLIMFinder = patterns_QSLIMFinder, print_table = F)
```

## Which motifs don't match ELM


```r
X4 = X[is.na(motif_mapping)]
X4[, .(NviralProt = data.table::uniqueN(viral_interactor)), by = domain_name]
```

```
##                                                    domain_name NviralProt
##  1:                                     Zinc finger, RING-type          2
##  2:                            Zinc finger, RING/FYVE/PHD-type          2
##  3:                                     Armadillo-like helical          2
##  4:                                        Armadillo-type fold          2
##  5:        P-loop containing nucleoside triphosphate hydrolase          1
##  6:                                      Zinc finger, LIM-type          1
##  7:                                   Immunoglobulin-like fold          1
##  8:                         Armadillo repeat-containing domain          1
##  9:                                         AAA+ ATPase domain          1
## 10:                                                 SAP domain          1
## 11:                                      Zinc finger, MIZ-type          1
## 12:                                               PINIT domain          1
## 13:                                     B-box-type zinc finger          2
## 14:                                          B30.2/SPRY domain          2
## 15:                                                SPRY domain          2
## 16:                             Butyrophylin-like, SPRY domain          2
## 17:                Concanavalin A-like lectin/glucanase domain          2
## 18:                                                  TRAF-like          1
## 19:                                                 PDZ domain          1
## 20:                           Importin-beta, N-terminal domain          1
## 21:                              Rho GTPase activation protein          1
## 22:                                       Immunoglobulin E-set          1
## 23:                                  FAD/NAD(P)-binding domain          1
## 24:                                      Zinc finger C2H2-type          1
## 25:                                      Protein kinase domain          1
## 26: Serine-threonine/tyrosine-protein kinase, catalytic domain          1
## 27:                                 Protein kinase-like domain          1
## 28:                  Tyrosine-protein kinase, catalytic domain          1
## 29:        S-adenosyl-L-methionine-dependent methyltransferase          1
## 30:                                                 SH2 domain          1
## 31:                                             PH domain-like          1
## 32:                                   GroEL-like apical domain          1
## 33:                               GroEL-like equatorial domain          1
## 34:                                                  C2 domain          1
## 35:                                   Calponin homology domain          1
##                                                    domain_name NviralProt
```

```r
X4[, .(NviralProt = data.table::uniqueN(viral_interactor)), by = Motif_Pattern]
```

```
##    Motif_Pattern NviralProt
## 1:      I..QV..T          1
## 2:          G.AG          1
## 3:      Q..TH..P          1
## 4:     W..G.G..Y          1
## 5:       P..N.QI          1
```

```r
DT::datatable(X4, escape = FALSE)
```

<!--html_preserve--><div id="htmlwidget-0171db64defae6d2ee11" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-0171db64defae6d2ee11">{"x":{"filter":"none","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44"],["IPR001841","IPR013083","IPR011989","IPR016024","IPR027417","IPR001781","IPR013783","IPR006911","IPR011989","IPR016024","IPR003593","IPR013083","IPR003034","IPR004181","IPR023321","IPR000315","IPR001870","IPR003877","IPR003879","IPR013320","IPR008974","IPR001478","IPR001841","IPR001494","IPR008936","IPR014756","IPR023753","IPR013087","IPR000719","IPR001245","IPR011009","IPR020635","IPR000315","IPR001870","IPR003877","IPR003879","IPR013320","IPR029063","IPR000980","IPR011993","IPR027409","IPR027413","IPR000008","IPR001715"],["Zinc finger, RING-type","Zinc finger, RING/FYVE/PHD-type","Armadillo-like helical","Armadillo-type fold","P-loop containing nucleoside triphosphate hydrolase","Zinc finger, LIM-type","Immunoglobulin-like fold","Armadillo repeat-containing domain","Armadillo-like helical","Armadillo-type fold","AAA+ ATPase domain","Zinc finger, RING/FYVE/PHD-type","SAP domain","Zinc finger, MIZ-type","PINIT domain","B-box-type zinc finger","B30.2/SPRY domain","SPRY domain","Butyrophylin-like, SPRY domain","Concanavalin A-like lectin/glucanase domain","TRAF-like","PDZ domain","Zinc finger, RING-type","Importin-beta, N-terminal domain","Rho GTPase activation protein","Immunoglobulin E-set","FAD/NAD(P)-binding domain","Zinc finger C2H2-type","Protein kinase domain","Serine-threonine/tyrosine-protein kinase, catalytic domain","Protein kinase-like domain","Tyrosine-protein kinase, catalytic domain","B-box-type zinc finger","B30.2/SPRY domain","SPRY domain","Butyrophylin-like, SPRY domain","Concanavalin A-like lectin/glucanase domain","S-adenosyl-L-methionine-dependent methyltransferase","SH2 domain","PH domain-like","GroEL-like apical domain","GroEL-like equatorial domain","C2 domain","Calponin homology domain"],["P03188","P03188","P03225","P03225","P03225","P03225","P03225","Q3KSQ2","Q3KSQ2","Q3KSQ2","P03225","P03225","P03116","P03116","P03116","P03188","P03188","P03188","P03188","P03188","P03188","P88980","P03225","P03225","P03225","P03225","P03225","P03225","P03225","P03225","P03225","P03225","P03225","P03225","P03225","P03225","P03225","P03225","P03225","P03225","P03225","P03225","P03225","P03225"],[0.000243859406798757,0.000243859406798757,0.00552849510611244,0.00552849510611244,0.00552849510611244,0.010784892409309,0.010784892409309,0.0130653186265909,0.0130653186265909,0.0130653186265909,0.0221336726712965,0.0221336726712965,0.0248104514326335,0.0248104514326335,0.0248104514326335,0.0355475694800765,0.0355475694800765,0.0355475694800765,0.0355475694800765,0.0355475694800765,0.0355475694800765,0.0420910073216716,0.0505135269345866,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803,0.147922622144803],[0.0302697448146913,0.0302697448146913,0.222646831058608,0.222646831058608,0.222646831058608,0.340003383020812,0.340003383020812,0.365058888194971,0.365058888194971,0.365058888194971,0.453073612873889,0.453073612873889,0.474894396547349,0.474894396547349,0.474894396547349,0.547690373090322,0.547690373090322,0.547690373090322,0.547690373090322,0.547690373090322,0.547690373090322,0.577655493669289,0.631522815389379,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678,0.985146190816678],[4,4,6,6,6,5,5,2,2,2,4,4,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],["<a href='https://www.ebi.ac.uk/interpro/entry/IPR001841'>IPR001841<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR013083'>IPR013083<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR011989'>IPR011989<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR016024'>IPR016024<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR027417'>IPR027417<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR001781'>IPR001781<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR013783'>IPR013783<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR006911'>IPR006911<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR011989'>IPR011989<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR016024'>IPR016024<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR003593'>IPR003593<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR013083'>IPR013083<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR003034'>IPR003034<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR004181'>IPR004181<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR023321'>IPR023321<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR000315'>IPR000315<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR001870'>IPR001870<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR003877'>IPR003877<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR003879'>IPR003879<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR013320'>IPR013320<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR008974'>IPR008974<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR001478'>IPR001478<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR001841'>IPR001841<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR001494'>IPR001494<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR008936'>IPR008936<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR014756'>IPR014756<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR023753'>IPR023753<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR013087'>IPR013087<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR000719'>IPR000719<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR001245'>IPR001245<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR011009'>IPR011009<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR020635'>IPR020635<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR000315'>IPR000315<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR001870'>IPR001870<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR003877'>IPR003877<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR003879'>IPR003879<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR013320'>IPR013320<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR029063'>IPR029063<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR000980'>IPR000980<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR011993'>IPR011993<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR027409'>IPR027409<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR027413'>IPR027413<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR000008'>IPR000008<\/a>","<a href='https://www.ebi.ac.uk/interpro/entry/IPR001715'>IPR001715<\/a>"],["<a href='http://www.uniprot.org/uniprot/P03188'>P03188<\/a>","<a href='http://www.uniprot.org/uniprot/P03188'>P03188<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/Q3KSQ2'>Q3KSQ2<\/a>","<a href='http://www.uniprot.org/uniprot/Q3KSQ2'>Q3KSQ2<\/a>","<a href='http://www.uniprot.org/uniprot/Q3KSQ2'>Q3KSQ2<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03116'>P03116<\/a>","<a href='http://www.uniprot.org/uniprot/P03116'>P03116<\/a>","<a href='http://www.uniprot.org/uniprot/P03116'>P03116<\/a>","<a href='http://www.uniprot.org/uniprot/P03188'>P03188<\/a>","<a href='http://www.uniprot.org/uniprot/P03188'>P03188<\/a>","<a href='http://www.uniprot.org/uniprot/P03188'>P03188<\/a>","<a href='http://www.uniprot.org/uniprot/P03188'>P03188<\/a>","<a href='http://www.uniprot.org/uniprot/P03188'>P03188<\/a>","<a href='http://www.uniprot.org/uniprot/P03188'>P03188<\/a>","<a href='http://www.uniprot.org/uniprot/P88980'>P88980<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>","<a href='http://www.uniprot.org/uniprot/P03225'>P03225<\/a>"],[4,4,6,6,6,5,5,2,2,2,4,4,2,2,2,2,2,2,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],[10,10,78,78,78,78,78,4,4,4,78,78,7,7,7,10,10,10,10,10,10,12,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78,78],[308,862,428,444,1055,93,828,9,428,444,157,862,35,11,6,86,89,83,56,194,26,208,308,23,105,127,58,712,637,170,698,91,86,89,83,56,194,129,151,562,30,31,171,111],[15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940,15940],["I..QV..T","I..QV..T","G.AG","G.AG","G.AG","G.AG","G.AG","Q..TH..P","Q..TH..P","Q..TH..P","G.AG","G.AG","W..G.G..Y","W..G.G..Y","W..G.G..Y","I..QV..T","I..QV..T","I..QV..T","I..QV..T","I..QV..T","I..QV..T","P..N.QI","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG","G.AG"],[0.091,0.091,0.07,0.07,0.07,0.07,0.07,0.049,0.049,0.049,0.07,0.07,0.16,0.16,0.16,0.091,0.091,0.091,0.091,0.091,0.091,0.2,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07],["IEEQVNKT","IEEQVNKT","GAAG","GAAG","GAAG","GAAG","GAAG","QVPTHWPP","QVPTHWPP","QVPTHWPP","GAAG","GAAG","WDSGLGCSY","WDSGLGCSY","WDSGLGCSY","IEEQVNKT","IEEQVNKT","IEEQVNKT","IEEQVNKT","IEEQVNKT","IEEQVNKT","PGENYQI","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG","GAAG"],[4,4,3,3,3,3,3,4,4,4,3,3,4,4,4,4,4,4,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3],[87,87,205,205,205,205,205,36,36,36,205,205,44,44,44,87,87,87,87,87,87,32,205,205,205,205,205,205,205,205,205,205,205,205,205,205,205,205,205,205,205,205,205,205],[5,5,47,47,47,47,47,4,4,4,47,47,3,3,3,5,5,5,5,5,5,3,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47,47],[66,66,107,107,107,107,107,30,30,30,107,107,28,28,28,66,66,66,66,66,66,26,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107],[5,5,26,26,26,26,26,4,4,4,26,26,3,3,3,5,5,5,5,5,5,3,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26,26],["interactors_of.Q8IUQ4.P03188","interactors_of.Q8IUQ4.P03188","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q96D09.Q3KSQ2","interactors_of.Q96D09.Q3KSQ2","interactors_of.Q96D09.Q3KSQ2","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.P63165.P03116","interactors_of.P63165.P03116","interactors_of.P63165.P03116","interactors_of.Q8IUQ4.P03188","interactors_of.Q8IUQ4.P03188","interactors_of.Q8IUQ4.P03188","interactors_of.Q8IUQ4.P03188","interactors_of.Q8IUQ4.P03188","interactors_of.Q8IUQ4.P03188","interactors_of.O94972.P88980","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225","interactors_of.Q9UBB9.P03225"],[null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null],[null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null,null]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>human_domain<\/th>\n      <th>domain_name<\/th>\n      <th>viral_interactor<\/th>\n      <th>p.value<\/th>\n      <th>fdr_pval<\/th>\n      <th>observed_statistic<\/th>\n      <th>human_domain_url<\/th>\n      <th>viral_interactor_url<\/th>\n      <th>domain_count_per_viral_interactor<\/th>\n      <th>viral_interactor_degree<\/th>\n      <th>total_domain_count<\/th>\n      <th>total_background_proteins<\/th>\n      <th>Motif_Pattern<\/th>\n      <th>Motif_pval<\/th>\n      <th>Motif_Match<\/th>\n      <th>Motif_IC<\/th>\n      <th>Motif_SeqNum<\/th>\n      <th>Motif_OccNum<\/th>\n      <th>Motif_UPNum<\/th>\n      <th>Motif_UPoccNum<\/th>\n      <th>Motif_Dataset<\/th>\n      <th>motif_mapping<\/th>\n      <th>motif..............................................pattern_mapping<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[4,5,6,9,10,11,12,14,16,17,18,19,20]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->
