# Identify if the domains we have identified tend to bind linear motifs (the domains we have identified are enriched in domains known to bind linear motifs)
Vitalii Kleshchevnikov  
17/08/2017  



## Setting up how many cores to use for analysis and fold-enrichment vs freqency


```r
# cores to use for multiprocessing: is the number of cores is small -> we are not on a cluster -> use all; if the number of cores is large -> we are on cluster -> use only 15 or 7 or 31 cores (n requested-1)
if(detectCores() <= 4) cores_to_use = detectCores() - 1
if(detectCores() > 4) cores_to_use = 31

# how many permutations?
N_permut = 30000
```

# How many of the domains identified using our approach are in the ELM database (as compared to the population of all domains we tested)?

## Calculate empirical pvalues using either frequency of a domain among interacting partners of viral protein or Fisher test p-value


```r
data = fread("./processed_data_files/viral_human_net_w_domains", sep = "\t", stringsAsFactors = F)
# frequency function: set up standard parameters
permuteFrequency = function(data, select_nodes = NULL, also_permuteYZ = F){
    res = permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human, # first set of interacting pairs (XY) that are to be permuted
                          associations2test = IDs_interactor_viral ~ IDs_domain_human, # set of interacting pairs to be tested (XZ), YZ interactions are assumed
                          node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
                                           IDs_domain_human ~ domain_count,
                                           IDs_interactor_viral + IDs_domain_human ~ domain_frequency_per_IDs_interactor_viral),
                          data = data,
                          statistic = IDs_interactor_viral + IDs_domain_human ~ .N / IDs_interactor_viral_degree,
                          select_nodes = select_nodes,
                          N = N_permut,
                          cores = cores_to_use, seed = 2, also_permuteYZ = also_permuteYZ)
    return(res)
}

# frequency: all proteins and domains - # permute IDs_interactor_viral ~ IDs_interactor_human
time = proc.time()
res = permuteFrequency(data, select_nodes = NULL)
proc.time() - time
```

```
##     user   system  elapsed 
##   11.931    1.704 4215.334
```

```r
plot(res, main = "frequency: all proteins and domains")
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/calculate_pvals_freq-1.png)<!-- -->

```r
# frequency: no low background count domains
res_low_back = permuteFrequency(data, select_nodes = IDs_domain_human ~ domain_count >= 3)
#res_low_back_alt = permuteFrequency(data[domain_count >= 3,])
#all.equal(res_low_back$data_with_pval[complete.cases(res_low_back$data_with_pval),p.value], res_low_back_alt$data_with_pval[complete.cases(res_low_back_alt$data_with_pval),p.value])
plot(res_low_back, main = "frequency: no low background count domains (>= 3)")
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/calculate_pvals_freq-2.png)<!-- -->

```r
# frequency: not considering (fixing interactions, degree of every node in the network stays the same, but only high degree proteins are taken into account, equivalent to  permuting only interactions of protein with the degree higher than 1) viral proteins with the degree of 1 - removing viral proteins with the degree of 1
res_low_deg = permuteFrequency(data, select_nodes = IDs_interactor_viral ~ IDs_interactor_viral_degree >= 2)
#res_low_deg_alt = permuteFrequency(data[IDs_interactor_viral_degree >= 2,])
#all.equal(res_low_deg$data_with_pval[complete.cases(res_low_deg$data_with_pval),p.value], res_low_deg_alt$data_with_pval[complete.cases(res_low_deg_alt$data_with_pval),p.value])

plot(res_low_deg, main = "frequency: not considering viral proteins with the degree of 1")
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/calculate_pvals_freq-3.png)<!-- -->

```r
# frequency: BOTH no low background count domains AND removing viral proteins with the degree of 1

res_low_deg_back = permuteFrequency(data, select_nodes = list(IDs_domain_human ~ domain_count >= 3,
                                                              IDs_interactor_viral ~ IDs_interactor_viral_degree >= 2))
#res_low_deg_back_alt = permuteFrequency(data[domain_count >= 3 & IDs_interactor_viral_degree >= 2,])
#all.equal(res_low_deg_back$data_with_pval[complete.cases(res_low_deg_back$data_with_pval),p.value], res_low_deg_back_alt$data_with_pval[complete.cases(res_low_deg_back_alt$data_with_pval),p.value])
plot(res_low_deg_back, main = "frequency: no low background count domains (>= 3)\nAND not considering viral proteins with the degree of 1")
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/calculate_pvals_freq-4.png)<!-- -->

```r
save(res, file="./processed_data_files/what_we_find_VS_ELM_output_freq.RData")
```


```r
# count function: set up standard parameters
permuteCount = function(data, select_nodes = NULL, also_permuteYZ = F){
    res = permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human, # first set of interacting pairs (XY) that are to be permuted
                          associations2test = IDs_interactor_viral ~ IDs_domain_human, # set of interacting pairs to be tested (XZ), YZ interactions are assumed
                          node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
                                           IDs_domain_human ~ domain_count,
                                           IDs_interactor_viral + IDs_domain_human ~ domain_frequency_per_IDs_interactor_viral),
                          data = data,
                          statistic = IDs_interactor_viral + IDs_domain_human ~ .N,
                          select_nodes = select_nodes,
                          N = N_permut,
                          cores = cores_to_use, seed = 2, also_permuteYZ = also_permuteYZ)
    return(res)
}
# count: all proteins and domains - # permute IDs_interactor_viral ~ IDs_interactor_human
time = proc.time()
res_count = permuteCount(data, select_nodes = NULL)
proc.time() - time
```

```
##     user   system  elapsed 
##    8.040    2.299 2450.669
```

```r
plot(res_count, main = "count: all proteins and domains")
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/calculate_pvals_count-1.png)<!-- -->

```r
save(list = ls(), file="./processed_data_files/what_we_find_VS_ELM_clust.RData")
```


```r
# Fisher test: set up standard parameters
permuteFisherTestPval = function(data, select_nodes = NULL, also_permuteYZ = F, N = N_permut){
    resFISHER = permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human, # first set of interacting pairs (XY) that are to be permuted
                                associations2test = IDs_interactor_viral ~ IDs_domain_human, # set of interacting pairs to be tested (XZ), YZ interactions are assumed
                                node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree, # attribute of X
                                                 IDs_domain_human ~ domain_count + N_prot_w_interactors, # attributes of Z
                                                 IDs_interactor_viral + IDs_domain_human ~ domain_count_per_IDs_interactor_viral), # attribute of both X and Z
                                data = data, # data.table containing data
                                statistic = IDs_interactor_viral + IDs_domain_human ~ fisher.test(
                                    matrix(c(domain_count_per_IDs_interactor_viral[1], 
                                             IDs_interactor_viral_degree[1] - domain_count_per_IDs_interactor_viral[1],
                                             domain_count[1] - domain_count_per_IDs_interactor_viral[1], 
                                             N_prot_w_interactors[1] - domain_count[1] - IDs_interactor_viral_degree[1] + domain_count_per_IDs_interactor_viral[1]),
                                           2,2), 
                                    alternative = "greater", conf.int = F)$p.value, # formula to calculate statisic by evaluating right-hand-side expression for each X and Z pair, right-hand-side expression is what is normally put in j in data.table DT[i, j, by], left-hand-side expression contains column names of X and Z which are used in by in data.table
                                select_nodes = select_nodes, # select a subset of the data, only nodes 
                                N = N, # number of permutations
                                cores = cores_to_use, seed = 1, also_permuteYZ = also_permuteYZ)
    # permutationPval returns the number of cases when permuted statitic is higher than the observed statistic (right tail of the distribution), in this case we are interested in the reverse - the lower tail, when p-values from permuted distribution that are lower than the observed p-value
    resFISHER$data_with_pval[, p.value := 1 - p.value]
    return(resFISHER)
}

# contingency matrix:
matrix(c("domain_count_per_IDs_interactor_viral[1]", 
         "IDs_interactor_viral_degree[1] - domain_count_per_IDs_interactor_viral[1]",
         "domain_count[1] - domain_count_per_IDs_interactor_viral[1]", 
         "N_prot_w_interactors[1] - domain_count[1] - IDs_interactor_viral_degree[1] + domain_count_per_IDs_interactor_viral[1]"),
       2,2)
```

```
##      [,1]                                                                       
## [1,] "domain_count_per_IDs_interactor_viral[1]"                                 
## [2,] "IDs_interactor_viral_degree[1] - domain_count_per_IDs_interactor_viral[1]"
##      [,2]                                                                                                                   
## [1,] "domain_count[1] - domain_count_per_IDs_interactor_viral[1]"                                                           
## [2,] "N_prot_w_interactors[1] - domain_count[1] - IDs_interactor_viral_degree[1] + domain_count_per_IDs_interactor_viral[1]"
```

```r
# permute IDs_interactor_viral ~ IDs_interactor_human
time = proc.time()
resFISHERpval = permuteFisherTestPval(data, select_nodes = NULL)
proc.time() - time
```

```
##     user   system  elapsed 
##   17.527    2.288 3176.940
```

```r
plot(resFISHERpval, main = "Fisher test P-value: all proteins and domains")
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/calculate_pvals_FisherTest-1.png)<!-- -->

```r
# Fisher test: set up standard parameters
permuteFisherTestOdds = function(data, select_nodes = NULL, also_permuteYZ = F){
    resFISHER = permutationPval(interactions2permute = IDs_interactor_viral ~ IDs_interactor_human, # first set of interacting pairs (XY) that are to be permuted
                                associations2test = IDs_interactor_viral ~ IDs_domain_human, # set of interacting pairs to be tested (XZ), YZ interactions are assumed
                                node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree, # attribute of X
                                                 IDs_domain_human ~ domain_count + N_prot_w_interactors, # attributes of Z
                                                 IDs_interactor_viral + IDs_domain_human ~ domain_count_per_IDs_interactor_viral), # attribute of both X and Z
                                data = data, # data.table containing data
                                statistic = IDs_interactor_viral + IDs_domain_human ~ fisher.test(
                                    matrix(c(domain_count_per_IDs_interactor_viral[1], 
                                             IDs_interactor_viral_degree[1] - domain_count_per_IDs_interactor_viral[1],
                                             domain_count[1] - domain_count_per_IDs_interactor_viral[1], 
                                             N_prot_w_interactors[1] - domain_count[1] - IDs_interactor_viral_degree[1] + domain_count_per_IDs_interactor_viral[1]),
                                           2,2), 
                                    alternative = "greater", conf.int = F)$estimate, # formula to calculate statisic by evaluating right-hand-side expression for each X and Z pair, right-hand-side expression is what is normally put in j in data.table DT[i, j, by], left-hand-side expression contains column names of X and Z which are used in by in data.table
                                select_nodes = select_nodes, # select a subset of the data, only nodes 
                                N = N_permut, # number of permutations
                                cores = cores_to_use, seed = 1, also_permuteYZ = also_permuteYZ)
    # permutationPval returns the number of cases when permuted statitic is higher than the observed statistic (right tail of the distribution), in this case we are interested in the reverse - the lower tail, when p-values from permuted distribution that are lower than the observed p-value
    resFISHER$data_with_pval[, p.value := 1 - p.value]
    return(resFISHER)
}

# permute IDs_interactor_viral ~ IDs_interactor_human: odds ratio
time = proc.time()
resFISHERodds = permuteFisherTestOdds(data, select_nodes = NULL)
proc.time() - time
```

```
##     user   system  elapsed 
##   17.305    2.588 2978.267
```

```r
plot(resFISHERodds, main = "Fisher test odds ratio: all proteins and domains")
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/calculate_pvals_FisherTest-2.png)<!-- -->

```r
# permute IDs_interactor_viral ~ IDs_interactor_human: odds ratio, remove all small values
time = proc.time()
resFISHERoddsBig = permuteFisherTestOdds(data, select_nodes = list(IDs_domain_human ~ domain_count >= 15,
                                                                   IDs_interactor_viral ~ IDs_interactor_viral_degree >= 15))
proc.time() - time
```

```
##     user   system  elapsed 
##    9.276    2.140 2304.369
```

```r
plot(resFISHERoddsBig, main = "Fisher test odds ratio: removed small (>= 15)\ndomain count and degree")
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/calculate_pvals_FisherTest-3.png)<!-- -->

```r
save(resFISHERodds, file="./processed_data_files/what_we_find_VS_ELM_output_FISHER.RData")

# Just FisherTest result (no permutation)
justFisherdata = copy(resFISHERodds$data_with_pval)
justFisherdata[, resFISHERodds_p.value := p.value][, p.value := NULL]
justFisherdata[, odds_ratio := observed_statistic][, observed_statistic := NULL]
justFisherdata[, resFISHERodds_YmissingZ_perX := YmissingZ_perX][, YmissingZ_perX := NULL]
justFisherdata[, resFISHERodds_higher_counts := higher_counts][, higher_counts := NULL]
justFisherdata[, resFISHERodds_not_missing := not_missing][, not_missing := NULL]

resJustFISHER = permuteFisherTestPval(justFisherdata, select_nodes = NULL, N = 1)

resJustFISHER$data_with_pval[, p.value := observed_statistic]
resJustFISHER$data_with_pval[, observed_statistic := odds_ratio]

my.p.adjust = function(res, adj_by = "p.value", ...) {
    if(class(res) != "XYZinteration_XZEmpiricalPval") stop("res should be the output of permutationPval()")
    if("fdr_pval" %in% colnames(res$data_with_pval)) stop("p value has already been fdr-corrected")
    columns = adj_by
    adj_by = formula(paste0("~",adj_by))[[2]]
    nodes = c(res$nodes$nodeX, res$nodes$nodeZ)
    columns = c(nodes, columns)
    temp = unique(res$data_with_pval[, c(columns), with = F])
    temp[, fdr_pval := p.adjust(eval(adj_by), ...)]
    res$data_with_pval = res$data_with_pval[temp, on = c(columns)]
    res
}
resJustFISHER = my.p.adjust(resJustFISHER, method = "fdr")
resJustFISHER005 = copy(resJustFISHER)
resJustFISHER005$data_with_pval = resJustFISHER$data_with_pval[fdr_pval < 0.05,]
resJustFISHER01 = copy(resJustFISHER)
resJustFISHER01$data_with_pval = resJustFISHER$data_with_pval[fdr_pval < 0.1,]
resJustFISHER001 = copy(resJustFISHER)
resJustFISHER001$data_with_pval = resJustFISHER$data_with_pval[fdr_pval < 0.01,]

# Fisher test gives the probability of the domain being enriched over the background P(D)
# If we do just the Fisher test we see that the top-scoring domains (overwhelming majority) are those seen at count 1, we need to account for the probability of seeing different counts while looking at viral proteins individually
# Permutations give the probability of seeing specific domain counts P(count)
# P(cound|D) are the frequencies of each count as is while looking at viral proteins individually
# We want to find P(D|count) = (P(count|D) * P(D)) / P(count)
# P.D_count. = P.count_D. * P.D. / P.count.
bayes = function(resJustFISHER, res_count, P.D. = "p.value", P.count. = "p.value", degree = "IDs_interactor_viral_degree", count = "domain_count_per_IDs_interactor_viral"){
    if(identical(c(resJustFISHER$nodes$nodeX, resJustFISHER$nodes$nodeZ),
                 c(res_count$nodes$nodeX, res_count$nodes$nodeZ))) {
        nodes = c(res_count$nodes$nodeX, res_count$nodes$nodeZ)
    } else stop("resJustFISHER and res_count contain data for different nodes")
    count_data = unique(res_count$data_with_pval[, c(nodes, P.count., degree, count,"observed_statistic"), with = F])
    count_data[, domains_per_prot := .N, by = .(eval(formula(paste0("~",res_count$nodes$nodeX))[[2]]))]
    count_data[, P.count_D. := .N / domains_per_prot, 
                   by = .(eval(formula("~observed_statistic")[[2]]),
                                          eval(formula(paste0("~",res_count$nodes$nodeX))[[2]]))]
    # count_data[, P.count_D. := eval(formula(paste0("~",count))[[2]]) / eval(formula(paste0("~",degree))[[2]])]
    count_data[, P.count. := p.value]
    count_data[, c("p.value", "observed_statistic", "domains_per_prot") := NULL]
    
    just_fisher = copy(resJustFISHER)
    just_fisher$data_with_pval = just_fisher$data_with_pval[count_data, on = nodes]
    just_fisher$data_with_pval[, P.D. := eval(formula(paste0("~",P.D.))[[2]])]
    
    just_fisher$data_with_pval[, P.D. := 1 - P.D.]
    just_fisher$data_with_pval[, P.count_D. := P.count_D.]
    just_fisher$data_with_pval[, P.count. := 1 - P.count.]
    just_fisher$data_with_pval[, P.D_count. := P.count_D. * P.D. / P.count.]
    
    just_fisher$data_with_pval[, p.value := 1 - P.D_count.]
    
    just_fisher
}

resBayes = bayes(resJustFISHER, res_count)
plot(resJustFISHER, IDs_interactor_viral + IDs_domain_human ~ fdr_pval)
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/calculate_pvals_FisherTest-4.png)<!-- -->

```r
plot(resBayes, IDs_interactor_viral + IDs_domain_human ~ P.count_D.)
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/calculate_pvals_FisherTest-5.png)<!-- -->


```r
# frequency function: set up standard parameters
permuteFrequencyRev = function(data, select_nodes = NULL, also_permuteYZ = F){
    res = permutationPval(interactions2permute = IDs_domain_human ~ IDs_interactor_human, # first set of interacting pairs (XY) that are to be permuted
                          associations2test = IDs_domain_human ~ IDs_interactor_viral, # set of interacting pairs to be tested (XZ), YZ interactions are assumed
                          node_attr = list(IDs_interactor_viral ~ IDs_interactor_viral_degree,
                                           IDs_domain_human ~ domain_count,
                                           IDs_interactor_viral + IDs_domain_human ~ domain_frequency_per_IDs_interactor_viral),
                          data = data,
                          statistic = IDs_interactor_viral + IDs_domain_human ~ .N / domain_count,
                          select_nodes = select_nodes,
                          N = N_permut,
                          cores = cores_to_use, seed = 2, also_permuteYZ = also_permuteYZ)
    return(res)
}

# frequency: all proteins and domains - # permute IDs_domain_human ~ IDs_interactor_human
time = proc.time()
resRev = permuteFrequencyRev(data, select_nodes = NULL)
proc.time() - time
```

```
##     user   system  elapsed 
##    9.274    2.337 1729.791
```

```r
plot(resRev, main = "frequency reverse: viral protein frequency among proteins with a domain")
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/calculate_pvals_freq_rev-1.png)<!-- -->

```r
# multiply one way and the reverse-way probabilities
mix = resRev
mix$data_with_pval = resRev$data_with_pval[res$data_with_pval, on = c("IDs_domain_human", "IDs_interactor_viral", "IDs_interactor_human", "domain_type", "N_prot_w_interactors", "domain_count", "domain_frequency", "Taxid_interactor_human", "Taxid_interactor_viral", "IDs_interactor_viral_degree", "domain_frequency_per_IDs_interactor_viral", "fold_enrichment", "IDs_interactor_viral_per_IDs_domain_human", "domain_count_per_IDs_interactor_viral", "IDs_interactor_human_degree", "IDs_domain_human_per_IDs_interactor_viral")]
mix$data_with_pval[, p.value := p.value * i.p.value]
plot(mix, main = "frequency mix: p-val viral protein-domain * p-val domain-viral protein")
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/calculate_pvals_freq_rev-2.png)<!-- -->

```r
save(mix, file="./processed_data_files/what_we_find_VS_ELM_output_freq_mix.RData")
```

## Viral protein degree and the background domain count of top-scoring proteins


```r
# frequency
PermutResult2D(res = res, N = 250) +
    ggtitle("2D-bin plots of 250 top-scoring viral protein - human domain pairs, \n statistic: frequency of a domain among interacting partners of a viral protein (viral-domain)")
PermutResult2D(res = res_low_back, N = 250) +
    ggtitle("2D-bin plots of 250 top-scoring viral protein - human domain pairs, \n statistic: frequency of a domain among interacting partners of a viral protein \n no low background count domains")
PermutResult2D(res = res_low_deg, N = 250) +
    ggtitle("2D-bin plots of 250 top-scoring viral protein - human domain pairs, \n statistic: frequency of a domain among interacting partners of a viral protein \n no viral proteins with the degree of 1")
PermutResult2D(res = res_low_deg_back, N = 250) +
    ggtitle("2D-bin plots of 250 top-scoring viral protein - human domain pairs, \n statistic: frequency of a domain among interacting partners of a viral protein \n no low background count domains \nAND no viral proteins with the degree of 1")

# count
PermutResult2D(res = res_count, N = 250) +
    ggtitle("2D-bin plots of 250 top-scoring viral protein - human domain pairs, \n statistic: count of a domain among interacting partners of a viral protein")

# reverse frequency
PermutResult2D(res = res, N = 250) +
    ggtitle("2D-bin plots of 250 top-scoring viral protein - human domain pairs, \n statistic: viral protein frequency among proteins with a domain")
PermutResult2D(res = mix, N = 250) +
    ggtitle("2D-bin plots of 250 top-scoring viral protein - human domain pairs, \n statistic: frequency p-val viral protein-domain * frequency p-val domain-viral protein")

# Fisher test p-value
PermutResult2D(res = resFISHER, N = 250) + 
    ggtitle("2D-bin plots of 250 top-scoring viral protein - human domain pairs, \n statistic: Fisher test p-value")
```

## Map domains known to interact with linear motifs from ELM to the domains we found


```r
interactiondomains = fread("http://elm.eu.org/interactiondomains.tsv")
interactiondomains[, pfam_id := `Interaction Domain Id`]

domains_known = interactiondomains[, unique(pfam_id)]

InterProScan_domains = readInterProGFF3("../viral_project/processed_data_files/all_human_viral_protein_domains.gff3.gz")
# get InterProID to member database ID mapping
InterPro2memberDB = getInterPro2memberDB(InterProScan_domains)
InterPro2memberDB = InterPro2memberDB[complete.cases(InterPro2memberDB)]
domains_known_mapped = unique(InterPro2memberDB[memberDBID %in% domains_known | InterProID %in% domains_known, InterProID])
domains_not_mapped = unique(domains_known[!(domains_known %in% InterPro2memberDB$memberDBID | domains_known %in% InterPro2memberDB$InterProID)])
```

I did Fisher test to evaluate if the domains that we find are enriched in domains known to interact with linear motifs (from ELM). I have picked some number of viral protein - human domain associations from the top (by p-value). Then I counted how many known domains we have found and did Fisher test. I decided to compare two statistic choices (frequency of a domain among interacting partners of a viral protein or Fisher test p-value) on how many of the known domains we tend to find. Finally, I was choosing different cutoffs (different number of top p-value pairs). 


```r
testEnrichment = function(N, res, rank_by = "p.value", domains_known_mapped, random = F, name = "", decreasing = F){
    if(random) {
        res$data_pval = unique(res$data_with_pval[,c("IDs_interactor_viral", "IDs_domain_human", rank_by, "domain_type", "domain_count", "IDs_interactor_viral_degree"), with = F])
        domains_found = res$data_pval[sample(1:nrow(res$data_with_pval), N), unique(IDs_domain_human)]
    } else {
        res$data_pval = unique(res$data_with_pval[,c("IDs_interactor_viral", "IDs_domain_human", rank_by, "domain_type", "domain_count", "IDs_interactor_viral_degree"), with = F])
        ind = order(unlist(res$data_pval[, c(rank_by), with = F]), decreasing = decreasing)[1:N] 
        domains_found = res$data_pval[ind, unique(IDs_domain_human)]
    }
    
    alldomains = res$data_with_pval[, unique(IDs_domain_human)]
    known = factor(alldomains %in% domains_known_mapped, levels = c("TRUE", "FALSE"))
    found = factor(alldomains %in% domains_found, levels = c("TRUE", "FALSE"))
    table_res = table(known, found)
    
    test = fisher.test(table(known, found), alternative = "greater", conf.int = T)
    
    return(c(pval = test$p.value, odds_ratio = as.vector(test$estimate), count = table_res["TRUE", "TRUE"], name = name))
}
runningTestEnrichment = function(res, name, rank_by = "p.value", decreasing = F){
    enrichment = sapply(Ns, testEnrichment, res = res, 
                        domains_known_mapped = domains_known_mapped, 
                        name = name, rank_by = rank_by, decreasing = decreasing)
    colnames(enrichment) = Ns
    return(enrichment)
}

Ns = seq(25, 500, 25)
# frequency
enrichment = runningTestEnrichment(res, name = "domain frequency among interactors of a viral protein, empirical pval")
enrichment_justfreq = runningTestEnrichment(res, rank_by = "observed_statistic", name = "domain frequency among interactors of a viral protein, frequency")
enrichment_low_back = runningTestEnrichment(res_low_back, name = "domain frequency: no low background, empirical pval")
enrichment_low_deg = runningTestEnrichment(res_low_deg, name = "domain frequency: no degree of 1, empirical pval")
enrichment_low_deg_back = runningTestEnrichment(res_low_deg_back, name = "domain frequency: no degree of 1 AND no low background, empirical pval")

# reverse frequency and mix
enrichmentRev = runningTestEnrichment(resRev, name = "viral protein frequency among proteins with a domain, empirical pval")
enrichmentMix = runningTestEnrichment(mix, name = "mix of viral protein -> domain and domain -> viral protein, empirical pval")

# count
enrichment_count = runningTestEnrichment(res_count, name = "domain count among interactors of a viral protein, empirical pval")

# Fisher test pval
enrichmentFISHERodds = runningTestEnrichment(resFISHERodds, name = "Fisher test odds ratio: domain overrepresentation over the background, empirical pval")
enrichmentFISHERpval = runningTestEnrichment(resFISHERpval, name = "Fisher test pval: domain overrepresentation over the background, empirical pval")
enrichmentFISHER_justodds = runningTestEnrichment(resFISHERodds, rank_by = "observed_statistic", name = "odds ratio: domain overrepresentation over the background, decreasing rank", decreasing = T)
enrichmentFISHER_justoddsIncreasing = runningTestEnrichment(resFISHERodds, rank_by = "observed_statistic", name = "odds ratio: domain overrepresentation over the background, increasing rank", decreasing = F)
enrichmentFISHERpval = runningTestEnrichment(resFISHERpval, name = "Fisher test pval: domain overrepresentation over the background, empirical pval")
enrichmentFISHER_justpval = runningTestEnrichment(resFISHERpval, rank_by = "observed_statistic", name = "Fisher test pval: domain overrepresentation over the background")
enrichmentFISHERoddsBig = runningTestEnrichment(resFISHERoddsBig, name = "Fisher test odds ratio: domain overrepresentation over the background \n(background count and degree >= 15), empirical pval")

enrichmentresJustFISHER005 = runningTestEnrichment(resJustFISHER005, rank_by = "observed_statistic", name = "odds ratio: domain overrepresentation over the background, filtered by Fisher test pval 0.05", decreasing = T)
enrichmentresJustFISHER01 = runningTestEnrichment(resJustFISHER01, rank_by = "observed_statistic", name = "odds ratio: domain overrepresentation over the background, filtered by Fisher test pval 0.1", decreasing = T)
enrichmentresJustFISHER001 = runningTestEnrichment(resJustFISHER001, rank_by = "observed_statistic", name = "odds ratio: domain overrepresentation over the background, filtered by Fisher test pval 0.01", decreasing = T)

# Bayes
enrichmentBayes = runningTestEnrichment(resBayes, name = "Bayesian approach: integrating count and domain probabilities P(D|count) = (P(count|D) * P(D)) / P(count)")


random_domains = function(N = 100, seed = seed, Ns = seq(25, 500, 25)){
    set.seed(seed)
    
    quantiles = c(0.975, 0.75, 0.5, 0.25, 0.025)
    quantile_names = c("97.5% quantile", "75% quantile", "median", "25% quantile", "2.5% quantile")
    
    pval_temp = replicate(N, {
        enrichmentRANDOM = sapply(Ns, testEnrichment, res = res, domains_known_mapped = domains_known_mapped, random = T, name = "N random proteins")[1,]
        names(enrichmentRANDOM) = Ns
        as.numeric(enrichmentRANDOM)
    })
    pval = apply(pval_temp, 1, quantile, probs = quantiles)
    rownames(pval) = quantile_names
    colnames(pval) = Ns
    
    odds_ratio_temp = replicate(N, {
        enrichmentRANDOM = sapply(Ns, testEnrichment, res = res, domains_known_mapped = domains_known_mapped, random = T, name = "N random proteins")[2,]
        names(enrichmentRANDOM) = Ns
        as.numeric(enrichmentRANDOM)
    })
    odds_ratio = apply(odds_ratio_temp, 1, quantile, probs = quantiles)
    rownames(odds_ratio) = quantile_names
    colnames(odds_ratio) = Ns
    
    count_temp = replicate(N, {
        enrichmentRANDOM = sapply(Ns, testEnrichment, res = res, domains_known_mapped = domains_known_mapped, random = T, name = "N random proteins")[3,]
        names(enrichmentRANDOM) = Ns
        as.numeric(enrichmentRANDOM)
    })
    count = apply(count_temp, 1, quantile, probs = quantiles)
    rownames(count) = quantile_names
    colnames(count) = Ns
    
    return(list(pval = pval, odds_ratio = odds_ratio, count = count))
}

enrichmentRANDOM = random_domains(1000, 1)
save(list = ls(), file="./processed_data_files/what_we_find_VS_ELM_clust.RData")
```

## As we include more proteins, the number of known domains we find increases and then levels off (probably because some of the known domains do not interact with viral proteins).


```r
plotEnrichment(enrichment, enrichmentBayes, enrichment_justfreq, enrichment_low_back, enrichment_low_deg, enrichment_low_deg_back,
               enrichment_count, enrichmentRev, enrichmentMix,
               enrichmentFISHERodds, enrichmentFISHER_justodds,
               enrichmentFISHERpval, enrichmentFISHER_justpval, enrichmentFISHERoddsBig,
               random_domains = enrichmentRANDOM, 
               domains_known_mapped = domains_known_mapped, type = "count", plot_type = "l")
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/domains_found-1.png)<!-- -->

## As we include more proteins, the Fisher test odds ratio decreases (we add more stuff that is not known). Odds ratio measures how much more likely are we to find a domain using our procedure if it’s a known domain as compared to if it’s not a known domain.


```r
plotEnrichment(enrichment, enrichment_justfreq, enrichment_low_back, enrichment_low_deg, enrichment_low_deg_back,
               enrichment_count, enrichmentRev, enrichmentMix,
               enrichmentFISHERodds, enrichmentFISHER_justodds,
               enrichmentFISHERpval, enrichmentFISHER_justpval, enrichmentFISHERoddsBig,
               random_domains = enrichmentRANDOM, 
               domains_known_mapped = domains_known_mapped, type = "odds_ratio", plot_type = "l")
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/fisher_odds_ratio-1.png)<!-- -->

## corresponding P-values from the Fisher test 


```r
plotEnrichment(enrichment, enrichment_justfreq, enrichment_low_back, enrichment_low_deg, enrichment_low_deg_back,
               enrichment_count, enrichmentRev, enrichmentMix,
               enrichmentFISHERodds, enrichmentFISHER_justodds,
               enrichmentFISHERpval, enrichmentFISHER_justpval, enrichmentFISHERoddsBig,
               random_domains = enrichmentRANDOM, 
               domains_known_mapped = domains_known_mapped, type = "pval", plot_type = "l")
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/fisher_pval-1.png)<!-- -->

## How many viral proteins are known per each of the domain instances in top 200 protein-domain pairs?


```r
selectTopHits = function(res, N){
    res$data_pval = unique(res$data_with_pval[,.(IDs_interactor_viral, IDs_domain_human, p.value, domain_type, domain_count, IDs_interactor_viral_degree)])
    pairs200pval = res$data_pval[order(p.value, decreasing = F)[1:N], max(p.value)]
    
    restop100 = res
    restop100$data_with_pval = restop100$data_with_pval[p.value <= pairs200pval, ]
    restop100$data_with_pval[, N_viral_per_human_w_domain := length(unique(IDs_interactor_viral)), by = .(IDs_interactor_human, IDs_domain_human)]
    return(restop100)
}
restop100 = selectTopHits(res, N = 250)
plot(restop100, IDs_interactor_human ~ N_viral_per_human_w_domain)
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/top200_n_inter-1.png)<!-- -->

```r
plot(restop100)
```

![](/hps/nobackup/research/petsalaki/users/vitalii/vitalii/viral_project/what_we_find_VS_ELM_clust_files/figure-html/top200_n_inter-2.png)<!-- -->

78 human proteins with enriched domains have 5 or more viral interacting partners.  
15 human proteins with enriched domains have 10 or more viral interacting partners.  

### what are those domains? are they known ELM-interacting domains? which proteins they are in? which viral taxons they interact with?


```r
restop100$data_with_pval[N_viral_per_human_w_domain >= 10, unique(IDs_domain_human)]
```

```
## [1] "IPR016024" "IPR009072" "IPR032454" "IPR032458" "IPR000504"
```

```r
restop100$data_with_pval[N_viral_per_human_w_domain >= 10, unique(IDs_domain_human)] %in% domains_known_mapped
```

```
## [1]  TRUE FALSE FALSE FALSE  TRUE
```

```r
restop100$data_with_pval[N_viral_per_human_w_domain >= 10 & IDs_domain_human == "IPR000504", unique(IDs_interactor_human)]
```

```
## [1] "P09651" "P22626" "Q99729"
```

```r
restop100$data_with_pval[N_viral_per_human_w_domain >= 10 & IDs_domain_human == "IPR000504", unique(Taxid_interactor_viral)]
```

```
##  [1]  10377  28344 211044  10254 380964  88776 130763 128952 680716  10376
## [11]  11097
```

```r
DT::datatable(restop100$data_with_pval[N_viral_per_human_w_domain >= 10,])
```

<!--html_preserve--><div id="htmlwidget-8355d2cc2279d72ee744" style="width:100%;height:auto;" class="datatables html-widget"></div>
<script type="application/json" data-for="htmlwidget-8355d2cc2279d72ee744">{"x":{"filter":"none","data":[["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65","66","67","68","69","70","71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86","87","88","89","90","91","92","93","94","95","96","97","98","99","100","101","102","103","104","105","106","107","108","109","110","111","112","113","114","115","116","117","118","119","120","121","122","123","124","125","126","127","128","129","130","131","132","133","134","135","136","137","138","139","140","141","142","143","144","145","146","147","148","149","150","151","152","153","154","155","156","157","158","159","160","161","162","163","164","165","166","167","168","169","170","171","172","173","174","175","176","177","178","179","180","181","182","183","184","185","186","187","188","189","190","191","192","193","194","195","196","197","198","199","200","201","202","203","204","205","206","207","208","209","210","211","212","213","214","215","216","217","218","219","220","221","222","223","224","225","226","227","228","229","230","231","232","233","234","235","236","237","238","239","240","241","242","243","244","245","246","247","248","249","250","251","252","253","254","255","256","257","258","259","260","261","262","263","264","265","266","267","268","269","270","271","272","273","274","275","276","277","278","279","280","281","282","283","284","285","286","287","288","289","290","291","292","293","294","295","296","297","298","299","300","301","302","303","304","305","306","307","308","309","310","311","312","313","314","315","316","317","318","319","320","321","322","323","324","325","326","327","328","329","330","331","332","333","334","335","336","337","338","339","340","341","342","343","344","345","346","347","348","349","350","351","352","353","354","355","356","357","358","359","360","361","362","363","364","365","366","367","368","369","370","371","372","373","374","375","376","377","378","379","380","381","382","383","384","385","386","387","388","389","390","391","392","393","394","395","396","397","398","399","400","401","402","403","404","405","406","407","408","409","410","411","412","413","414","415","416","417","418","419","420","421","422","423","424","425","426","427","428","429","430","431","432","433","434","435","436","437","438","439","440","441","442","443","444","445","446","447","448","449","450","451","452","453","454","455"],["P04012","P03101","P03107","P06821","I6T1Z2","C5E519","C5E526","Q20MH8","Q1K9H5","B4URF7","P04012","P03101","P03107","P06821","I6T1Z2","C5E519","C5E526","Q20MH8","Q1K9H5","B4URF7","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","O56264","P04487","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q8AZK7","Q6QDQ4","P03496","P21605","Q9WPI5","Q99AU3","O56264","Q05127","D1LN35","Q5MJ03","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","O56264","P04487","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","O56264","P04487","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q8AZK7","Q6QDQ4","P03496","Q05127","P21605","Q9WPI5","Q99AU3","O56264","D1LN35","Q5MJ03","P19712-PRO_0000038050","Q8AZK7","Q6QDQ4","P03496","Q05127","P21605","Q9WPI5","Q99AU3","O56264","D1LN35","Q5MJ03","P19712-PRO_0000038050","Q9WMX2","P03495","Q997F2","Q84940","Q05127","Q9WPI5","O56264","P04487","P0C1C7","P15059","P35256","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","O56264","P04487","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","O56264","P04487","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","O56264","P04487","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","O56264","P04487","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","O56264","P04487","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q8AZK7","P03496","P21605","Q9WPI5","Q99AU3","O56264","Q05127","D1LN35","Q5MJ03","P19712-PRO_0000038050","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","O56264","P04487","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","O56264","P04487","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466","Q9WMX2","P03495","Q997F2","Q0HD54","Q84940","P23057","Q05127","Q9WPI5","P0C1C7","P15059","P0C1C6","P68466"],["IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR016024","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR000504","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458"],[1,1,0.444444444444444,0.178217821782178,0.108695652173913,0.152173913043478,0.129032258064516,0.184466019417476,0.114583333333333,0.0952380952380952,1,1,0.444444444444444,0.178217821782178,0.108695652173913,0.152173913043478,0.129032258064516,0.184466019417476,0.114583333333333,0.0952380952380952,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481481,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.113333333333333,0.482142857142857,0.141509433962264,0.214285714285714,0.229508196721311,0.458333333333333,0.25,0.166666666666667,0.452380952380952,0.470588235294118,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481481,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481481,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.113333333333333,0.482142857142857,0.141509433962264,0.166666666666667,0.214285714285714,0.229508196721311,0.458333333333333,0.25,0.452380952380952,0.470588235294118,0.161290322580645,0.113333333333333,0.482142857142857,0.141509433962264,0.166666666666667,0.214285714285714,0.229508196721311,0.458333333333333,0.25,0.452380952380952,0.470588235294118,0.161290322580645,0.180327868852459,0.333333333333333,0.323529411764706,0.578947368421053,0.433333333333333,0.180327868852459,0.231481481481481,0.2,0.379310344827586,0.666666666666667,0.260869565217391,0.358490566037736,0.392857142857143,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481481,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481481,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481481,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481481,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481481,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.113333333333333,0.141509433962264,0.214285714285714,0.229508196721311,0.458333333333333,0.25,0.166666666666667,0.452380952380952,0.470588235294118,0.161290322580645,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481481,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481481,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857],[0,0,2,26,18,20,12,12,29,34,0,0,2,26,18,20,12,12,29,34,6,4,1,7,1,0,6,2,3,9,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,8,0,18,6,2,2,3,6,2,2,6,4,1,7,1,0,6,2,3,9,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,3,9,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,8,0,18,6,6,2,2,3,2,2,5,8,0,18,6,6,2,2,3,2,2,5,6,4,1,1,6,2,3,9,1,1,0,3,0,6,4,1,7,1,0,6,2,3,9,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,3,9,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,3,9,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,3,9,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,3,9,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,8,18,6,2,2,3,6,2,2,5,6,4,1,7,1,0,6,2,3,9,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,3,9,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0,6,4,1,7,1,0,6,2,1,1,3,0],[72,0,472,126,26445,1540,5760,57,6710,47280,72,0,472,126,26445,1540,5760,57,6710,47280,1100,0,0,8390,0,0,0,1265,0,0,0,0,0,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,14263,0,0,648,42,0,0,2810,0,0,1100,0,0,8390,0,0,0,1265,0,0,0,0,0,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,1100,0,0,8390,0,0,0,1265,0,0,0,0,0,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,14263,0,0,2810,648,42,0,0,0,0,915,14263,0,0,2810,648,42,0,0,0,0,915,1100,0,0,0,0,1265,0,0,0,0,1008,0,0,1100,0,0,8390,0,0,0,1265,0,0,0,0,0,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,1100,0,0,8390,0,0,0,1265,0,0,0,0,0,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,1100,0,0,8390,0,0,0,1265,0,0,0,0,0,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,1100,0,0,8390,0,0,0,1265,0,0,0,0,0,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,1100,0,0,8390,0,0,0,1265,0,0,0,0,0,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,14263,0,648,42,0,0,2810,0,0,915,1100,0,0,8390,0,0,0,1265,0,0,0,0,0,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,1100,0,0,8390,0,0,0,1265,0,0,0,0,0,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0,3000,2740,0,8390,0,0,2810,3310,10,40,870,0],[658353,1171376,2567792,101404386,108337485,73134740,106433008,108790485,204616698,198852720,658353,1171376,2567792,101404386,108337485,73134740,106433008,108790485,204616698,198852720,40775240,73169800,24420110,41674850,14347949,8452140,95053790,40786592,148816125,218001207,21153693,56137008,62393169,20513746,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,131154558,92901249,300292170,24171633,51910208,17814093,160721415,36559150,50975708,9405872,40775240,73169800,24420110,41674850,14347949,8452140,95053790,40786592,148816125,218001207,21153693,56137008,62393169,20513746,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,40775240,73169800,24420110,41674850,14347949,8452140,95053790,40786592,148816125,218001207,21153693,56137008,62393169,20513746,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,131154558,92901249,300292170,36559150,24171633,51910208,17814093,160721415,50975708,9405872,79019340,131154558,92901249,300292170,36559150,24171633,51910208,17814093,160721415,50975708,9405872,79019340,40775240,73169800,24420110,14347949,95053790,40786592,148816125,218001207,21153693,56137008,9342072,62393169,20513746,40775240,73169800,24420110,41674850,14347949,8452140,95053790,40786592,148816125,218001207,21153693,56137008,62393169,20513746,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,40775240,73169800,24420110,41674850,14347949,8452140,95053790,40786592,148816125,218001207,21153693,56137008,62393169,20513746,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,40775240,73169800,24420110,41674850,14347949,8452140,95053790,40786592,148816125,218001207,21153693,56137008,62393169,20513746,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,40775240,73169800,24420110,41674850,14347949,8452140,95053790,40786592,148816125,218001207,21153693,56137008,62393169,20513746,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,40775240,73169800,24420110,41674850,14347949,8452140,95053790,40786592,148816125,218001207,21153693,56137008,62393169,20513746,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,131154558,300292170,24171633,51910208,17814093,160721415,36559150,50975708,9405872,79019340,40775240,73169800,24420110,41674850,14347949,8452140,95053790,40786592,148816125,218001207,21153693,56137008,62393169,20513746,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,40775240,73169800,24420110,41674850,14347949,8452140,95053790,40786592,148816125,218001207,21153693,56137008,62393169,20513746,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860,37068400,36584900,22200100,41674850,13043590,8452140,36559150,37078720,19230630,23390420,32838510,18648860],[0.000109363821536471,0,0.00018381551153676,1.24254980450254e-06,0.00024409833770832,2.10570243361773e-05,5.41185493883627e-05,5.2394287974725e-07,3.2793022590952e-05,0.000237763908886939,0.000109363821536471,0,0.00018381551153676,1.24254980450254e-06,0.00024409833770832,2.10570243361773e-05,5.41185493883627e-05,5.2394287974725e-07,3.2793022590952e-05,0.000237763908886939,2.69771557445155e-05,0,0,0.000201320460661526,0,0,0,3.10150943721898e-05,0,0,0,0,0,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,0.000108749556382173,0,0,2.68082839086627e-05,8.0908941840495e-07,0,0,7.6861743229807e-05,0,0,2.69771557445155e-05,0,0,0.000201320460661526,0,0,0,3.10150943721898e-05,0,0,0,0,0,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,2.69771557445155e-05,0,0,0.000201320460661526,0,0,0,3.10150943721898e-05,0,0,0,0,0,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,0.000108749556382173,0,0,7.6861743229807e-05,2.68082839086627e-05,8.0908941840495e-07,0,0,0,0,1.15794437159308e-05,0.000108749556382173,0,0,7.6861743229807e-05,2.68082839086627e-05,8.0908941840495e-07,0,0,0,0,1.15794437159308e-05,2.69771557445155e-05,0,0,0,0,3.10150943721898e-05,0,0,0,0,0.000107898975730438,0,0,2.69771557445155e-05,0,0,0.000201320460661526,0,0,0,3.10150943721898e-05,0,0,0,0,0,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,2.69771557445155e-05,0,0,0.000201320460661526,0,0,0,3.10150943721898e-05,0,0,0,0,0,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,2.69771557445155e-05,0,0,0.000201320460661526,0,0,0,3.10150943721898e-05,0,0,0,0,0,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,2.69771557445155e-05,0,0,0.000201320460661526,0,0,0,3.10150943721898e-05,0,0,0,0,0,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,2.69771557445155e-05,0,0,0.000201320460661526,0,0,0,3.10150943721898e-05,0,0,0,0,0,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,0.000108749556382173,0,2.68082839086627e-05,8.0908941840495e-07,0,0,7.6861743229807e-05,0,0,1.15794437159308e-05,2.69771557445155e-05,0,0,0.000201320460661526,0,0,0,3.10150943721898e-05,0,0,0,0,0,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,2.69771557445155e-05,0,0,0.000201320460661526,0,0,0,3.10150943721898e-05,0,0,0,0,0,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0,8.09314672335466e-05,7.48942869872543e-05,0,0.000201320460661526,0,0,7.6861743229807e-05,8.92695324973462e-05,5.20003764827257e-07,1.71010182801335e-06,2.64932848658481e-05,0],["O00410","O00410","O00410","O00410","O00410","O00410","O00410","O00410","O00410","O00410","O00410","O00410","O00410","O00410","O00410","O00410","O00410","O00410","O00410","O00410","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P04908","P09651","P09651","P09651","P09651","P09651","P09651","P09651","P09651","P09651","P09651","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P0C0S8","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P20671","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P22626","P62805","P62805","P62805","P62805","P62805","P62805","P62805","P62805","P62805","P62805","P62805","P62805","P62805","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q16777","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q6FI13","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q7L7L0","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q93077","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q96KK5","Q99729","Q99729","Q99729","Q99729","Q99729","Q99729","Q99729","Q99729","Q99729","Q99729","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q99878","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1","Q9BTM1"],["IPR011989|IPR034085|IPR016024","IPR011989|IPR034085|IPR016024","IPR011989|IPR034085|IPR016024","IPR011989|IPR034085|IPR016024","IPR011989|IPR034085|IPR016024","IPR011989|IPR034085|IPR016024","IPR011989|IPR034085|IPR016024","IPR011989|IPR034085|IPR016024","IPR011989|IPR034085|IPR016024","IPR011989|IPR034085|IPR016024","IPR011989|IPR016024|IPR016024","IPR011989|IPR016024|IPR016024","IPR011989|IPR016024|IPR016024","IPR011989|IPR016024|IPR016024","IPR011989|IPR016024|IPR016024","IPR011989|IPR016024|IPR016024","IPR011989|IPR016024|IPR016024","IPR011989|IPR016024|IPR016024","IPR011989|IPR016024|IPR016024","IPR011989|IPR016024|IPR016024","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR000504|IPR000504|IPR000504|IPR034803|IPR000504|IPR000504|IPR034845|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR034803|IPR000504|IPR000504|IPR034845|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR034803|IPR000504|IPR000504|IPR034845|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR034803|IPR000504|IPR000504|IPR034845|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR034803|IPR000504|IPR000504|IPR034845|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR034803|IPR000504|IPR000504|IPR034845|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR034803|IPR000504|IPR000504|IPR034845|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR034803|IPR000504|IPR000504|IPR034845|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR034803|IPR000504|IPR000504|IPR034845|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR034803|IPR000504|IPR000504|IPR034845|IPR000504|IPR000504|IPR000504","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR000504|IPR000504|IPR034486|IPR000504|IPR000504","IPR000504|IPR000504|IPR034486|IPR000504|IPR000504","IPR000504|IPR000504|IPR034486|IPR000504|IPR000504","IPR000504|IPR000504|IPR034486|IPR000504|IPR000504","IPR000504|IPR000504|IPR034486|IPR000504|IPR000504","IPR000504|IPR000504|IPR034486|IPR000504|IPR000504","IPR000504|IPR000504|IPR034486|IPR000504|IPR000504","IPR000504|IPR000504|IPR034486|IPR000504|IPR000504","IPR000504|IPR000504|IPR034486|IPR000504|IPR000504","IPR000504|IPR000504|IPR034486|IPR000504|IPR000504","IPR000504|IPR000504|IPR034486|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR000504|IPR000504","IPR035425|IPR004823|IPR009072","IPR035425|IPR004823|IPR009072","IPR035425|IPR004823|IPR009072","IPR035425|IPR004823|IPR009072","IPR035425|IPR004823|IPR009072","IPR035425|IPR004823|IPR009072","IPR035425|IPR004823|IPR009072","IPR035425|IPR004823|IPR009072","IPR035425|IPR004823|IPR009072","IPR035425|IPR004823|IPR009072","IPR035425|IPR004823|IPR009072","IPR035425|IPR004823|IPR009072","IPR035425|IPR004823|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR000504|IPR000504|IPR034846|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR034846|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR034846|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR034846|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR034846|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR034846|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR034846|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR034846|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR034846|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504","IPR000504|IPR000504|IPR034846|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504|IPR000504","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR009072|IPR007125","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR007125|IPR009072","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032454","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458","IPR032458"],["Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Domain","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site","Conserved_site"],[16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415,16415],[501,501,501,501,501,501,501,501,501,501,501,501,501,501,501,501,501,501,501,501,80,80,80,80,80,80,80,80,80,80,80,80,80,80,20,20,20,20,20,20,20,20,20,20,20,20,16,16,16,16,16,16,16,16,16,16,16,16,395,395,395,395,395,395,395,395,395,395,80,80,80,80,80,80,80,80,80,80,80,80,80,80,20,20,20,20,20,20,20,20,20,20,20,20,16,16,16,16,16,16,16,16,16,16,16,16,80,80,80,80,80,80,80,80,80,80,80,80,80,80,20,20,20,20,20,20,20,20,20,20,20,20,16,16,16,16,16,16,16,16,16,16,16,16,395,395,395,395,395,395,395,395,395,395,395,395,395,395,395,395,395,395,395,395,395,395,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,20,20,20,20,20,20,20,20,20,20,20,20,16,16,16,16,16,16,16,16,16,16,16,16,80,80,80,80,80,80,80,80,80,80,80,80,80,80,20,20,20,20,20,20,20,20,20,20,20,20,16,16,16,16,16,16,16,16,16,16,16,16,80,80,80,80,80,80,80,80,80,80,80,80,80,80,20,20,20,20,20,20,20,20,20,20,20,20,16,16,16,16,16,16,16,16,16,16,16,16,80,80,80,80,80,80,80,80,80,80,80,80,80,80,20,20,20,20,20,20,20,20,20,20,20,20,16,16,16,16,16,16,16,16,16,16,16,16,80,80,80,80,80,80,80,80,80,80,80,80,80,80,20,20,20,20,20,20,20,20,20,20,20,20,16,16,16,16,16,16,16,16,16,16,16,16,395,395,395,395,395,395,395,395,395,395,80,80,80,80,80,80,80,80,80,80,80,80,80,80,20,20,20,20,20,20,20,20,20,20,20,20,16,16,16,16,16,16,16,16,16,16,16,16,80,80,80,80,80,80,80,80,80,80,80,80,80,80,20,20,20,20,20,20,20,20,20,20,20,20,16,16,16,16,16,16,16,16,16,16,16,16],[0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.0305208650624429,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.024063356685958,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00487359122753579,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.00121839780688395,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158,0.000974718245507158],[9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606,9606],[10580,333760,333760,211044,387139,643960,643960,381518,381518,381518,10580,333760,333760,211044,387139,643960,643960,381518,381518,381518,333284,381517,121791,383543,10915,11214,128952,380964,130763,10299,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,10377,28344,211044,10254,380964,88776,130763,128952,680716,10376,333284,381517,121791,383543,10915,11214,128952,380964,130763,10299,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,130763,10299,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,10377,28344,211044,128952,10254,380964,88776,130763,680716,10376,11097,10377,28344,211044,128952,10254,380964,88776,130763,680716,10376,11097,333284,381517,121791,10915,128952,380964,130763,10299,121791,10254,33727,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,130763,10299,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,130763,10299,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,130763,10299,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,130763,10299,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,130763,10299,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,10377,211044,10254,380964,88776,130763,128952,680716,10376,11097,333284,381517,121791,383543,10915,11214,128952,380964,130763,10299,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,130763,10299,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254,333284,381517,121791,383543,10915,11214,128952,380964,121791,10254,928303,10254],[3,4,9,101,138,92,124,103,192,210,3,4,9,101,138,92,124,103,192,210,61,60,34,70,19,12,60,61,108,145,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,150,56,212,42,61,24,108,60,42,17,61,60,34,70,19,12,60,61,108,145,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,108,145,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,150,56,212,60,42,61,24,108,42,17,93,150,56,212,60,42,61,24,108,42,17,93,61,60,34,19,60,61,108,145,29,36,23,53,28,61,60,34,70,19,12,60,61,108,145,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,108,145,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,108,145,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,108,145,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,108,145,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,150,212,42,61,24,108,60,42,17,93,61,60,34,70,19,12,60,61,108,145,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,108,145,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28,61,60,34,70,19,12,60,61,29,36,53,28],[23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,19,19,19,19,19,19,19,19,19,19,19,19,19,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,11,11,11,11,11,11,11,11,11,11,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14],[6,7,18,129,218,144,219,204,294,313,6,7,18,129,218,144,219,204,294,313,89,69,42,92,20,10,59,74,104,186,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,254,52,249,49,74,39,104,59,54,23,89,69,42,92,20,10,59,74,104,186,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,104,186,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,254,52,249,59,49,74,39,104,54,23,150,254,52,249,59,49,74,39,104,54,23,150,89,69,42,20,59,74,104,186,30,30,31,47,51,89,69,42,92,20,10,59,74,104,186,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,104,186,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,104,186,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,104,186,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,104,186,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,254,249,49,74,39,104,59,54,23,150,89,69,42,92,20,10,59,74,104,186,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,104,186,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51,89,69,42,92,20,10,59,74,30,30,47,51],[220,220,220,220,220,220,220,220,220,220,220,220,220,220,220,220,220,220,220,220,37,37,37,37,37,37,37,37,37,37,37,37,37,37,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,154,154,154,154,154,154,154,154,154,154,37,37,37,37,37,37,37,37,37,37,37,37,37,37,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,37,37,37,37,37,37,37,37,37,37,37,37,37,37,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,154,154,154,154,154,154,154,154,154,154,154,154,154,154,154,154,154,154,154,154,154,154,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,37,37,37,37,37,37,37,37,37,37,37,37,37,37,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,37,37,37,37,37,37,37,37,37,37,37,37,37,37,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,37,37,37,37,37,37,37,37,37,37,37,37,37,37,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,37,37,37,37,37,37,37,37,37,37,37,37,37,37,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,154,154,154,154,154,154,154,154,154,154,37,37,37,37,37,37,37,37,37,37,37,37,37,37,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,37,37,37,37,37,37,37,37,37,37,37,37,37,37,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18,18],[3,4,4,18,15,14,16,19,22,20,3,4,4,18,15,14,16,19,22,20,11,20,11,10,11,10,26,11,25,29,11,24,19,11,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,17,27,30,9,14,11,27,10,19,8,11,20,11,10,11,10,26,11,25,29,11,24,19,11,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,20,11,10,11,10,26,11,25,29,11,24,19,11,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,17,27,30,10,9,14,11,27,19,8,15,17,27,30,10,9,14,11,27,19,8,15,11,20,11,11,26,11,25,29,11,24,6,19,11,11,20,11,10,11,10,26,11,25,29,11,24,19,11,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,20,11,10,11,10,26,11,25,29,11,24,19,11,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,20,11,10,11,10,26,11,25,29,11,24,19,11,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,20,11,10,11,10,26,11,25,29,11,24,19,11,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,20,11,10,11,10,26,11,25,29,11,24,19,11,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,17,30,9,14,11,27,10,19,8,15,11,20,11,10,11,10,26,11,25,29,11,24,19,11,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,11,20,11,10,11,10,26,11,25,29,11,24,19,11,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10],[1,1,0.444444444444444,0.178217821782178,0.108695652173913,0.152173913043478,0.129032258064516,0.184466019417476,0.114583333333333,0.0952380952380952,1,1,0.444444444444444,0.178217821782178,0.108695652173913,0.152173913043478,0.129032258064516,0.184466019417476,0.114583333333333,0.0952380952380952,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481482,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.113333333333333,0.482142857142857,0.141509433962264,0.214285714285714,0.229508196721311,0.458333333333333,0.25,0.166666666666667,0.452380952380952,0.470588235294118,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481482,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481482,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.113333333333333,0.482142857142857,0.141509433962264,0.166666666666667,0.214285714285714,0.229508196721311,0.458333333333333,0.25,0.452380952380952,0.470588235294118,0.161290322580645,0.113333333333333,0.482142857142857,0.141509433962264,0.166666666666667,0.214285714285714,0.229508196721311,0.458333333333333,0.25,0.452380952380952,0.470588235294118,0.161290322580645,0.180327868852459,0.333333333333333,0.323529411764706,0.578947368421053,0.433333333333333,0.180327868852459,0.231481481481482,0.2,0.379310344827586,0.666666666666667,0.260869565217391,0.358490566037736,0.392857142857143,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481482,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481482,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481482,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481482,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481482,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.113333333333333,0.141509433962264,0.214285714285714,0.229508196721311,0.458333333333333,0.25,0.166666666666667,0.452380952380952,0.470588235294118,0.161290322580645,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481482,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.180327868852459,0.333333333333333,0.323529411764706,0.142857142857143,0.578947368421053,0.833333333333333,0.433333333333333,0.180327868852459,0.231481481481482,0.2,0.379310344827586,0.666666666666667,0.358490566037736,0.392857142857143,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857,0.163934426229508,0.166666666666667,0.294117647058824,0.142857142857143,0.526315789473684,0.833333333333333,0.166666666666667,0.163934426229508,0.344827586206897,0.277777777777778,0.188679245283019,0.357142857142857],[32.7644710578842,32.7644710578842,14.5619871368374,5.83921266378135,3.56135554977003,4.98589776967804,4.22767368488829,6.04393155436699,3.7542623087159,3.12042581503659,32.7644710578842,32.7644710578842,14.5619871368374,5.83921266378135,3.56135554977003,4.98589776967804,4.22767368488829,6.04393155436699,3.7542623087159,3.12042581503659,37.0010245901639,68.3958333333333,66.3841911764706,29.3125,118.792763157895,170.989583333333,88.9145833333333,37.0010245901639,47.4971064814815,41.0375,77.8297413793103,136.791666666667,73.5577830188679,80.609375,134.549180327869,136.791666666667,241.397058823529,117.25,431.973684210526,683.958333333333,136.791666666667,134.549180327869,283.01724137931,227.986111111111,154.858490566038,293.125,168.186475409836,170.989583333333,301.746323529412,146.5625,539.967105263158,854.947916666667,170.989583333333,168.186475409836,353.771551724138,284.982638888889,193.573113207547,366.40625,4.70978902953587,20.0363924050633,5.88070217339384,8.90506329113924,9.5376634156464,19.04694092827,10.3892405063291,6.92616033755274,18.7995780590717,19.5562174236783,37.0010245901639,68.3958333333333,66.3841911764706,29.3125,118.792763157895,170.989583333333,88.9145833333333,37.0010245901639,47.4971064814815,41.0375,77.8297413793103,136.791666666667,73.5577830188679,80.609375,134.549180327869,136.791666666667,241.397058823529,117.25,431.973684210526,683.958333333333,136.791666666667,134.549180327869,283.01724137931,227.986111111111,154.858490566038,293.125,168.186475409836,170.989583333333,301.746323529412,146.5625,539.967105263158,854.947916666667,170.989583333333,168.186475409836,353.771551724138,284.982638888889,193.573113207547,366.40625,37.0010245901639,68.3958333333333,66.3841911764706,29.3125,118.792763157895,170.989583333333,88.9145833333333,37.0010245901639,47.4971064814815,41.0375,77.8297413793103,136.791666666667,73.5577830188679,80.609375,134.549180327869,136.791666666667,241.397058823529,117.25,431.973684210526,683.958333333333,136.791666666667,134.549180327869,283.01724137931,227.986111111111,154.858490566038,293.125,168.186475409836,170.989583333333,301.746323529412,146.5625,539.967105263158,854.947916666667,170.989583333333,168.186475409836,353.771551724138,284.982638888889,193.573113207547,366.40625,4.70978902953587,20.0363924050633,5.88070217339384,6.92616033755274,8.90506329113924,9.5376634156464,19.04694092827,10.3892405063291,18.7995780590717,19.5562174236783,6.70273581053491,4.70978902953587,20.0363924050633,5.88070217339384,6.92616033755274,8.90506329113924,9.5376634156464,19.04694092827,10.3892405063291,18.7995780590717,19.5562174236783,6.70273581053491,37.0010245901639,68.3958333333333,66.3841911764706,118.792763157895,88.9145833333333,37.0010245901639,47.4971064814815,41.0375,77.8297413793103,136.791666666667,53.5271739130435,73.5577830188679,80.609375,37.0010245901639,68.3958333333333,66.3841911764706,29.3125,118.792763157895,170.989583333333,88.9145833333333,37.0010245901639,47.4971064814815,41.0375,77.8297413793103,136.791666666667,73.5577830188679,80.609375,134.549180327869,136.791666666667,241.397058823529,117.25,431.973684210526,683.958333333333,136.791666666667,134.549180327869,283.01724137931,227.986111111111,154.858490566038,293.125,168.186475409836,170.989583333333,301.746323529412,146.5625,539.967105263158,854.947916666667,170.989583333333,168.186475409836,353.771551724138,284.982638888889,193.573113207547,366.40625,37.0010245901639,68.3958333333333,66.3841911764706,29.3125,118.792763157895,170.989583333333,88.9145833333333,37.0010245901639,47.4971064814815,41.0375,77.8297413793103,136.791666666667,73.5577830188679,80.609375,134.549180327869,136.791666666667,241.397058823529,117.25,431.973684210526,683.958333333333,136.791666666667,134.549180327869,283.01724137931,227.986111111111,154.858490566038,293.125,168.186475409836,170.989583333333,301.746323529412,146.5625,539.967105263158,854.947916666667,170.989583333333,168.186475409836,353.771551724138,284.982638888889,193.573113207547,366.40625,37.0010245901639,68.3958333333333,66.3841911764706,29.3125,118.792763157895,170.989583333333,88.9145833333333,37.0010245901639,47.4971064814815,41.0375,77.8297413793103,136.791666666667,73.5577830188679,80.609375,134.549180327869,136.791666666667,241.397058823529,117.25,431.973684210526,683.958333333333,136.791666666667,134.549180327869,283.01724137931,227.986111111111,154.858490566038,293.125,168.186475409836,170.989583333333,301.746323529412,146.5625,539.967105263158,854.947916666667,170.989583333333,168.186475409836,353.771551724138,284.982638888889,193.573113207547,366.40625,37.0010245901639,68.3958333333333,66.3841911764706,29.3125,118.792763157895,170.989583333333,88.9145833333333,37.0010245901639,47.4971064814815,41.0375,77.8297413793103,136.791666666667,73.5577830188679,80.609375,134.549180327869,136.791666666667,241.397058823529,117.25,431.973684210526,683.958333333333,136.791666666667,134.549180327869,283.01724137931,227.986111111111,154.858490566038,293.125,168.186475409836,170.989583333333,301.746323529412,146.5625,539.967105263158,854.947916666667,170.989583333333,168.186475409836,353.771551724138,284.982638888889,193.573113207547,366.40625,37.0010245901639,68.3958333333333,66.3841911764706,29.3125,118.792763157895,170.989583333333,88.9145833333333,37.0010245901639,47.4971064814815,41.0375,77.8297413793103,136.791666666667,73.5577830188679,80.609375,134.549180327869,136.791666666667,241.397058823529,117.25,431.973684210526,683.958333333333,136.791666666667,134.549180327869,283.01724137931,227.986111111111,154.858490566038,293.125,168.186475409836,170.989583333333,301.746323529412,146.5625,539.967105263158,854.947916666667,170.989583333333,168.186475409836,353.771551724138,284.982638888889,193.573113207547,366.40625,4.70978902953587,5.88070217339384,8.90506329113924,9.5376634156464,19.04694092827,10.3892405063291,6.92616033755274,18.7995780590717,19.5562174236783,6.70273581053491,37.0010245901639,68.3958333333333,66.3841911764706,29.3125,118.792763157895,170.989583333333,88.9145833333333,37.0010245901639,47.4971064814815,41.0375,77.8297413793103,136.791666666667,73.5577830188679,80.609375,134.549180327869,136.791666666667,241.397058823529,117.25,431.973684210526,683.958333333333,136.791666666667,134.549180327869,283.01724137931,227.986111111111,154.858490566038,293.125,168.186475409836,170.989583333333,301.746323529412,146.5625,539.967105263158,854.947916666667,170.989583333333,168.186475409836,353.771551724138,284.982638888889,193.573113207547,366.40625,37.0010245901639,68.3958333333333,66.3841911764706,29.3125,118.792763157895,170.989583333333,88.9145833333333,37.0010245901639,47.4971064814815,41.0375,77.8297413793103,136.791666666667,73.5577830188679,80.609375,134.549180327869,136.791666666667,241.397058823529,117.25,431.973684210526,683.958333333333,136.791666666667,134.549180327869,283.01724137931,227.986111111111,154.858490566038,293.125,168.186475409836,170.989583333333,301.746323529412,146.5625,539.967105263158,854.947916666667,170.989583333333,168.186475409836,353.771551724138,284.982638888889,193.573113207547,366.40625],[10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,14,14,14,14,14,14,14,14,14,14,14,14,14,14,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,10,10,10,10,10,10,10,10,10,10,14,14,14,14,14,14,14,14,14,14,14,14,14,14,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,14,14,14,14,14,14,14,14,14,14,14,14,14,14,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,13,13,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,14,14,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,14,14,14,14,14,14,14,14,14,14,14,14,14,14,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,14,14,14,14,14,14,14,14,14,14,14,14,14,14,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,14,14,14,14,14,14,14,14,14,14,14,14,14,14,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,14,14,14,14,14,14,14,14,14,14,14,14,14,14,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,10,10,10,10,10,10,10,10,10,10,14,14,14,14,14,14,14,14,14,14,14,14,14,14,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,14,14,14,14,14,14,14,14,14,14,14,14,14,14,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12]],"container":"<table class=\"display\">\n  <thead>\n    <tr>\n      <th> <\/th>\n      <th>IDs_interactor_viral<\/th>\n      <th>IDs_domain_human<\/th>\n      <th>observed_statistic<\/th>\n      <th>YmissingZ_perX<\/th>\n      <th>higher_counts<\/th>\n      <th>not_missing<\/th>\n      <th>p.value<\/th>\n      <th>IDs_interactor_human<\/th>\n      <th>all_IDs_domain<\/th>\n      <th>domain_type<\/th>\n      <th>N_prot_w_interactors<\/th>\n      <th>domain_count<\/th>\n      <th>domain_frequency<\/th>\n      <th>Taxid_interactor_human<\/th>\n      <th>Taxid_interactor_viral<\/th>\n      <th>IDs_interactor_viral_degree<\/th>\n      <th>IDs_interactor_human_degree<\/th>\n      <th>IDs_domain_human_per_IDs_interactor_viral<\/th>\n      <th>IDs_interactor_viral_per_IDs_domain_human<\/th>\n      <th>domain_count_per_IDs_interactor_viral<\/th>\n      <th>domain_frequency_per_IDs_interactor_viral<\/th>\n      <th>fold_enrichment<\/th>\n      <th>N_viral_per_human_w_domain<\/th>\n    <\/tr>\n  <\/thead>\n<\/table>","options":{"columnDefs":[{"className":"dt-right","targets":[3,4,5,6,7,11,12,13,14,15,16,17,18,19,20,21,22,23]},{"orderable":false,"targets":0}],"order":[],"autoWidth":false,"orderClasses":false}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

## R session information


```r
save(list = ls(), file="./processed_data_files/what_we_find_VS_ELM_clust11092017.RData")
#R.utils::gzip(filename = "./processed_data_files/what_we_find_VS_ELM_clust11092017.RData",
#              destname = "./processed_data_files/what_we_find_VS_ELM_clust11092017.RData.gz",
#              remove = T, overwrite = T)
Sys.Date()
```

```
## [1] "2017-09-27"
```

```r
sessionInfo()
```

```
## R version 3.4.1 (2017-06-30)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Red Hat Enterprise Linux Server 7.3 (Maipo)
## 
## Matrix products: default
## BLAS: /hps/nobackup/research/petsalaki/users/vitalii/R/lib64/R/lib/libRblas.so
## LAPACK: /hps/nobackup/research/petsalaki/users/vitalii/R/lib64/R/lib/libRlapack.so
## 
## locale:
##  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
##  [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
##  [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    parallel  methods   stats     graphics  grDevices utils    
## [8] datasets  base     
## 
## other attached packages:
##  [1] RColorBrewer_1.1-2   GGally_1.3.2         ggplot2_2.2.1       
##  [4] rtracklayer_1.36.3   GenomicRanges_1.28.3 GenomeInfoDb_1.12.2 
##  [7] MItools_0.1.23       Biostrings_2.44.1    XVector_0.16.0      
## [10] data.table_1.10.4    PSICQUIC_1.14.0      plyr_1.8.4          
## [13] httr_1.3.1           biomaRt_2.32.1       IRanges_2.10.2      
## [16] S4Vectors_0.14.3     BiocGenerics_0.22.0  rmarkdown_1.6       
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.12               lattice_0.20-35           
##  [3] Rsamtools_1.28.0           rprojroot_1.2             
##  [5] digest_0.6.12              R6_2.2.2                  
##  [7] backports_1.1.0            RSQLite_2.0               
##  [9] evaluate_0.10.1            zlibbioc_1.22.0           
## [11] rlang_0.1.2                lazyeval_0.2.0            
## [13] blob_1.1.0                 R.utils_2.5.0             
## [15] R.oo_1.21.0                Matrix_1.2-11             
## [17] DT_0.2                     qvalue_2.8.0              
## [19] gsubfn_0.6-6               proto_1.0.0               
## [21] splines_3.4.1              BiocParallel_1.10.1       
## [23] downloader_0.4             stringr_1.2.0             
## [25] htmlwidgets_0.9            RCurl_1.95-4.8            
## [27] bit_1.1-12                 munsell_0.4.3             
## [29] DelayedArray_0.2.7         compiler_3.4.1            
## [31] htmltools_0.3.6            SummarizedExperiment_1.6.3
## [33] tibble_1.3.4               GenomeInfoDbData_0.99.0   
## [35] matrixStats_0.52.2         XML_3.98-1.9              
## [37] reshape_0.8.7              GenomicAlignments_1.12.1  
## [39] bitops_1.0-6               R.methodsS3_1.7.1         
## [41] grid_3.4.1                 jsonlite_1.5              
## [43] gtable_0.2.0               DBI_0.7                   
## [45] magrittr_1.5               scales_0.5.0              
## [47] stringi_1.1.5              reshape2_1.4.2            
## [49] tools_3.4.1                bit64_0.9-7               
## [51] Biobase_2.36.2             yaml_2.1.14               
## [53] AnnotationDbi_1.38.1       colorspace_1.3-2          
## [55] memoise_1.1.0              knitr_1.17
```
