---
title: 'Degree distribution in the human-human, human-viral and viral-viral network'
author: "Vitalii Kleshchevnikov"
date: "04/09/2017"
output: 
  html_document: 
    keep_md: yes
    toc: yes
---



### Load PPI data and show degree distribution


```r
# Load data.
HumanViralPPI = loadHumanViralPPI(directory = "./data_files/",
                                  loadIntActFTP_dir = "./data_files/IntActRelease_2017Nov13/",
                                  release = "2017Nov13")
```

```
## ... loading local copy ...
```

```
## Read 0.0% of 790620 rowsRead 6.3% of 790620 rowsRead 16.4% of 790620 rowsRead 22.8% of 790620 rowsRead 31.6% of 790620 rowsRead 36.7% of 790620 rowsRead 44.3% of 790620 rowsRead 51.9% of 790620 rowsRead 59.4% of 790620 rowsRead 69.6% of 790620 rowsRead 79.7% of 790620 rowsRead 91.1% of 790620 rowsRead 790620 rows and 42 (of 42) columns from 2.995 GB file in 00:00:28
```

```r
plotHumanViralDegree = function(degree_distributions, y_text = 2, median_text_size = 5.5){
    ggplot(degree_distributions, aes(x = N, fill = legend, color = legend)) + geom_density() + facet_grid(legend ~ data_name) + theme_light() +
        theme(strip.text.y = element_text(size = 18,angle = 0),
              strip.text.x = element_text(size = 18,angle = 0),
              panel.grid.major =  element_line(color = 'grey', size = 0.2, linetype = 'solid'),
              panel.background = element_rect(fill = '#FFFFFF', colour = 'grey'),
              legend.position = "none",
              axis.title = element_text(size = 18),
              axis.text.x = element_text(size = 14)) +
        scale_x_log10(breaks = c(1, 10, 100, 1000),
                      labels = c(1, 10, 100, 1000)) +
        xlab("N interacting partners per protein")+
        geom_vline(aes(xintercept = medianN, color = legend)) + 
        geom_text(aes(x = medianN*9.5,label = medianN_lab, y = y_text), color = "black", size = median_text_size)+
        annotation_logticks(sides = "b", color = "#444444")
}

degree_distributions_all = humanViralDegree(data = HumanViralPPI,
  directory = "./data_files/", data_name = "Full IntAct")
```

```
## loading local copy of MI ontology
```

```
## loading local copy of MI ontology
## loading local copy of MI ontology
```

```r
plotHumanViralDegree(degree_distributions_all)
```

![](Degree_distribution_in_the_network_files/figure-html/load-1.png)<!-- -->

#### Degree distribution based on all data, two-hybrid or affinity purification - mass spectrometry


```r
degree_distributions_two = humanViralDegree(data = HumanViralPPI,
  directory = "./data_files/", Interaction_detection_methods = "MI:0018",
  data_name = "two-hybrid data")
```

```
## loading local copy of MI ontology
## loading local copy of MI ontology
## loading local copy of MI ontology
```

```r
plotHumanViralDegree(degree_distributions_two)
```

![](Degree_distribution_in_the_network_files/figure-html/degree_all_data-1.png)<!-- -->

```r
degree_distributionsAPMS = humanViralDegree(data = HumanViralPPI,
  directory = "./data_files/", Interaction_detection_methods = "MI:0004",
  Identification_method = "MI:0433", PMIDs = NULL, inverse_filter = F,
  data_name = "AP-MS data")
```

```
## loading local copy of MI ontology
## loading local copy of MI ontology
## loading local copy of MI ontology
```

```r
plotHumanViralDegree(degree_distributionsAPMS)
```

![](Degree_distribution_in_the_network_files/figure-html/degree_all_data-2.png)<!-- -->

```r
degree_distributions = rbind(degree_distributions_all,degree_distributions_two,degree_distributionsAPMS)
degree_distributions[, legend := gsub(":.+$","",legend)]
p1 = plotHumanViralDegree(degree_distributions, y_text = 3)
p1
```

![](Degree_distribution_in_the_network_files/figure-html/degree_all_data-3.png)<!-- -->

#### Degree distribution within human network based on specific datasets (Vidal's and Mann's)


```r
degree_distributionsVidal = humanViralDegree(data = HumanViralPPI,
  directory = "./data_files/", PMIDs = c("25416956", "unassigned1304"),
  data_name = "Vidal's data \nincluding unpublished")
```

```
## loading local copy of MI ontology
## loading local copy of MI ontology
## loading local copy of MI ontology
```

```r
plotHumanViralDegree(degree_distributionsVidal)
```

![](Degree_distribution_in_the_network_files/figure-html/degree_Vidal_data-1.png)<!-- -->

```r
degree_distributionsMann = humanViralDegree(data = HumanViralPPI,
  directory = "./data_files/", PMIDs = "26496610",
  data_name = "Mann's \n2015 paper")
```

```
## loading local copy of MI ontology
## loading local copy of MI ontology
## loading local copy of MI ontology
```

```r
plotHumanViralDegree(degree_distributionsMann)
```

![](Degree_distribution_in_the_network_files/figure-html/degree_Vidal_data-2.png)<!-- -->


```r
degree_distributions = rbind(degree_distributions,degree_distributionsVidal,degree_distributionsMann)
degree_distributions[, legend := gsub(":.+$","",legend)]
degree_distributions = degree_distributions[legend %in% c("full human-human", "viral-interacting proteins, human-human"),]
degree_distributions[legend == "viral-interacting proteins, human-human",
                     legend := "viral-interacting proteins,\nhuman-human"]
p4 = plotHumanViralDegree(degree_distributions, y_text = 1, median_text_size = 7)
p4
```

![](Degree_distribution_in_the_network_files/figure-html/degree_Vidal_data_comb-1.png)<!-- -->
