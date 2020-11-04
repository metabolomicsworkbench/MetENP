# MetENP
Metabolite enrichment analysis and their associated enriched pathways.

## Introduction to Enrichment Pipeline

MetENP is a R package that enables detection of significant metabolites from metabolite information 
(names or names and concentration along with metadata information) and provides

1. Enrichment score of metabolite class,
2. Maps to pathway of the species of choice,
3. Calculate enrichment score of pathways,
4. Plots the pathways and shows the metabolite increase or decrease
5. Gets gene info, reaction info, enzyme info

For more info, check out the vignette.


## Installation

install.packages("devtools")
library("devtools")

MetENP package depends on following Bioconductor packages to function properly: KEGGREST, KEGGgraph, pathview and KEGG.db. 
 You may need to install these via:


if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")    
 BiocManager::install("KEGGREST")<br/>
 BiocManager::install("KEGGgraph")<br/>
 BiocManager::install("pathview")<br/>
 BiocManager::install("KEGG.db")<br/>
 
 #### Now proceed with installation
 
 devtools::install("MetENP")
 
 If above steps gives error:
Install other dependencies and then try installing again: plyr,dplyr,tidyr,purrr,tidygraph,reshape2,ggplot2,ggrepel,
    igraph,ggraph,httr,stringr,jsonlite,rjson,tidyverse,magrittr
    
    
#### If you do not wish to install, alternatively, download from github(https://github.com/metabolomicsworkbench/MetENP) and load libraries and functions

suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(tidygraph))
suppressMessages(library(KEGGREST))
suppressMessages(library(KEGGgraph))
suppressMessages(library(KEGG.db))
suppressMessages(library(pathview))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(igraph))
suppressMessages(library(ggraph))
suppressMessages(library(httr))
suppressMessages(library(stringr))
suppressMessages(library(jsonlite))
suppressMessages(library(rjson))
suppressMessages(library(tidyverse))

### And load all the function with appropriate path (replace 'path' to your own path). 
### Please note this step is needed only when you do not wish to download or are hving difficulty in downloading the package

source('path/compoundinfo.R')
source('path/anova_ana.R')
source('path/met_pathways.R')
source('path/mapspspath.R')
source('path/metclassenrichment.R')
source('path/metcountplot.R')
source('path/getmwstudies.R')
source('path/path_enrichmentscore.R')
source('path/pathinfo.R')
source('path/plot_met_enrichment.R')
source('path/plot_volcano.R')
source('path/rxninfo.R')
source('path/significant_met.R')
source('path/significant_met_own.R')
source('path/enzyme_gene_info.R')
source('path/plot_heatmap.R')
source('path/plot_pathway_networks.R')
source('path/react_substrate.R')
source('path/dotplot_met_class_path.R')
source('path/convert_refmet.R')
source('path/map_keggid.R')
source('path/partial_join.R')
source('path/getExtension.R')
source('path/separate_data.R')

Now please follow example in the vignette
