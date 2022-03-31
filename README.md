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
Contact: biosonal@gmail.com; kschoudhary@eng.ucsd.edu


## Installation

install.packages("devtools")<br/>
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
1) <strong>Through devtools </strong></br>

 devtools::install("MetENP")
 
 If above steps gives error:
Install other dependencies and then try installing again: plyr,dplyr,tidyr,purrr,tidygraph,reshape2,ggplot2,ggrepel,
    igraph,ggraph,httr,stringr,jsonlite,rjson,tidyverse,magrittr

2) <strong>Through Anaconda </strong></br>

   A conda environment.yml file is present in the repository. </br>
   This file can be used be with Anaconda to install all of the R requirements. </br>

   Create the conda environment and activate the environment by running:

```
   conda env create -n metenp -f environment.yml
   conda activate metenp
```

#### If you do not wish to install, alternatively, download from github(https://github.com/metabolomicsworkbench/MetENP) and load libraries and functions

suppressMessages(library(plyr))<br/>
suppressMessages(library(dplyr))<br/>
suppressMessages(library(tidyr))<br/>
suppressMessages(library(tidygraph))<br/>
suppressMessages(library(KEGGREST))<br/>
suppressMessages(library(KEGGgraph))<br/>
suppressMessages(library(KEGG.db))<br/>
suppressMessages(library(pathview))<br/>
suppressMessages(library(reshape2))<br/>
suppressMessages(library(ggplot2))<br/>
suppressMessages(library(ggrepel))<br/>
suppressMessages(library(igraph))<br/>
suppressMessages(library(ggraph))<br/>
suppressMessages(library(httr))<br/>
suppressMessages(library(stringr))<br/>
suppressMessages(library(jsonlite))<br/>
suppressMessages(library(rjson))<br/>
suppressMessages(library(tidyverse))<br/>

#### And load all the function with appropriate path (replace 'path' to your own path). 
#### Please note this step is needed only when you do not wish to download or are hving difficulty in downloading the package

source('path/compoundinfo.R')<br/>
source('path/anova_ana.R')<br/>
source('path/met_pathways.R')<br/>
source('path/mapspspath.R')<br/>
source('path/metclassenrichment.R')<br/>
source('path/metcountplot.R')<br/>
source('path/getmwstudies.R')<br/>
source('path/path_enrichmentscore.R')<br/>
source('path/pathinfo.R')<br/>
source('path/plot_met_enrichment.R')<br/>
source('path/plot_volcano.R')<br/>
source('path/rxninfo.R')<br/>
source('path/significant_met.R')<br/>
source('path/significant_met_own.R')<br/>
source('path/enzyme_gene_info.R')<br/>
source('path/plot_heatmap.R')<br/>
source('path/plot_pathway_networks.R')<br/>
source('path/react_substrate.R')<br/>
source('path/dotplot_met_class_path.R')<br/>
source('path/convert_refmet.R')<br/>
source('path/map_keggid.R')<br/>
source('path/partial_join.R')<br/>
source('path/getExtension.R')<br/>
source('path/separate_data.R')<br/>

Now please follow example in the vignette

Run the vignette Jupyter Notebook on the web using My Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/metabolomicsworkbench/MetENP/test01?filepath=vignettes%2FMetENP_vignette_Jupyter_notebook.ipynb)
