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
Contact: biosonal@gmail.com; kschoudhary@eng.ucsd.edu, Mano Maurya (mano@sdsc.edu) or Sumana Srinivasan (susrinivasan@eng.ucsd.edu).

A BioRxiv preprint on MetENP is available at: https://www.biorxiv.org/content/10.1101/2020.11.20.391912v1 (DOI: https://doi.org/10.1101/2020.11.20.391912). Please follow the details at the GitHub page (this page) to ensure you use the latest version.

## Installation

Download the MetENP R package: On linux/unix command prompt (or download a zip file by clicking on "Code" and unpack the zip file to a local folder):

$git clone https://github.com/metabolomicsworkbench/MetENP.git MetENP

Start R and install devtools.

install.packages("devtools")<br/>
library("devtools")

MetENP package depends on following Bioconductor packages to function properly: KEGGREST, KEGGgraph, and pathview. 
 You may need to install these via:


if (!requireNamespace("BiocManager", quietly = TRUE))
 install.packages("BiocManager")    
 BiocManager::install("KEGGREST")<br/>
 BiocManager::install("KEGGgraph")<br/>
 BiocManager::install("pathview")<br/>
 
 #### Now proceed with installation
1) <strong>Through devtools </strong></br>

 devtools::install("MetENP")
 
 If you do not have admin priviligages, you can install the packages in the user area, e.g., /home/username/.local/R. Please see detailed instructions below.
 
 If above steps give error:
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
#suppressMessages(library(KEGG.db)); # only for BiocManager version < 3.13<br/>
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
suppressMessages(library(tictoc))<br/>

#### And load all the function with appropriate path (replace 'path' to your own path, e.g., they are inside R folder relative to the MetENP folder).

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

Run the vignette Jupyter Notebook on the web using My Binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/metabolomicsworkbench/MetENP/HEAD?labpath=vignettes%2FMetENP_vignette_Jupyter_notebook.ipynb)

## Detailed instructions 

The syntax of paths is for linux/unix operating system. It can be adjusted for Windows.

## (1) If one needs to install Jupyter</br>

If already installed, this section can be ignored.</br>

[[ for general information only:</br>
basic commands: whereis python3 # can check version by starting python3</br>
system python packages go to: /usr/lib/python3.9 /usr/lib64/python3.9 [if installed via sudo dnf] /usr/local/lib/python3.9  /usr/local/lib64/python3.9 [if installed via sudo pip3]</br>
user/local python typically goes into $HOME/.local/various-folders</br>
system R is at /usr/bin/R /usr/lib64/R</br>
system R packages go into /usr/lib64/R/library</br>
]]
#### # is comment, $ at start indicates linux command
Be in your home folder, /home/username: </br>
&#35; install jupyter in user area, in the folder /home/username/.local/bunch-of-folders</br>
$pip3 install --user --no-cache-dir jupyter</br>
&#35; test it</br>
$ jupyter notebook --ip=your_ip_address_format_123.249.124.012 --port=8080</br>
Go to the page listed, e.g.,</br>
http://123.249.124.012:8080/?token=4228fsdrjh346t3fdgve716452997a25f3e36b0dc2c3f02a3a0aa34</br>
User can try to open any existing jupyter notebook if they are in the folder/subfolder</br>
ctrl-C to stop</br>
Install R package IRkernel # need to install R kernel for jupyter in user area; set libloc to user area, e.g., ${HOME}/.local/R if you do not have admin privileges.
&#35; start R, being in home area ; > indicates R prompt</br>
$R</br>
&#62;reposlink = 'http://cran.r-project.org'; libloc = "/usr/lib64/R/library/";</br>
&#62;pkgnames = c('IRkernel'); install.packages(pkgnames, repos=reposlink, lib=libloc);</br>
&#62;IRkernel::installspec() # for only current user</br>
&#62;q()</br>
$ls -al .local/share/jupyter/kernels/</br>
The output will appear similar to the lines below:</br>
total 0</br>
drwxrwxr-x. 4 username username  31 Nov  1 00:08 .</br>
drwxrwxr-x. 7 username username 140 Nov  2 10:45 ..</br>
drwxr-xr-x. 2 username username  64 Nov  1 00:08 ir</br>
drwxrwxr-x. 2 username username  69 Oct 27 22:40 python3</br>
 
#### # now jupyter notebooks based on R code should work after you select R kernel after starting jupyter
$ jupyter notebook --ip=123.249.124.012 --port=8080</br>

## (2)	How to Install MetENP R package through R devtools in user area</br>

If already installed, this section can be ignored.</br>
Download or clone MetENP folder from github. Be careful not to overwrite existing folders [create and be in a different folder as needed].</br>
$git clone https://github.com/metabolomicsworkbench/MetENP.git MetENP</br>
The necessary files for installing MetENP R package are located inside the sub-folder R inside the MetENP folder.</br>
#### # install MetENP R package in user area (e.g., /home/username/.local/R, or /home/username/some-folder/R):</br>

$cd ~/.local; mkdir R;</br>

#### Go to the folder one-level up of MetENP folder.
If you want to install MetENP R package in the folder USER_PWD/R below, create the R folder.</br>
[username@server one-level-up-of-MetENP]$mkdir R</br>
[username@server one-level-up-of-MetENP]$R</br></br>

&#35; If devtools is not already installed for all, install it in system R or user R area (see how to set libloc below)</br>
&#62;USER_HOME=Sys.getenv("HOME"); # so that we don’t need to hard code /home/username</br>
&#62;reposlink = 'http://cran.r-project.org'; libloc = paste0(USER_HOME, “/.local/R/");</br>
&#62;#pkgnames = c("devtools"); install.packages(pkgnames, repos=reposlink, lib=libloc);</br>
&#62;library("devtools");</br>
</br>
&#62;devtools::install("MetENP", args = paste0("--library=", USER_HOME, "/.local/R")); # for unix local account # uses R CMD INSTALL</br>
#### To install MetENP in the folder USER_PWD/R
&#62;USER_PWD=Sys.getenv("PWD");</br>
&#62;devtools::install("MetENP", args = paste0("--library=", USER_PWD, "/R")); # for unix local account # uses R CMD INSTALL</br>
&#62;q()</br>

One-liner for the above:</br>
&#62;library("devtools"); USER_HOME=Sys.getenv("HOME");devtools::install("MetENP", args = paste0("--library=", USER_HOME, "/.local/R")); # for unix local account # uses R CMD INSTALL</br>
&#62;library("devtools"); USER_PWD=Sys.getenv("PWD");devtools::install("MetENP", args = paste0("--library=", USER_PWD, "/R")); # for unix local account # uses R CMD INSTALL</br>

#### # if all went well, this would have installed MetENP in /home/username/.local/R (or current-folder/R)
$ ls -al /home/username/.local/R</br>
 
&#35; to check if MetENP can be loaded</br>
$R</br>
&#35; modify .libPaths so that it can find R package MetENP</br>
&#62;USER_HOME=Sys.getenv("HOME");</br>
&#62;USER_PWD=Sys.getenv("PWD");</br>
&#62;.libPaths( c( .libPaths(), paste0(USER_HOME, "/.local/R") )); # since MetENP installed in user area, need to include that in path</br>
&#62;#.libPaths( c( .libPaths(), paste0(USER_PWD, "/R") )); # If installed in the R folder inside PWD</br>
&#62;library("MetENP") # should load without errors</br>
&#35; Now ready to run jupyter, being in a folder containing *.ipynb file, e.g., </br>
/path-to-MetENP-folder/vignettes/</br>
$ jupyter notebook --ip=123.249.124.012 --port=8080</br>
Go to webpage listed and open a MetENP jupyter notebook </br>
Near top in that file, insert the lines, or some of these lines to set .libPaths and load MetENP R library.</br>
USER_HOME=Sys.getenv("HOME");</br>
USER_PWD=Sys.getenv("PWD");</br>
.libPaths( c( .libPaths(), paste0(USER_HOME, "/.local/R") ))</br>
#.libPaths( c( .libPaths(), paste0(USER_PWD, "/../../R") ))</br>
library("MetENP") # Now, it should load without errors</br>

## (3) For developers: more advanced features
#### # To make changes in the MetENP functions, such as to adjust plot features [2023/02/09]
When possible, it is best to adjust them using arguments passed to function. If that is not possible, e.g., because the corresponding variable cannot be set through the arguments passed, then make a copy of the relevant R file in the folder ..../MetENP/R, edit the file/function suitably and save. Then source the file, e.g.:</br>
&#62;source("PATH-TO-MetENP-FOLDER/R/significant_met_edited.R");</br>
You may also have to load some libraries (see the lines with @importFrom in the R file) if such error messages are displayed. Once the edits are final, the edited file can be renamed (please keep a backup of the original file first) to the original file name and the MetENP package can be rebuilt using the steps above.

#### To update the documentation automatically after significant change such as adding or removing arguments to a function [2024/03/27]
Assuming that roxygen tags are also updated accordingly for adding/removing arguments (e.g., in the line(s) #'@param and the line after #'@examples), then use the commmand:

&#62;devtools::document('MetENP')

This will update the .rd files inside MetENP/man; reinstall the MetENP package again using the procedure explained above.

#### To check the package structure, etc. [2024/03/27]

Use the command rcmdcheck<br>
&#62;rcmdcheck("MetENP")

to identify errors/warning/notes.



