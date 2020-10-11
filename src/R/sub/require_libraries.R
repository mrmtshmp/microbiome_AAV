# options(BioC_mirror="https://bioconductor.org/")
# source("https://bioconductor.org/biocLite.R")


packages <- c(
  'apex',
  'beeswarm',
  #  "brms",
  "coin",
  "caret",
  'dplyr', # progress bar
  'ExploratoryDataAnalysis',
  #  'Formula',
  "gridExtra",
  'GGally',
  'ggplot2',
  'gplots',
  "ggrepel",
#  'GMD',
  "GUniFrac",
  "Matrix",
  "metagenomeSeq",
  "phylotools",
  "phyloseq",
  #  'pvclust',
  'pander',
  "phangorn",
  'plyr',  # progress bar
  "pscl",    # zeroinfl(...)
  "RColorBrewer",
  'readxl',
  'reshape2',
  "rms",
  'robustbase',
  "seqinr",
  "stringr",
  'tidyr',
  'tidyverse',
  "vegan"
)

new.packages <-
  packages[!(packages %in% installed.packages())] 
# installed.packages() returns installed packages 

if(length(new.packages) > 0){ 
  install.packages(new.packages, repos=c('http://cran.us.r-project.org'))
#  biocLite(new.packages)
}
require('beeswarm')
require("coin")
require("caret")
require("gridExtra")
require("ggrepel")
require("ExploratoryDataAnalysis")
require("phylotools")
require("apex")
require("seqinr")
require("stringr")
require("phangorn")
require("phyloseq")
require("GUniFrac")
require("rms")
#require("pscl")
#require("brms")
require("RColorBrewer")
require('robustbase')
require("Matrix")
require("metagenomeSeq")
require('readxl')
require('plyr')  # progress bar
require('dplyr') # progress bar
require('tidyr')
require('tidyverse')
require('GGally')
require('ggplot2')
require('gplots')
#require('GMD')
require('vegan')
#require('pvclust')
#require('reshape2')
#require('pander')


if(Bibtex){
  write(toBibtex(citation()),file="CRAN")
  for(i in 1:length(packages)){
    write(toBibtex(citation(packages[i])),file=sprintf("../Biblio/%s%s.bib",packages[i],"_CRAN"))
  }
}

