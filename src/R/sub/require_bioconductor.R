

#source("https://bioconductor.org/biocLite.R")

#biocLite(suppressUpdates = TRUE)

#biocLite("ShortRead", suppressUpdates = TRUE)

if (
  !requireNamespace(
    "BiocManager",
    quietly = TRUE
    )
  )
  install.packages("BiocManager")

# BiocManager::install("phyloseq")

require(BiocManager)
packages <- c(
#  "DECIPHER",
  "msa",
#  "dada2",
#  "phyloseq",
  "Biostrings",
#  "ggtree",
#  "ComplexHeatmap",
  "rjson",
  "limma",
  "metagenomeSeq",
  "interactiveDisplay",
  "microbiome"
  ) 

new.packages <-
  packages[!(packages %in% installed.packages())] 
# installed.packages() returns installed packages 

if(length(new.packages) > 0){ 
  BiocManager::install(new.packages, suppressUpdates = TRUE)
  }

#require("DECIPHER")
require("msa")
require("microbiome")
#require("dada2")
#require("phyloseq")
require("Biostrings")
#require("ggtree")
#require("ComplexHeatmap")
require("biomformat")
require("limma")
require("metagenomeSeq")
require("interactiveDisplay")

if(Bibtex){
  write(toBibtex(citation()),file="CRAN")
  for(i in 1:length(packages)){
    write(
      toBibtex(
        citation(packages[i])
        ),
      file=sprintf("../Biblio/%s%s.bib",packages[i],"_CRAN"))
  }
}

