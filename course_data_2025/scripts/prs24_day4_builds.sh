#!/bin/bash

## List of dependencies for Day 4
## Installing cmake
sudo apt-get install build-essential libssl-dev
cd /tmp
wget https://github.com/Kitware/CMake/releases/download/v3.20.0/cmake-3.20.0.tar.gz
gunzip  Downloads/cmake-3.29.4.tar.gz 
sudo apt-get -y install cmake

##R packages
echo "install.packages(“Rcpp”) # Required package before “reshape2 package” installation
install.packages(“reshape2”)
install.packages(“viridisLite”) # Required before installing “viridis”
install.packages(“viridis”)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GenomicRanges")
#install Bioconduct # Required package before “GenomicRanges” package installation
#install.packages(“GenomicRanges”)
install.packages(“dplyr”)
install.packages(“tidyr”)
install.packages(“vcfR”)
install.packages(“memuse”)
#install.packages(“panelr”)
install.packages(“data.table”)
install.packages(“panelr”)" > rscript_packages.R

#run R install
Rscript rscript_packages.R

#intalling bcf tools

sudo apt-get install autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 -O samtools.tar.bz2
tar -xjvf samtools.tar.bz2
cd samtools-{version}
cd samtools-1.3.1/
make
sudo make prefix=/usr/local/bin install
wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 -O bcftools.tar.bz2
tar -xjvf bcftools.tar.bz2
cd bcftools-1.3.1/
make
sudo make prefix=/usr/local/bin install
sudo ln -s /usr/local/bin/bin/bcftools /usr/bin/bcftools

