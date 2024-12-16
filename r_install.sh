#!/bin/bash

# Install Every Package into ./R_PACKAGES/
while read r; do
    R -e "install.packages(${r}, repos='http://cran.us.r-project.org', lib='./R_PACKAGES/')"
done <R_requirements.txt

while read bio; do
    R -e "BiocManager::install(${bio}, force=TRUE, lib='./R_PACKAGES/')"
done <R_BiocManager.txt