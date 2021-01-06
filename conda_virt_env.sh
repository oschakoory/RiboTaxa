#!/bin/sh

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" 
echo "Creating virtual environment and installing RiboTaxa... " `date` 
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" 
echo ""

conda env create -n RiboTaxa_py27 --file RiboTaxa_py27_requirements.txt

conda install -n RiboTaxa_py27 -c bioconda/label/cf201901 multiqc

conda env create -n RiboTaxa_py36 --file RiboTaxa_py36_requirements.yml

conda install -n RiboTaxa_py36 -c bioconda matam
