#!/bin/sh

# Program configuration
__author__='Oshma Chakoory'
__email__='oshma.chakoory@uca.fr'
__credits__=["Oshma"]
__status__='Development'
__version__='1.5'

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" 
echo "Creating virtual environment and installing RiboTaxa... " `date` 
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" 
echo ""

conda env create -n RiboTaxa_py27 --file RiboTaxa_py27_requirements.yml

conda env create -n mm --file mm.yml

conda install -c bioconda matam -n mm -y

conda env create -n RiboTaxa_py36 --file RiboTaxa_py36_requirements.yml
