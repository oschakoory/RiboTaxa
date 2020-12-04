#!/bin/sh

conda env create -n MicroTaxa_py27 --file MicroTaxa_py27_requirements.txt

conda install -n MicroTaxa_py27 -c bioconda/label/cf201901 multiqc

conda env create -n MicroTaxa_py36 --file MicroTaxa_py36_requirements.yml

conda install -n MicroTaxa_py36 -c bioconda matam
