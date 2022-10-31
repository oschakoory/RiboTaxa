#!/bin/bash

__author__='Oshma Chakoory'
__email__='oshma.chakoory@uca.fr'
__credits__=["Oshma"]
__status__='Development'
__version__='1.4'

# Handling errors
#set -x # debug mode on
set -o errexit # ensure script will stop in case of ignored error
set -e

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#		Taxonomic classification using sklearn classifier of QIIME2020.8
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#activate virtual environment for RiboTaxa
source activate RiboTaxa_py36
#echo "Qiime2 virtual environment has been activated successfully..." >&2

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >&2
echo ">Grouping multiple taxonomic tables starts..." `date`>&2
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" >&2
echo ""

#echo "Set up directories ..." >&2

INPUT=$1
echo "Input files = " >&2
printf '%s\n' "$INPUT" >&2

#set up output directory
OUTPUT=$2
mkdir -p "$OUTPUT"
echo "Output files = " >&2
printf '%s\n' "$OUTPUT" >&2

./group_taxonomy.sh $INPUT $OUTPUT

echo ">>Grouping multiple taxonomic tables ends successfully on : "`date` >&2

conda deactivate
#echo "Qiime2 virtual environment has been deactivated successfully..." >&2

