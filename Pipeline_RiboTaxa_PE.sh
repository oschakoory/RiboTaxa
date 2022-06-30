#!/bin/sh

# Program configuration
__author__='Oshma Chakoory'
__email__='oshma.chakoory@uca.fr'
__credits__=["Oshma"]
__status__='Development'
__version__='1.4'


# Handling errors
#set -x # debug mode on
set -o errexit # ensure script will stop in case of ignored error
set -e


CONFIG_PATH=$1
CONFIG="${CONFIG_PATH[@]}"

#Jobname=$(awk '/^JOB_NAME/{print $3}' "${CONFIG}")


DATA_DIR=$(awk '/^DATA_DIR/{print $3}' "${CONFIG}")
#echo "Data path = $DATA_DIR" | tee /dev/fd/3

OUTPUT=$(awk '/^OUTPUT/{print $3}' "${CONFIG}")

RiboTaxa_DIR=$(awk '/^RiboTaxa_DIR/{print $3}' "${CONFIG}")

FORMAT=$(awk '/^FORMAT/{print $3}' "${CONFIG}")
#echo $FORMAT

for NAME in `ls "$DATA_DIR"/*_R1."$FORMAT"` 
do
SHORTNAME=$(basename ""${NAME[@]}"" | sed 's/_R1.'$FORMAT'//') 
#SHORTNAME=$(basename ""${NAME[@]}"") 
#echo $SHORTNAME

RESULTS="$OUTPUT/$SHORTNAME"
mkdir -p $RESULTS

exec 3>&1 1>"$RESULTS"/RiboTaxa_"$SHORTNAME".log 2>&3 2>"$RESULTS"/RiboTaxa_"$SHORTNAME".stderr
echo " "
echo "RiboTaxa -- A complete pipeline from raw metagenomics to species-level identification" | tee /dev/fd/3
echo "By Oshma Chakoory, Sophie Marre & Pierre Peyret" | tee /dev/fd/3
echo "University Clermont Auvergne, France " | tee /dev/fd/3
echo "Version: 1.4" | tee /dev/fd/3

echo "This program is distributed under the AGPL-3.0 License. See LICENSE for more information." | tee /dev/fd/3

echo "" | tee /dev/fd/3
echo "***********************************************************************************************" | tee /dev/fd/3
echo ">Running RiboTaxa on $NAME"  | tee /dev/fd/3
echo "***********************************************************************************************" | tee /dev/fd/3
echo "" | tee /dev/fd/3
#run RiboTaxa.sh script to perform
### Perform a quality control on the (meta)genomics data using BBTOOLS
### Filter 16S/18S reads using SORTMERNA
### Reconstruct full length 16S/18S rRNA sequences using EMIRGE

source "$RiboTaxa_DIR"/scripts/RiboTaxa_PE.sh $CONFIG_PATH $SHORTNAME

#run sklearn_classfier.sh script to perform
### Taxonomic classification of full length reconstrcuted sequences using sklearn classifier of Qiime2

source "$RiboTaxa_DIR"/scripts/sklearn_classifier.sh $CONFIG_PATH $SHORTNAME

done

#mv RiboTaxa_"$Jobname".log "$RESULTS"
#mv RiboTaxa_"$Jobname".stderr "$RESULTS"
#mv "$RESULTS"/output_MetaRib "$RESULTS"/SSU_sequences
