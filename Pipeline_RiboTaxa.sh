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

Jobname=$(awk '/^JOB_NAME/{print $3}' "${CONFIG}")

exec 3>&1 1>RiboTaxa_"$Jobname".log 2>&3 2>RiboTaxa_"$Jobname".stderr
echo " "
echo "RiboTaxa -- A complete pipeline from raw metagenomics to species-level identification" | tee /dev/fd/3
echo "By Oshma Chakoory, Sophie Marre & Pierre Peyret" | tee /dev/fd/3
echo "University Clermont Auvergne, France " | tee /dev/fd/3
echo "Version: 1.4" | tee /dev/fd/3

echo "This program is distributed under the AGPL-3.0 License. See LICENSE for more information." | tee /dev/fd/3


DATA_DIR=$(awk '/^DATA_DIR/{print $3}' "${CONFIG}")
#echo "Data path = $DATA_DIR" | tee /dev/fd/3

OUTPUT=$(awk '/^OUTPUT/{print $3}' "${CONFIG}")

RiboTaxa_DIR=$(awk '/^RiboTaxa_DIR/{print $3}' "${CONFIG}")

FORMAT=$(awk '/^FORMAT/{print $3}' "${CONFIG}")
#echo $FORMAT

for NAME in `ls "$DATA_DIR"/*_R1.$FORMAT | sed 's/_R1.'$FORMAT'//'` 
do
SHORTNAME=$(basename ""${NAME[@]}"") 
#echo $SHORTNAME
echo "" | tee /dev/fd/3
echo "***********************************************************************************************" | tee /dev/fd/3
echo ">Running RiboTaxa on $NAME"  | tee /dev/fd/3
echo "***********************************************************************************************" | tee /dev/fd/3
echo "" | tee /dev/fd/3
#run RiboTaxa.sh script to perform
### Perform a quality control on the (meta)genomics data using BBTOOLS
### Filter 16S/18S reads using SORTMERNA
### Reconstruct full length 16S/18S rRNA sequences using EMIRGE

source "$RiboTaxa_DIR"/RiboTaxa.sh $CONFIG_PATH 

#run sklearn_classfier.sh script to perform
### Taxonomic classification of full length reconstrcuted sequences using sklearn classifier of Qiime2

source "$RiboTaxa_DIR"/sklearn_classifier.sh $CONFIG_PATH 

done

mv RiboTaxa_"$Jobname".log "$OUTPUT"
mv RiboTaxa_"$Jobname".stderr "$OUTPUT"
#mv "$OUTPUT"/output_MetaRib "$OUTPUT"/SSU_sequences