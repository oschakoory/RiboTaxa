#!/bin/sh

exec 3>&1 1>MicroTaxa.log 2>&3

# Program configuration
__author__='Oshma Chakoory'
__email__='oshma.chakoory@uca.fr'
__credits__=["Oshma"]
__status__='Development'
__version__='0.0.1'


# Handling errors
#set -x # debug mode on
set -o errexit # ensure script will stop in case of ignored error

CONFIG_PATH=$1
CONFIG="${CONFIG_PATH[@]}"

DATA_DIR=$(awk '/^DATA_DIR/{print $3}' "${CONFIG}")
#echo "Data path = $DATA_DIR" | tee /dev/fd/3

OUTPUT=$(awk '/^OUTPUT/{print $3}' "${CONFIG}")

MicroTaxa_DIR=$(awk '/^MicroTaxa_DIR/{print $3}' "${CONFIG}")

for NAME in `ls "$DATA_DIR"/*_1.fastq | sed 's/_1.fastq//'` 
do
SHORTNAME=$(basename ""${NAME[@]}"") 
echo "" | tee /dev/fd/3
echo "***********************************************************************************************" | tee /dev/fd/3
echo "               Running MicroTaxa on $NAME"  | tee /dev/fd/3
echo "***********************************************************************************************" | tee /dev/fd/3
echo "" | tee /dev/fd/3
#run MicroTaxa.sh script to perform
### Perform a quality control on the (meta)genomics data using BBTOOLS
### Filter 16S/18S reads using SORTMERNA
### Reconstruct full length 16S/18S rRNA sequences using EMIRGE

source "$MicroTaxa_DIR"/MicroTaxa.sh $CONFIG_PATH 

#run sklearn_classfier.sh script to perform
### Taxonomic classification of full length reconstrcuted sequences using sklearn classifier of Qiime2

source "$MicroTaxa_DIR"/sklearn_classifier.sh $CONFIG_PATH 

done

mv MicroTaxa.log "$OUTPUT"