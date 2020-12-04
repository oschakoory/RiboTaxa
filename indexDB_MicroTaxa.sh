#!/bin/sh

echo -n "" > indexdb_MicroTaxa.log
exec 3>&1 1>> indexdb_MicroTaxa.log 2>&3

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

MicroTaxa_DIR=$(awk '/^MicroTaxa_DIR/{print $3}' "${CONFIG}")

DB_DIR=$(awk '/^DB_DIR/{print $3}' "${CONFIG}")
echo "DB_DIR_path = $DB_DIR" | tee /dev/fd/3

NAME=($(ls "$DB_DIR" | sed 's/.fasta//')) 
#echo "name = " | tee /dev/fd/3
#printf '%s\n' "${NAME[@]}" | tee /dev/fd/3

SHORTNAME=$(basename ""${NAME[@]}"")  
#echo "shortname = $SHORTNAME" | tee /dev/fd/3

OUTPUT=$(awk '/^OUTPUT/{print $3}' "${CONFIG}")
echo "OUTPUT_path = $OUTPUT" | tee /dev/fd/3

mkdir -p "$OUTPUT/sortmerna_indexed_DB"
mkdir -p "$OUTPUT/bowtie_indexed_DB"

THREAD=$(awk '/^THREAD/{print $3}' "${CONFIG}")

CLUSTER_ID=$(awk '/^CLUSTER_ID/{print $3}' "${CONFIG}")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#			Indexing Databases for sortmeRNA
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

conda activate MicroTaxa_py36
#echo "MicroTaxa virtual environment has been activated successfully..." | tee /dev/fd/3

echo "Starting indexing database for sortmerna... " `date` | tee /dev/fd/3

matam_db_preprocessing.py -i "$DB_DIR"  \
	--clustering_id_threshold $CLUSTER_ID \
	-d "$OUTPUT"/sortmerna_indexed_DB \
	--out_db_name "$SHORTNAME"_indexed \
	--cpu "$THREAD" \
	-v

echo "Finished indexing database for sortmerna... " `date` | tee /dev/fd/3

conda deactivate
#echo "MicroTaxa virtual environment has been deactivated successfully..." | tee /dev/fd/3


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#			Fixing and Indexing Database for EMIRGE
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

conda activate MicroTaxa_py27
#echo "MicroTaxa virtual environment has been activated successfully..." | tee /dev/fd/3

echo "Remove ambiguous characters and replacing non-standard base characters in reference database..." `date` | tee /dev/fd/3

python "$MicroTaxa_DIR"/fix_nonstandard_chars.py < "$OUTPUT"/sortmerna_indexed_DB/*.clustered.fasta* > "$OUTPUT"/bowtie_indexed_DB/"$SHORTNAME"_clustered_fixed.fasta 

echo "Finished fixing database for emirge..."`date` | tee /dev/fd/3


echo "starting indexing the newly fixed and clustered database with bowtie..."`date` | tee /dev/fd/3
#Create bowtie-index for your new reference database
bowtie-build "$OUTPUT"/bowtie_indexed_DB/"$SHORTNAME"_clustered_fixed.fasta  "$OUTPUT"/bowtie_indexed_DB/"$SHORTNAME"_bowtie_indexed \
	-p "$THREAD"

echo "Finished indexing the newly fixed and clustred database with bowtie..."`date` | tee /dev/fd/3

echo "Saving results..."

echo "Job ends successfully on : "`date` | tee /dev/fd/3

conda deactivate
#echo "MicroTaxa virtual environment has been deactivated successfully..." | tee /dev/fd/3
