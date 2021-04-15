#!/bin/bash

__author__='Oshma Chakoory'
__email__='oshma.chakoory@uca.fr'
__credits__=["Oshma"]
__status__='Development'
__version__='1.4'

# Handling errors
#set -x # debug mode on
set -o errexit # ensure script will stop in case of ignored error

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#		Taxonomic classification using sklearn classifier of QIIME2020.8
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#activate virtual environment for RiboTaxa
conda activate RiboTaxa_py36
#echo "Qiime2 virtual environment has been activated successfully..." | tee /dev/fd/3

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo "Taxonomic classification using sklearn classifier of QIIME2020.8 ..." `date`| tee /dev/fd/3
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo ""

#echo "Set up directories ..." | tee /dev/fd/3

CONFIG_PATH=$1
CONFIG="${CONFIG_PATH[@]}"

RiboTaxa_DIR=$(awk '/^RiboTaxa_DIR/{print $3}' "${CONFIG}")

#Set up database directory
SKLEARN_DB=$(awk '/^SKLEARN_DB/{print $3}' "${CONFIG}")
echo "Sklearn Database = $SKLEARN_DB" | tee /dev/fd/3

DIR=$(awk '/^OUTPUT/{print $3}' "${CONFIG}")
mkdir -p "$DIR/data"
export TMPDIR="$DIR/data"

echo $TMPDIR | tee /dev/fd/3

#set up output directory
OUTPUT=$(awk '/^OUTPUT/{print $3}' "${CONFIG}")
mkdir -p "$OUTPUT/tax_class_sklearn_qiime2"

THREAD=$(awk '/^THREAD/{print $3}' "${CONFIG}")
#echo "Number of threads used = $THREAD" | tee /dev/fd/3

CONFIDENCE=$(awk '/^CONFIDENCE/{print $3}' "${CONFIG}")
#echo "Confidence Threshold for taxonomic classification = $CONFIDENCE" | tee /dev/fd/3

BATCH=$(awk '/^BATCH/{print $3}' "${CONFIG}")
#echo "Number of reads to process per batch = $BATCH" | tee /dev/fd/3

echo "Importing data ..." | tee /dev/fd/3
#--------Importation des données dans QIIME 2
qiime tools import \
  --input-path "$OUTPUT"/SSU_sequences/"$SHORTNAME"_SSU_sequences.fasta \
  --output-path "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_SSU_sequences_qiime2.qza \
  --type FeatureData[Sequence]

echo "Classifying reconstructed sequences using sklearn_classifer..." | tee /dev/fd/3
qiime feature-classifier classify-sklearn \
  --i-classifier $SKLEARN_DB \
  --i-reads "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_SSU_sequences_qiime2.qza  \
  --o-classification "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy \
  --p-n-jobs 1 \
  --p-confidence $CONFIDENCE \
  --p-reads-per-batch $BATCH \
  --verbose

rm -rf taxonomy.tsv

#handling taxonomy table from qiime2
for zipfile in `ls "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy.qza` 
do 
  unzip -j "$zipfile" '*/*/*.tsv'
done

mv taxonomy.tsv "$OUTPUT"/tax_class_sklearn_qiime2/taxonomy.tsv

cat "$OUTPUT"/tax_class_sklearn_qiime2/taxonomy.tsv | sed 1d |sort -k1 -k2 | tr ';' \\t > "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy.tsv

#handling length of reconstructed sequences
cat "$OUTPUT"/SSU_sequences/"$SHORTNAME"_covstats.txt| sed 1d | awk '{ print $1,$3 }' |sort -k1 -k2 | tr ' ' \\t > "$OUTPUT"/SSU_sequences/"$SHORTNAME"_covstats_readlength.tsv


#handling relative abundance table from bbmap
cat "$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats.txt | sed 1d | tr ',' \\t | awk '{ print $1,$8 }' | sort -k1 -k2 > "$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats_readsCount.tsv

cat "$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats_readsCount.tsv | tr ' ' \\t > "$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats_readsCount1.tsv


total=$(awk '{s+=$2}END{print s}' "$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats_readsCount1.tsv)

awk -v total=$total '{ printf ("%s\t%s\t%.2f\n", $1, $2, ($2/total)*100)}' "$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats_readsCount1.tsv > "$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats_RA.tsv

join "$OUTPUT"/SSU_sequences/"$SHORTNAME"_covstats_readlength.tsv "$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats_RA.tsv |tr  ' ' \\t > "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_SSU_length_abundance.tsv

join "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy.tsv "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_SSU_length_abundance.tsv |tr  ' ' \\t > "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_SSU_taxonomy_abundance_renamed.tsv

awk 'BEGIN{print "ID\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tConfidence\tLength\tAssigned_reads\tRelative_abundance "}1' "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_SSU_taxonomy_abundance_renamed.tsv > "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_SSU_taxonomy_abundance.tsv

#cleaning
rm "$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats_readsCount.tsv
rm "$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats_readsCount1.tsv
rm "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_SSU_taxonomy_abundance_renamed.tsv
rm "$OUTPUT"/SSU_sequences/"$SHORTNAME"_covstats_readlength.tsv 
rm "$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats_RA.tsv
rm "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_SSU_length_abundance.tsv
rm "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy.tsv

echo "Taxonomic classification using sklearn_classifer ends successfully on : "`date` | tee /dev/fd/3

conda deactivate
#echo "Qiime2 virtual environment has been deactivated successfully..." | tee /dev/fd/3

