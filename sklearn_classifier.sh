#!/bin/bash

__author__='Oshma Chakoory'
__email__='oshma.chakoory@uca.fr'
__credits__=["Oshma"]
__status__='Development'
__version__='0.0.1'

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
  --input-path "$OUTPUT"/output_emirge/"$SHORTNAME"_renamed_16S18S_recons.fasta \
  --output-path "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_renamed_16S18S_recons_qiime2.qza \
  --type FeatureData[Sequence]

echo "Classifying reconstructed sequences using sklearn_classifer..." | tee /dev/fd/3
qiime feature-classifier classify-sklearn \
  --i-classifier $SKLEARN_DB \
  --i-reads "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_renamed_16S18S_recons_qiime2.qza  \
  --o-classification "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy \
  --p-n-jobs 1 \
  --p-confidence $CONFIDENCE \
  --p-reads-per-batch $BATCH \
  --verbose

rm -rf taxonomy.tsv

for zipfile in `ls "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy.qza` 
do 
  unzip -j "$zipfile" '*/*/*.tsv'
done

mv taxonomy.tsv "$OUTPUT"/tax_class_sklearn_qiime2/taxonomy.tsv

cat "$OUTPUT"/tax_class_sklearn_qiime2/taxonomy.tsv | sed 1d |sort -k1,1n -k2,2 | tr ';' \\t > "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy.tsv
cat  "$OUTPUT"/output_emirge/"$SHORTNAME"_renamed_16S18S_recons.fasta |grep '^>'|tr -d '>' | sort -n | awk '{sub("NormPrior=", "",$4);print}'|tr ' ' \\t | cut -f 1,4 > "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_emirge_abundance.tsv
join "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy.tsv "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_emirge_abundance.tsv |tr  ' ' \\t > "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_SSU_taxonomy_abundance.tsv

rm "$OUTPUT"/tax_class_sklearn_qiime2/"$SHORTNAME"_emirge_abundance.tsv
rm "$OUTPUT"/tax_class_sklearn_qiime2/taxonomy.tsv

echo "Taxonomic classification using sklearn_classifer ends successfully on : "`date` | tee /dev/fd/3

conda deactivate
#echo "Qiime2 virtual environment has been deactivated successfully..." | tee /dev/fd/3

