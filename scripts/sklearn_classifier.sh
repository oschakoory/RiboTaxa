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
conda activate RiboTaxa_py36
#echo "Qiime2 virtual environment has been activated successfully..." | tee /dev/fd/3

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo ">Taxonomic classification using sklearn classifier of QIIME2020.8 ..." `date`| tee /dev/fd/3
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo ""

#echo "Set up directories ..." | tee /dev/fd/3

CONFIG_PATH=$1
CONFIG="${CONFIG_PATH[@]}"

SHORTNAME=$2
#SHORTNAME=$(basename ""${NAME[@]}"" | sed 's/_R1.'$FORMAT'//') 
echo "SHORTNAME = " >&2
printf '%s\n' "$SHORTNAME" >&2

RiboTaxa_DIR=$(awk '/^RiboTaxa_DIR/{print $3}' "${CONFIG}")

#Set up database directory
SKLEARN_DB=$(awk '/^SKLEARN_DB/{print $3}' "${CONFIG}")
#echo "Sklearn Database = $SKLEARN_DB" | tee /dev/fd/3

#DIR=$(awk '/^OUTPUT/{print $3}' "${CONFIG}")
#mkdir -p "$DIR/data"
#export TMPDIR="$DIR/data"

#echo $TMPDIR | tee /dev/fd/3

#set up output directory
OUTPUT=$(awk '/^OUTPUT/{print $3}' "${CONFIG}")
#mkdir -p "$OUTPUT/Taxonomy"

RESULTS="$OUTPUT/$SHORTNAME"
mkdir -p $RESULTS
mkdir -p "$RESULTS/Taxonomy"

mkdir -p "$RESULTS/data"
export TMPDIR="$RESULTS/data"

THREAD=$(awk '/^THREAD/{print $3}' "${CONFIG}")
#echo "Number of threads used = $THREAD" | tee /dev/fd/3

CONFIDENCE=$(awk '/^CONFIDENCE/{print $3}' "${CONFIG}")
#echo "Confidence Threshold for taxonomic classification = $CONFIDENCE" | tee /dev/fd/3

BATCH=$(awk '/^BATCH/{print $3}' "${CONFIG}")
#echo "Number of reads to process per batch = $BATCH" | tee /dev/fd/3

TAXONOMY="$RESULTS"/Taxonomy/"$SHORTNAME"_SSU_taxonomy_abundance.tsv


if [ -f "$TAXONOMY" ];
then
  echo "Taxonomoy file : '${TAXONOMY}' is present."
else
echo "Importing data ..." | tee /dev/fd/3
echo "command line :" | tee /dev/fd/3
echo "qiime tools import \
  --input-path $RESULTS/SSU_sequences/"$SHORTNAME"_SSU_sequences.fasta \
  --output-path $RESULTS/Taxonomy/"$SHORTNAME"_SSU_sequences_qiime2.qza \
  --type FeatureData[Sequence]" | tee /dev/fd/3
  
#--------Importation des données dans QIIME 2
qiime tools import \
  --input-path "$RESULTS"/SSU_sequences/"$SHORTNAME"_SSU_sequences.fasta \
  --output-path "$RESULTS"/Taxonomy/"$SHORTNAME"_SSU_sequences_qiime2.qza \
  --type FeatureData[Sequence]


echo ">Classifying reconstructed sequences using sklearn_classifer..." | tee /dev/fd/3

echo "command line :" | tee /dev/fd/3
echo "qiime feature-classifier classify-sklearn \
  --i-classifier $SKLEARN_DB \
  --i-reads $RESULTS/Taxonomy/"$SHORTNAME"_SSU_sequences_qiime2.qza  \
  --o-classification $RESULTS/Taxonomy/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy \
  --p-n-jobs 1 \
  --p-confidence $CONFIDENCE \
  --p-reads-per-batch $BATCH \
  --verbose \
  --p-read-orientation same" | tee /dev/fd/3

qiime feature-classifier classify-sklearn \
  --i-classifier $SKLEARN_DB \
  --i-reads "$RESULTS"/Taxonomy/"$SHORTNAME"_SSU_sequences_qiime2.qza  \
  --o-classification "$RESULTS"/Taxonomy/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy \
  --p-n-jobs 1 \
  --p-confidence $CONFIDENCE \
  --p-reads-per-batch $BATCH \
  --verbose \
  --p-read-orientation same

rm -rf taxonomy.tsv

#handling taxonomy table from qiime2
for zipfile in `ls "$RESULTS"/Taxonomy/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy.qza` 
do 
  unzip -j "$zipfile" '*/*/*.tsv'
  mv taxonomy.tsv "$SHORTNAME"_taxonomy.tsv
done

mv "$SHORTNAME"_taxonomy.tsv "$RESULTS"/Taxonomy/"$SHORTNAME"_taxonomy.tsv

if [[ ${SKLEARN_DB[@]} == *GTDB* ]]; then
    #echo "$FILE is gzipped"
    cat "$RESULTS"/Taxonomy/"$SHORTNAME"_taxonomy.tsv | sed 1d |sort -k1 | tr ' ' '_' | awk '{ print $1,$2 }'| tr ' ' \\t  > "$RESULTS"/Taxonomy/"$SHORTNAME"_sorted_taxonomy.tsv
else
   cat "$RESULTS"/Taxonomy/"$SHORTNAME"_taxonomy.tsv | sed 1d |sort -k1 | tr -d ' ' | awk '{ print $1,$2 }'| tr ' ' \\t  > "$RESULTS"/Taxonomy/"$SHORTNAME"_sorted_taxonomy.tsv
fi

#cat "$RESULTS"/Taxonomy/"$SHORTNAME"_taxonomy.tsv | sed 1d |sort -k1 | tr ' ' ';' | awk '{ print $1,$2 }'| tr ' ' \\t  > "$RESULTS"/Taxonomy/"$SHORTNAME"_sorted_taxonomy.tsv
cat "$RESULTS"/Taxonomy/"$SHORTNAME"_sorted_taxonomy.tsv | tr ';' \\t |sed 's/\t\+/\t/g;s/^\t//' > "$RESULTS"/Taxonomy/"$SHORTNAME"_sorted_tab_taxonomy.tsv
rm "$RESULTS"/Taxonomy/"$SHORTNAME"_sorted_taxonomy.tsv
awk -v c=8 'BEGIN{FS=OFS="\t"} {for(i=NF+1; i<=c; i++) $i="\t"} 1' "$RESULTS"/Taxonomy/"$SHORTNAME"_sorted_tab_taxonomy.tsv > "$RESULTS"/Taxonomy/"$SHORTNAME"_sorted_taxonomy.tsv

cat "$RESULTS"/Taxonomy/"$SHORTNAME"_taxonomy.tsv | sed 1d |sort -k1 | tr ' ' ';' | awk '{ print $1,$3 }' > "$RESULTS"/Taxonomy/"$SHORTNAME"_sorted_confidence.tsv

join "$RESULTS"/Taxonomy/"$SHORTNAME"_sorted_confidence.tsv "$RESULTS"/Taxonomy/"$SHORTNAME"_sorted_taxonomy.tsv  |sort -k1| tr ' ' \\t > "$RESULTS"/Taxonomy/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy.tsv
rm "$RESULTS"/Taxonomy/"$SHORTNAME"_sorted_taxonomy.tsv 
rm "$RESULTS"/Taxonomy/"$SHORTNAME"_sorted_confidence.tsv
rm "$RESULTS"/Taxonomy/"$SHORTNAME"_sorted_tab_taxonomy.tsv

cat "$RESULTS"/SSU_sequences/"$SHORTNAME"_Abundance.tsv | sed 1d |sort -k1| tr ' ' \\t > "$RESULTS"/Taxonomy/"$SHORTNAME"_abundance.tsv


join "$RESULTS"/Taxonomy/"$SHORTNAME"_abundance.tsv "$RESULTS"/Taxonomy/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy.tsv |tr  ' ' \\t | awk -v c=12 'BEGIN{FS=OFS="\t"} {for(i=NF+1; i<=c; i++) $i="\t"} 1' | awk '{ print $1,$6,$7,$8,$9,$10,$11,$12,$5,$2,$3,$4}' |awk 'BEGIN{print "ID\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tConfidence\tLength\tAssigned_reads\tRelative_abundance "}1' |tr  ' ' \\t  > "$RESULTS"/Taxonomy/"$SHORTNAME"_SSU_taxonomy_abundance.tsv


#cat "$RESULTS"/Taxonomy/taxonomy.tsv | sed 1d |sort -k1| tr ';' \\t > "$RESULTS"/Taxonomy/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy.tsv

#cat "$RESULTS"/SSU_sequences/"$SHORTNAME"_Abundance.tsv | sed 1d |sort -k1| tr ' ' \\t > "$RESULTS"/Taxonomy/abundance.tsv

#handling length of reconstructed sequences
#cat "$RESULTS"/SSU_sequences/"$SHORTNAME"_covstats.txt| sed 1d | awk '{ print $1,$3 }' |sort -k1 -k2 | tr ' ' \\t > "$RESULTS"/SSU_sequences/"$SHORTNAME"_covstats_readlength.tsv


#handling relative abundance table from bbmap
#cat "$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats.txt | sed 1d | awk '{ print $1,$8 }' | sort -k1 -k2 | tr ' ' \\t > "$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats_readsCount.tsv

#cat "$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats_readsCount.tsv | tr ' ' \\t > "$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats_readsCount1.tsv


#total=$(awk '{s+=$2}END{print s}' "$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats_readsCount.tsv)

#awk -v total=$total '{ printf ("%s\t%s\t%.6f\n", $1, $2, ($2/total)*100)}' "$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats_readsCount.tsv > "$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats_RA.tsv

#join "$RESULTS"/SSU_sequences/"$SHORTNAME"_covstats_readlength.tsv "$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats_RA.tsv |tr  ' ' \\t > "$RESULTS"/Taxonomy/"$SHORTNAME"_SSU_length_abundance.tsv

#join "$RESULTS"/Taxonomy/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy.tsv "$RESULTS"/Taxonomy/abundance.tsv |tr  ' ' \\t |awk '!($11==0){print}'|awk 'BEGIN{print "ID\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tConfidence\tLength\tAssigned_reads\tRelative_abundance "}1' > "$RESULTS"/Taxonomy/"$SHORTNAME"_SSU_taxonomy_abundance.tsv

#cat "$RESULTS"/Taxonomy/"$SHORTNAME"_SSU_taxonomy_abundance.tsv | awk '!($11==0){print}' > "$RESULTS"/Taxonomy/"$SHORTNAME"_filtered_SSU_taxonomy_abundance.tsv

#cleaning
#rm "$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats_readsCount.tsv
#rm "$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats_readsCount1.tsv
#rm "$RESULTS"/Taxonomy/"$SHORTNAME"_SSU_taxonomy_abundance_renamed.tsv
#rm "$RESULTS"/SSU_sequences/"$SHORTNAME"_covstats_readlength.tsv 
#rm "$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats_RA.tsv
#rm "$RESULTS"/Taxonomy/"$SHORTNAME"_SSU_length_abundance.tsv
rm "$RESULTS"/Taxonomy/"$SHORTNAME"_renamed_16S18S_recons_qiime2_taxonomy.tsv
rm "$RESULTS"/Taxonomy/"$SHORTNAME"_taxonomy.tsv
rm "$RESULTS"/Taxonomy/"$SHORTNAME"_SSU_sequences_qiime2.qza
rm "$RESULTS"/Taxonomy/"$SHORTNAME"_abundance.tsv

rm -rf "$RESULTS"/data

fi 

echo "Finalising abundance table..." | tee /dev/fd/3


cp $RiboTaxa_DIR/scripts/PostRiboTaxa.sh "$RESULTS"/Taxonomy

cd "$RESULTS"/Taxonomy

./PostRiboTaxa.sh "$SHORTNAME"

rm PostRiboTaxa.sh

echo ">Taxonomic classification using sklearn_classifer ends successfully on : "`date` | tee /dev/fd/3

conda deactivate
#echo "Qiime2 virtual environment has been deactivated successfully..." | tee /dev/fd/3

