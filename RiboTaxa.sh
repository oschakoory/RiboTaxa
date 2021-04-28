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

#activate virtual environment for RiboTaxa
conda activate RiboTaxa_py27
#echo "RiboTaxa virtual environment has been activated successfully..." | tee /dev/fd/3

#echo "Setting up directories for quality control..." | tee /dev/fd/3

#set up 
CONFIG_PATH=$1
CONFIG="${CONFIG_PATH[@]}"

RiboTaxa_DIR=$(awk '/^RiboTaxa_DIR/{print $3}' "${CONFIG}")

#Set up data directory
DATA_DIR=$(awk '/^DATA_DIR/{print $3}' "${CONFIG}")
#echo "Data path = $DATA_DIR" | tee /dev/fd/3

#Set up output directory
OUTPUT=$(awk '/^OUTPUT/{print $3}' "${CONFIG}")
#rm -rf "$OUTPUT"
mkdir -p "$OUTPUT"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#				Filtering 16S/18S reads using sortemrna
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#echo "Setting up directories for sortmerna..." | tee /dev/fd/3
echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo "Filtering 16S/18S using sortmerna starting on : "`date` | tee /dev/fd/3
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo ""

mkdir -p "$OUTPUT/output_sortmerna"

#echo "Setting up indexed database for sortmerna..." | tee /dev/fd/3
SORTMERNA_DB=$(awk '/^SORTMERNA_DB/{print $3}' "${CONFIG}")
echo "Database directory for sortmerna = $SORTMERNA_DB" | tee /dev/fd/3

NAME=($(ls "$SORTMERNA_DB"/*.clustered.fasta* | sed 's/.fasta//')) 
#echo "name = " | tee /dev/fd/3
#printf '%s\n' "${NAME[@]}" | tee /dev/fd/3

SORTME_NAME=$(basename ""${NAME[@]}"")  

THREAD=$(awk '/^THREAD/{print $3}' "${CONFIG}")
#echo "Number of threads used = $THREAD" | tee /dev/fd/3

echo "Merging paired files into single files... " | tee /dev/fd/3
#bash merge-paired-reads.sh "$OUTPUT"/quality_control/"$SHORTNAME"_1_trimmed.fastq "$OUTPUT"/quality_control/"$SHORTNAME"_2_trimmed.fastq "$OUTPUT"/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq
reformat.sh in1="$OUTPUT"/quality_control/"$SHORTNAME"_1.fastq in2="$OUTPUT"/quality_control/"$SHORTNAME"_2.fastq out="$OUTPUT"/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq overwrite=t

echo "Filtering 16S/18S reads...." | tee /dev/fd/3
sortmerna --ref "$SORTMERNA_DB"/"$SORTME_NAME".fasta,"$SORTMERNA_DB"/$SORTME_NAME \
	--reads "$OUTPUT"/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq \
	--fastx \
	--aligned "$OUTPUT"/output_sortmerna/"$SHORTNAME"_16S18S \
	--paired_in \
	-a "$THREAD" \
	--log \
	-v

echo "Unmerging single files into paired files...." | tee /dev/fd/3
#bash unmerge-paired-reads.sh "$OUTPUT"/output_sortmerna/"$SHORTNAME"_16S18S.fastq "$OUTPUT"/output_sortmerna/"$SHORTNAME"_R1_16S18Sreads.fastq "$OUTPUT"/output_sortmerna/"$SHORTNAME"_R2_16S18Sreads.fastq
reformat.sh in="$OUTPUT"/output_sortmerna/"$SHORTNAME"_16S18S.fastq out1="$OUTPUT"/output_sortmerna/"$SHORTNAME"_R1_16S18Sreads.fastq out2="$OUTPUT"/output_sortmerna/"$SHORTNAME"_R2_16S18Sreads.fastq overwrite=t

echo "Saving results..." | tee /dev/fd/3

echo "Filtering 16S/18S using sortmerna ends successfully on : "`date` | tee /dev/fd/3

rm "$OUTPUT"/output_sortmerna/"$SHORTNAME"_16S18S.fastq

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#			Reconstructing 16S/18S full length sequences using EMIRGE
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo "Reconstructing 16S18S using EMIRGE starting on : "`date` | tee /dev/fd/3
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo ""

#echo "Setting up directories for EMIRGE..." | tee /dev/fd/3
mkdir -p "$OUTPUT/SSU_sequences/output_emirge"
rm -rf "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/*

#echo "Setting up indexed database for EMIRGE..." | tee /dev/fd/3
EMIRGE_DB=$(awk '/^EMIRGE_DB/{print $3}' "${CONFIG}")
echo "Database directory for EMIRGE = $EMIRGE_DB" | tee /dev/fd/3

NAME=($(ls "$EMIRGE_DB"/*.fasta*)) 
#echo "name = " | tee /dev/fd/3
#printf '%s\n' "${NAME[@]}" | tee /dev/fd/3

REF_NAME=$(basename ""${NAME[@]}"")  

NAME=($(ls "$EMIRGE_DB"/*.1.ebwt* | sed 's/.1.ebwt//')) 
#echo "name = " | tee /dev/fd/3
#printf '%s\n' "${NAME[@]}" | tee /dev/fd/3

BWT_NAME=$(basename ""${NAME[@]}"")  

MAX_LENGTH=$(awk '/^MAX_LENGTH/{print $3}' "${CONFIG}")
IDENTITY=$(awk '/^IDENTITY/{print $3}' "${CONFIG}")
NUM_ITERATION=$(awk '/^NUM_ITERATION/{print $3}' "${CONFIG}")
MEAN_INSERT_SIZE=$(awk '/^MEAN_INSERT_SIZE/{print $3}' "${CONFIG}")
STD_DEV=$(awk '/^STD_DEV/{print $3}' "${CONFIG}")


echo "Running emirge amplicon to reconstruct 16S/18S full length sequences..." | tee /dev/fd/3
emirge_amplicon.py \
	-1 "$OUTPUT"/output_sortmerna/"$SHORTNAME"_R1_16S18Sreads.fastq \
	-f "$EMIRGE_DB"/"$REF_NAME" \
	-b "$EMIRGE_DB"/$BWT_NAME \
	-l "$MAX_LENGTH" \
	-2 "$OUTPUT"/output_sortmerna/"$SHORTNAME"_R2_16S18Sreads.fastq \
	-j "$IDENTITY" \
	-n "$NUM_ITERATION" \
	-i "$MEAN_INSERT_SIZE" \
	-s "$STD_DEV" \
	-a "$THREAD" \
	--phred33 "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons

#-j 1 : 100 % identity
# -i : mean insert size of 300 bp
#-s : standard deviation of 100 bp 
# -l : max. read length of 300 bp (careful, if not given correctly, command will stop!!!!!)
# -n : number of iterations
#-a : number of cores

echo "Merging iterations into fasta sequence... " | tee /dev/fd/3

emirge_rename_fasta.py --no_N "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/iter.$NUM_ITERATION > "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_renamed_16S18S_recons.fasta

#rm -rv "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/!("initial_mapping"|"iter.$NUM_ITERATION")

#zip -jr  "$SHORTNAME"_amplicon_16S18S_recons.zip "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/*

cd "$OUTPUT"/SSU_sequences/output_emirge && zip -rm "$SHORTNAME"_amplicon_16S18S_recons.zip "$SHORTNAME"_amplicon_16S18S_recons/ && cd -

#rm -r "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons


echo "Finalising reconstructed sequences..." | tee /dev/fd/3

cat "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_renamed_16S18S_recons.fasta "$OUTPUT"/SSU_sequences/output_MetaRib/"$SHORTNAME"_MetaRib_SSU.fasta > "$OUTPUT"/SSU_sequences/emirge_metarib_SSU_sequences.fasta

#clustering at 97%
vsearch --cluster_fast "$OUTPUT"/SSU_sequences/emirge_metarib_SSU_sequences.fasta --centroids "$OUTPUT"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.fasta --id 0.97

#covert small letters into capital letters
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' "$OUTPUT"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.fasta > "$OUTPUT"/SSU_sequences/all_SSU_sequences.fasta

#rm -rv "$OUTPUT"/output_MetaRib/"$SHORTNAME"/Iteration/iter_1/emirge_amp/!("initial_mapping"|"iter.$NUM_ITERATION")

#cd "$OUTPUT"/SSU_sequences/output_MetaRib/"$SHORTNAME" && zip -rm Iteration.zip Iteration/ && cd -

#zip -r Iteration.zip "$OUTPUT"/output_MetaRib/"$SHORTNAME"/Iteration/

#rm -r "$OUTPUT"/output_MetaRib/"$SHORTNAME"/Iteration

#mv Iteration.zip "$OUTPUT"/output_MetaRib/"$SHORTNAME"

echo "Saving results..." | tee /dev/fd/3

echo "Reconstructing 16S/18S sequences ends successfully on : "`date` | tee /dev/fd/3

echo "Calculating relative abundances of reconstructed sequences..." | tee /dev/fd/3

bbmap.sh -Xmx3g in="$OUTPUT"/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq \
	ref="$OUTPUT"/SSU_sequences/all_SSU_sequences.fasta \
	vslow \
	k=12 \
	minid=0.96 \
	maxindel=1 \
	minhits=2 \
	idfilter=0.98 \
	scafstats="$OUTPUT"/SSU_sequences/all_scafstats.txt \
	covstats="$OUTPUT"/SSU_sequences/all_covstats.txt

cat "$OUTPUT"/SSU_sequences/all_scafstats.txt | sed 1d | tr ',' \\t | awk '{ print $1,$8 }' | sort -k1 -k2 | awk '!($2==0){print}' | awk '{print $1}' > "$OUTPUT"/SSU_sequences/id_file.txt

cat "$OUTPUT"/SSU_sequences/all_scafstats.txt | sed 1d | tr ',' \\t | awk '{ print $1,$8 }' | sort -k1 -k2 | awk '!($2==0){print}' > "$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats.txt

awk 'NR==FNR{ids[$0]; next} ($1 in ids){ printf ">" $0 }' "$OUTPUT"/SSU_sequences/id_file.txt RS='>' "$OUTPUT"/SSU_sequences/all_SSU_sequences.fasta > "$OUTPUT"/SSU_sequences/"$SHORTNAME"_SSU_sequences.fasta

awk -F '\t' 'NR==FNR {id[$1]; next} $1 in id' "$OUTPUT"/SSU_sequences/id_file.txt "$OUTPUT"/SSU_sequences/all_covstats.txt > "$OUTPUT"/SSU_sequences/"$SHORTNAME"_covstats.txt

rm "$OUTPUT"/SSU_sequences/all_scafstats.txt
rm "$OUTPUT"/SSU_sequences/all_covstats.txt
rm "$OUTPUT"/SSU_sequences/id_file.txt
rm "$OUTPUT"/SSU_sequences/all_SSU_sequences.fasta

rm "$OUTPUT"/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq
rm "$OUTPUT"/SSU_sequences/emirge_metarib_SSU_sequences.fasta
rm "$OUTPUT"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.fasta



echo "Saving results..." | tee /dev/fd/3

echo "Relative abundance calculation ends successfully on : "`date` | tee /dev/fd/3

conda deactivate
#echo "RiboTaxa virtual environment has been deactivated successfully..." | tee /dev/fd/3
