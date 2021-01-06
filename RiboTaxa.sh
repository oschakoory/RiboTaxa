#!/bin/sh

# Program configuration
__author__='Oshma Chakoory'
__email__='oshma.chakoory@uca.fr'
__credits__=["Oshma"]
__status__='Development'
__version__='0.0.1'


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
mkdir -p "$OUTPUT"

#Set up results sub-directories
mkdir -p "$OUTPUT/quality_control"
mkdir -p "$OUTPUT/quality_control/before_fastqc"
mkdir -p "$OUTPUT/quality_control/before_fastqc/multiqc"
mkdir -p "$OUTPUT/quality_control/after_fastqc"
mkdir -p "$OUTPUT/quality_control/after_fastqc/multiqc"

NAME=($(ls "$DATA_DIR"/*"$SHORTNAME"_1*))
FILE=$(basename ""${NAME[@]}"") 
#echo "FILE=$FILE"

#echo "Verifying the format of input files..." | tee /dev/fd/3
for NAME in `ls "$DATA_DIR"/*"$SHORTNAME"*`; do
	FILE=($(basename ""${NAME[@]}""))
	if [[ ${NAME[@]} == *.gz* ]]; then
		echo "$FILE is gzipped"
		gunzip ${NAME[@]}
	else
		echo "$FILE is not gzipped"
	fi
done

echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo "Quality control starting on : "`date` | tee /dev/fd/3
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo ""

KTRIM=$(awk '/^ktrim/{print $3}' "${CONFIG}")

KMER=$(awk '/^kmer/{print $3}' "${CONFIG}")

MINLENGTH=$(awk '/^minlength/{print $3}' "${CONFIG}")

TRIMQ=$(awk '/^trimq/{print $3}' "${CONFIG}")

QTRIM=$(awk '/^qtrim/{print $3}' "${CONFIG}")

MAXNS=$(awk '/^maxns/{print $3}' "${CONFIG}")

echo "Run the First quality control on raw data... " | tee /dev/fd/3

fastqc "$DATA_DIR"/"$SHORTNAME"_1.fastq "$DATA_DIR"/"$SHORTNAME"_2.fastq -dir $OUTPUT -o "$OUTPUT"/quality_control/before_fastqc

echo "Removing adapters from sequences..." | tee /dev/fd/3
bbduk.sh -Xmx1g in1="$DATA_DIR"/"$SHORTNAME"_1.fastq in2="$DATA_DIR"/"$SHORTNAME"_2.fastq out1="$OUTPUT"/quality_control/"$SHORTNAME"_1_noadapt.fastq out2="$OUTPUT"/quality_control/"$SHORTNAME"_2_noadapt.fastq ref="$RiboTaxa_DIR"/adapters/TruSeq3-PE.fa ktrim=$KTRIM k=$KMER mink=11 tpe tbo

echo "Trimming sequences..." | tee /dev/fd/3
bbduk.sh -Xmx2g \
	in1="$OUTPUT"/quality_control/"$SHORTNAME"_1_noadapt.fastq \
	in2="$OUTPUT"/quality_control/"$SHORTNAME"_2_noadapt.fastq \
	out1="$OUTPUT"/quality_control/"$SHORTNAME"_1_trimmed.fastq \
	out2="$OUTPUT"/quality_control/"$SHORTNAME"_2_trimmed.fastq \
	minlen=$MINLENGTH \
	qtrim=$QTRIM \
	trimq=$TRIMQ \
	maxns=$MAXNS


echo "Run the second quality control on trimmed data..." | tee /dev/fd/3
fastqc "$OUTPUT"/quality_control/"$SHORTNAME"_1_trimmed.fastq "$OUTPUT"/quality_control/"$SHORTNAME"_2_trimmed.fastq -dir $OUTPUT -o "$OUTPUT"/quality_control/after_fastqc

echo "Running multiQC..." | tee /dev/fd/3

multiqc -f "$OUTPUT"/quality_control/before_fastqc/* --outdir "$OUTPUT"/quality_control/before_fastqc/multiqc/

multiqc -f "$OUTPUT"/quality_control/after_fastqc/* --outdir "$OUTPUT"/quality_control/after_fastqc/multiqc/

echo "Saving results for quality control..." | tee /dev/fd/3

echo "Quality control ends successfully on : "`date` | tee /dev/fd/3

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
bash merge-paired-reads.sh "$OUTPUT"/quality_control/"$SHORTNAME"_1_trimmed.fastq "$OUTPUT"/quality_control/"$SHORTNAME"_2_trimmed.fastq "$OUTPUT"/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq

echo "Filtering 16S/18S reads...." | tee /dev/fd/3
sortmerna --ref "$SORTMERNA_DB"/"$SORTME_NAME".fasta,"$SORTMERNA_DB"/$SORTME_NAME \
	--reads "$OUTPUT"/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq \
	--fastx \
	--aligned "$OUTPUT"/output_sortmerna/"$SHORTNAME"_16S18S \
	--other "$OUTPUT"/output_sortmerna/"$SHORTNAME"_other_than_16S18S \
	--paired_in \
	-a "$THREAD" \
	--log \
	-v


echo "Unmerging single files into paired files...." | tee /dev/fd/3
bash unmerge-paired-reads.sh "$OUTPUT"/output_sortmerna/"$SHORTNAME"_16S18S.fastq "$OUTPUT"/output_sortmerna/"$SHORTNAME"_R1_16S18Sreads.fastq "$OUTPUT"/output_sortmerna/"$SHORTNAME"_R2_16S18Sreads.fastq

echo "Saving results..." | tee /dev/fd/3

echo "Filtering 16S/18S using sortmerna ends successfully on : "`date` | tee /dev/fd/3

rm "$OUTPUT"/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq
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
mkdir -p "$OUTPUT/output_emirge"
rm -rf "$OUTPUT"/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/*

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


echo "Running emirge amplicon to reconstrct 16S/18S full length sequences..." | tee /dev/fd/3
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
	--phred33 "$OUTPUT"/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons

#-j 1 : 100 % identity
# -i : mean insert size of 300 bp
#-s : standard deviation of 100 bp 
# -l : max. read length of 300 bp (careful, if not given correctly, command will stop!!!!!)
# -n : number of iterations
#-a : number of cores

echo "Merging iterations into fasta sequence... " | tee /dev/fd/3

emirge_rename_fasta.py "$OUTPUT"/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/iter.$NUM_ITERATION > "$OUTPUT"/output_emirge/"$SHORTNAME"_renamed_16S18S_recons.fasta

echo "Saving results..." | tee /dev/fd/3

echo "Reconstructing 16S18S using EMIRGE ends successfully on : "`date` | tee /dev/fd/3

conda deactivate
#echo "RiboTaxa virtual environment has been deactivated successfully..." | tee /dev/fd/3
