#!/bin/sh

exec 3>&1 1>RiboTaxa.log 2>&3

# Program configuration
__author__='Oshma Chakoory'
__email__='oshma.chakoory@uca.fr'
__credits__=["Oshma"]
__status__='Development'
__version__='1.4'

echo " "
echo "RiboTaxa -- A complete pipeline from raw metagenomics to species-level identification" | tee /dev/fd/3
echo "By Oshma Chakoory, Sophie Marre & Pierre Peyret" | tee /dev/fd/3
echo "University Clermont Auvergne, France " | tee /dev/fd/3
echo "Version: 1.4"

echo "This program is distributed under the AGPL-3.0 License. See LICENSE for more information."


# Handling errors
#set -x # debug mode on
set -o errexit # ensure script will stop in case of ignored error

CONFIG_PATH=$1
CONFIG="${CONFIG_PATH[@]}"

DATA_DIR=$(awk '/^DATA_DIR/{print $3}' "${CONFIG}")
#echo "Data path = $DATA_DIR" | tee /dev/fd/3

OUTPUT=$(awk '/^OUTPUT/{print $3}' "${CONFIG}")

RiboTaxa_DIR=$(awk '/^RiboTaxa_DIR/{print $3}' "${CONFIG}")

FORMAT=$(awk '/^FORMAT/{print $3}' "${CONFIG}")
#echo $FORMAT

THREAD=$(awk '/^THREAD/{print $3}' "${CONFIG}")
#echo "Number of threads used = $THREAD" | tee /dev/fd/3

#Set up results sub-directories
mkdir -p "$OUTPUT/quality_control"
mkdir -p "$OUTPUT/quality_control/before_fastqc"
mkdir -p "$OUTPUT/quality_control/before_fastqc/multiqc"
mkdir -p "$OUTPUT/quality_control/after_fastqc"
mkdir -p "$OUTPUT/quality_control/after_fastqc/multiqc"

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

#activate virtual environment for RiboTaxa
conda activate RiboTaxa_py27

for NAME in `ls "$DATA_DIR"/*_1.$FORMAT | sed 's/_1.'$FORMAT'//'` 
do
SHORTNAME=$(basename ""${NAME[@]}"") 
echo $SHORTNAME

echo "Run the First quality control on raw data... " | tee /dev/fd/3

fastqc "$DATA_DIR"/"$SHORTNAME"_1.$FORMAT  "$DATA_DIR"/"$SHORTNAME"_2.$FORMAT  -dir $OUTPUT -o "$OUTPUT"/quality_control/before_fastqc

echo "Removing adapters from sequences..." | tee /dev/fd/3
bbduk.sh -Xmx1g \
	in1="$DATA_DIR"/"$SHORTNAME"_1.$FORMAT  \
	in2="$DATA_DIR"/"$SHORTNAME"_2.$FORMAT  \
	out1="$OUTPUT"/quality_control/"$SHORTNAME"_1_noadapt.$FORMAT  \
	out2="$OUTPUT"/quality_control/"$SHORTNAME"_2_noadapt.$FORMAT  \
	ref="$RiboTaxa_DIR"/adapters/TruSeq3-PE.fa \
	ktrim=$KTRIM \
	k=$KMER \
	mink=11 \
	tpe \
	tbo

echo "Trimming sequences..." | tee /dev/fd/3
bbduk.sh -Xmx2g \
	in1="$OUTPUT"/quality_control/"$SHORTNAME"_1_noadapt.$FORMAT  \
	in2="$OUTPUT"/quality_control/"$SHORTNAME"_2_noadapt.$FORMAT  \
	out1="$OUTPUT"/quality_control/"$SHORTNAME"_1.$FORMAT  \
	out2="$OUTPUT"/quality_control/"$SHORTNAME"_2.$FORMAT  \
	minlen=$MINLENGTH \
	qtrim=$QTRIM \
	trimq=$TRIMQ \
	maxns=$MAXNS


echo "Run the second quality control on trimmed data..." | tee /dev/fd/3
fastqc "$OUTPUT"/quality_control/"$SHORTNAME"_1.$FORMAT  "$OUTPUT"/quality_control/"$SHORTNAME"_2.$FORMAT  -dir $OUTPUT -o "$OUTPUT"/quality_control/after_fastqc

echo "Running multiQC..." | tee /dev/fd/3

multiqc -f "$OUTPUT"/quality_control/before_fastqc/* --outdir "$OUTPUT"/quality_control/before_fastqc/multiqc/

multiqc -f "$OUTPUT"/quality_control/after_fastqc/* --outdir "$OUTPUT"/quality_control/after_fastqc/multiqc/

echo "Saving results for quality control..." | tee /dev/fd/3

echo "Quality control ends successfully on : "`date` | tee /dev/fd/3
 
rm "$OUTPUT"/quality_control/*_noadapt.$FORMAT 
done

for NAME in `ls "$OUTPUT"/quality_control/*fastq*`; do
	if [[ ${NAME[@]} == *.gz* ]]; then
		gunzip -f $NAME
#	else
#		echo "files are unzipped"
	fi
done


for NAME in `ls "$OUTPUT"/quality_control/*_1.fastq |sed 's/_1.fastq//'` 
do
FILE=$(basename $NAME)
echo $FILE >> "$OUTPUT"/quality_control/samples.list.txt
done

cat "$OUTPUT"/quality_control/*_1.fastq > "$OUTPUT"/quality_control/all_1.fastq 
cat "$OUTPUT"/quality_control/*_2.fastq > "$OUTPUT"/quality_control/all_2.fastq

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo "Running MetaRib on raw reads to reconstruct 16S/18S full length sequences..."`date` | tee /dev/fd/3
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3

#SAMPLE=$(awk '{s++}END{print s/4}' "$OUTPUT"/output_sortmerna/"$SHORTNAME"_R1_16S18Sreads.fastq)

mkdir -p "$OUTPUT/SSU_sequences"

#echo "Setting up indexed database for EMIRGE..." | tee /dev/fd/3
EMIRGE_DB=$(awk '/^EMIRGE_DB/{print $3}' "${CONFIG}")
#echo "Database directory for EMIRGE = $EMIRGE_DB" | tee /dev/fd/3
NAME=($(ls "$EMIRGE_DB"/*.fasta*)) 
#echo "name = " | tee /dev/fd/3
#printf '%s\n' "${NAME[@]}" | tee /dev/fd/3
REF_NAME=$(basename ""${NAME[@]}"")  
NAME=($(ls "$EMIRGE_DB"/*.1.ebwt* | sed 's/.1.ebwt//')) 
#echo "name = " | tee /dev/fd/3
#printf '%s\n' "${NAME[@]}" | tee /dev/fd/3
BWT_NAME=$(basename ""${NAME[@]}"")  

python2 "$RiboTaxa_DIR"/run_MetaRib.py -cfg "$CONFIG_PATH" -p "$OUTPUT"/quality_control -b "$EMIRGE_DB"/$BWT_NAME -l "$EMIRGE_DB"/"$REF_NAME" >> RiboTaxa.log


for NAME in `cat "$OUTPUT"/quality_control/samples.list.txt`
do
echo $NAME
cat "$OUTPUT"/SSU_sequences/output_MetaRib/Abundance/all.dedup.filtered.est.ab.txt | tr -s '\t' ',' | csvcut -c Contig_ID,"$NAME"_estab| tr ',' '\t' | sed 1d | awk '!($2=="0.00000"){print}' |awk '{print $1}' > "$OUTPUT"/SSU_sequences/output_MetaRib/file.txt
awk 'NR==FNR{ids[$0]; next} ($1 in ids){ printf ">" $0 }' "$OUTPUT"/SSU_sequences/output_MetaRib/file.txt RS='>' "$OUTPUT"/SSU_sequences/output_MetaRib/Abundance/all.dedup.filtered.fasta > "$OUTPUT"/SSU_sequences/output_MetaRib/"$NAME"_MetaRib_SSU.fasta
rm "$OUTPUT"/SSU_sequences/output_MetaRib/file.txt
done

conda deactivate

rm "$OUTPUT"/quality_control/all_1.fastq
rm "$OUTPUT"/quality_control/all_2.fastq
rm "$OUTPUT"/quality_control/samples.list.txt
rm -r "$OUTPUT"/SSU_sequences/output_MetaRib/Abundance

for NAME in `ls "$OUTPUT"/quality_control/*_1.fastq | sed 's/_1.fastq//'` 
do
SHORTNAME=$(basename ""${NAME[@]}"") 
#run RiboTaxa.sh script to perform
### Perform a quality control on the (meta)genomics data using BBTOOLS
### Filter 16S/18S reads using SORTMERNA
### Reconstruct full length 16S/18S rRNA sequences using EMIRGE

source "$RiboTaxa_DIR"/RiboTaxa.sh $CONFIG_PATH 

#run sklearn_classfier.sh script to perform
### Taxonomic classification of full length reconstrcuted sequences using sklearn classifier of Qiime2

source "$RiboTaxa_DIR"/sklearn_classifier.sh $CONFIG_PATH 

done

mv RiboTaxa.log "$OUTPUT"
#mv "$OUTPUT"/output_MetaRib "$OUTPUT"/SSU_sequences