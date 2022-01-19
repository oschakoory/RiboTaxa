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

#activate virtual environment for RiboTaxa
conda activate RiboTaxa_py27
#echo "RiboTaxa virtual environment has been activated successfully..." | tee /dev/fd/3

#echo "Setting up directories for quality control..." | tee /dev/fd/3

#set up 
CONFIG_PATH=$1
CONFIG="${CONFIG_PATH[@]}"

SHORTNAME=$2
#SHORTNAME=$(basename ""${NAME[@]}"" | sed 's/.'$FORMAT'//') 
echo "SHORTNAME = " >&2
printf '%s\n' "$SHORTNAME" >&2

RiboTaxa_DIR=$(awk '/^RiboTaxa_DIR/{print $3}' "${CONFIG}")

#Set up data directory
DATA_DIR=$(awk '/^DATA_DIR/{print $3}' "${CONFIG}")
#echo "Data path = $DATA_DIR" | tee /dev/fd/3

#Set up output directory
OUTPUT=$(awk '/^OUTPUT/{print $3}' "${CONFIG}")
#rm -rf "$OUTPUT"
mkdir -p "$OUTPUT"

RESULTS="$OUTPUT/$SHORTNAME"
mkdir -p $RESULTS

#fomrat of files (can either be fastq or fastq.gz)
FORMAT=$(awk '/^FORMAT/{print $3}' "${CONFIG}")

#Set up results sub-directories
mkdir -p "$RESULTS/quality_control"
mkdir -p "$RESULTS/quality_control/before_fastqc"
mkdir -p "$OUTPUT/multiqc"
mkdir -p "$RESULTS/quality_control/after_fastqc"
mkdir -p "$OUTPUT/multiqc/before_qc"
mkdir -p "$OUTPUT/multiqc/after_qc"

echo "" | tee /dev/fd/3
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo ">Quality control starting on : "`date` | tee /dev/fd/3
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo "" | tee /dev/fd/3

RAM=$(awk '/^RAM/{print $3}' "${CONFIG}")

KTRIM=$(awk '/^ktrim/{print $3}' "${CONFIG}")

KMER=$(awk '/^kmer/{print $3}' "${CONFIG}")

MINLENGTH=$(awk '/^minlength/{print $3}' "${CONFIG}")

TRIMQ=$(awk '/^trimq/{print $3}' "${CONFIG}")

QTRIM=$(awk '/^qtrim/{print $3}' "${CONFIG}")

MAXNS=$(awk '/^maxns/{print $3}' "${CONFIG}")

THREAD=$(awk '/^THREAD/{print $3}' "${CONFIG}")
#echo "Number of threads used = $THREAD" | tee /dev/fd/3

FORWARD="$RESULTS/quality_control/"$SHORTNAME"_trimmed.fastq"

if [ -f "$FORWARD" ];
then
   echo "Forward file : '${FORWARD}' is present."
else
#echo "File '${FORWARD}' not found."

echo "Run the First quality control on raw data... " | tee /dev/fd/3

echo "command line :" | tee /dev/fd/3
echo "fastqc $DATA_DIR/"$SHORTNAME".$FORMAT -dir $RESULTS -o $RESULTS/quality_control/before_fastqc" | tee /dev/fd/3

fastqc "$DATA_DIR"/"$SHORTNAME".$FORMAT -dir $RESULTS -o "$RESULTS"/quality_control/before_fastqc

echo "Removing adapters from sequences..." | tee /dev/fd/3
echo "command line :" | tee /dev/fd/3
echo "bbduk.sh -Xmx${RAM}g \
	in=$DATA_DIR/"$SHORTNAME".$FORMAT  \
	out=$RESULTS/quality_control/"$SHORTNAME"_noadapt.$FORMAT  \
	ref=$RiboTaxa_DIR/adapters/TruSeq3-SE.fa \
	ktrim=$KTRIM \
	k=$KMER \
	mink=11 \
	threads=$THREAD \
	tpe \
	tbo" | tee /dev/fd/3

bbduk.sh -Xmx${RAM}g \
	in="$DATA_DIR"/"$SHORTNAME".$FORMAT  \
	out="$RESULTS"/quality_control/"$SHORTNAME"_noadapt.$FORMAT  \
	ref="$RiboTaxa_DIR"/adapters/TruSeq3-SE.fa \
	ktrim=$KTRIM \
	k=$KMER \
	mink=11 \
	threads=$THREAD \
	tpe \
	tbo


echo "Trimming sequences..." | tee /dev/fd/3
echo "command line :" | tee /dev/fd/3
echo "bbduk.sh -Xmx${RAM}g \
	in=$RESULTS/quality_control/"$SHORTNAME"_noadapt.$FORMAT  \
	out=$RESULTS/quality_control/"$SHORTNAME"_trimmed.$FORMAT  \
	minlen=$MINLENGTH \
	qtrim=$QTRIM \
	trimq=$TRIMQ \
	threads=$THREAD \
	maxns=$MAXNS" | tee /dev/fd/3

bbduk.sh -Xmx${RAM}g \
	in="$RESULTS"/quality_control/"$SHORTNAME"_noadapt.$FORMAT  \
	out="$RESULTS"/quality_control/"$SHORTNAME"_trimmed.$FORMAT  \
	minlen=$MINLENGTH \
	qtrim=$QTRIM \
	trimq=$TRIMQ \
	threads=$THREAD \
	maxns=$MAXNS


echo "Run the second quality control on trimmed data..." | tee /dev/fd/3

echo "command line :" | tee /dev/fd/3
echo "fastqc $RESULTS/quality_control/"$SHORTNAME"_trimmed.$FORMAT -dir $RESULTS -o $RESULTS/quality_control/after_fastqc" | tee /dev/fd/3

fastqc "$RESULTS"/quality_control/"$SHORTNAME"_trimmed.$FORMAT  -dir $RESULTS -o "$RESULTS"/quality_control/after_fastqc

echo "Running multiQC..." | tee /dev/fd/3

multiqc -f "$OUTPUT"/*/quality_control/before_fastqc/* --outdir "$OUTPUT/multiqc/before_qc"

multiqc -f "$OUTPUT"/*/quality_control/after_fastqc/* --outdir "$OUTPUT/multiqc/after_qc"

echo "Saving results for quality control..." | tee /dev/fd/3

echo ">Quality control ends successfully on : "`date` | tee /dev/fd/3
 
rm "$RESULTS"/quality_control/*_noadapt.$FORMAT 

fi

#unzip fastq files if compressed as sortmerna takes uncompressed files only
for IDNAME in `ls "$RESULTS"/quality_control/*"$SHORTNAME"*`; do
	FILE=($(basename ""${IDNAME[@]}""))
	#echo $FILE
	if [[ ${IDNAME[@]} == *.gz* ]]; then
		#echo "$FILE is gzipped"
		gunzip -f $IDNAME
#	else
#		echo "$FILE is not gzipped"
	fi
done

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#				Filtering 16S/18S reads using sortemrna
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#echo "Setting up directories for sortmerna..." | tee /dev/fd/3
echo ""
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo ">Filtering 16S/18S using sortmerna starting on : "`date` | tee /dev/fd/3
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo ""

mkdir -p "$RESULTS/output_sortmerna"

#echo "Setting up indexed database for sortmerna..." | tee /dev/fd/3
SORTMERNA_DB=$(awk '/^SORTMERNA_DB/{print $3}' "${CONFIG}")
echo "Database directory for sortmerna = $SORTMERNA_DB" | tee /dev/fd/3

SMRNAME=($(ls "$SORTMERNA_DB"/*.clustered.fasta* | sed 's/.fasta//')) 
#echo "name = " | tee /dev/fd/3
#printf '%s\n' "${NAME[@]}" | tee /dev/fd/3

SORTME_NAME=$(basename ""${SMRNAME[@]}"")  


FORWARD="$RESULTS/output_sortmerna/"$SHORTNAME"_16S18S.fastq"

if [ -f "$FORWARD" ];
then
	echo "Forward file : '${FORWARD}' is present."
else
#echo "Merging paired files into single files... " | tee /dev/fd/3
#echo "command line :" | tee /dev/fd/3
#echo "reformat.sh in1=$RESULTS/quality_control/"$SHORTNAME"_1_trimmed.fastq in2=$RESULTS/quality_control/"$SHORTNAME"_2_trimmed.fastq out=$RESULTS/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq overwrite=t" | tee /dev/fd/3

#bash merge-paired-reads.sh "$RESULTS"/quality_control/"$SHORTNAME"_1_trimmed.fastq "$RESULTS"/quality_control/"$SHORTNAME"_2_trimmed.fastq "$RESULTS"/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq
##reformat.sh in1="$RESULTS"/quality_control/"$SHORTNAME"_1_trimmed.fastq in2="$RESULTS"/quality_control/"$SHORTNAME"_2_trimmed.fastq out="$RESULTS"/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq overwrite=t


echo "Filtering 16S/18S reads...." | tee /dev/fd/3
echo "command line :" | tee /dev/fd/3
echo "sortmerna --ref $SORTMERNA_DB/$SORTME_NAME.fasta,$SORTMERNA_DB/$SORTME_NAME \
	--reads "$RESULTS"/quality_control/"$SHORTNAME"_trimmed.fastq \
	--fastx \
	--aligned $RESULTS/output_sortmerna/"$SHORTNAME"_16S18S \
	-a $THREAD \
	--log \
	-v" | tee /dev/fd/3

sortmerna --ref $SORTMERNA_DB/$SORTME_NAME.fasta,$SORTMERNA_DB/$SORTME_NAME \
	--reads $RESULTS/quality_control/"$SHORTNAME"_trimmed.fastq \
	--fastx \
	--aligned $RESULTS/output_sortmerna/"$SHORTNAME"_16S18S \
	-a $THREAD \
	--log \
	-v

#echo "Unmerging single files into paired files...." | tee /dev/fd/3
#bash unmerge-paired-reads.sh "$RESULTS"/output_sortmerna/"$SHORTNAME"_16S18S.fastq "$RESULTS"/output_sortmerna/"$SHORTNAME"_R1_16S18Sreads.fastq "$RESULTS"/output_sortmerna/"$SHORTNAME"_R2_16S18Sreads.fastq
#echo "command line :" | tee /dev/fd/3
#echo "reformat.sh in=$RESULTS/output_sortmerna/"$SHORTNAME"_16S18S.fastq out1=$RESULTS/output_sortmerna/"$SHORTNAME"_R1_16S18Sreads.fastq out2=$RESULTS/output_sortmerna/"$SHORTNAME"_R2_16S18Sreads.fastq overwrite=t" | tee /dev/fd/3

#reformat.sh in1="$RESULTS"/output_sortmerna/"$SHORTNAME"_R1_16S18Sreads.fastq in2="$RESULTS"/output_sortmerna/"$SHORTNAME"_R2_16S18Sreads.fastq overwrite=t out="$RESULTS"/output_sortmerna/"$SHORTNAME"_16S18S.fastq 

#reformat.sh in="$RESULTS"/output_sortmerna/"$SHORTNAME"_16S18S.fastq out1="$RESULTS"/output_sortmerna/"$SHORTNAME"_R1_16S18Sreads.fastq out2="$RESULTS"/output_sortmerna/"$SHORTNAME"_R2_16S18Sreads.fastq overwrite=t reads=31013

echo "Saving results..." | tee /dev/fd/3

#rm "$RESULTS"/output_sortmerna/"$SHORTNAME"_16S18S.fastq
#rm -f "$RESULTS"/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq

fi

echo ">Filtering 16S/18S using sortmerna ends successfully on : "`date` | tee /dev/fd/3


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#			Reconstructing 16S/18S full length sequences using EMIRGE
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

echo "" | tee /dev/fd/3
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo ">Reconstructing 16S18S using EMIRGE starting on : "`date` | tee /dev/fd/3
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo "" | tee /dev/fd/3

#echo "Setting up directories for EMIRGE..." | tee /dev/fd/3
mkdir -p "$RESULTS/SSU_sequences/output_emirge"

#echo "Setting up indexed database for EMIRGE..." | tee /dev/fd/3
EMIRGE_DB=$(awk '/^EMIRGE_DB/{print $3}' "${CONFIG}")
echo "Database directory for EMIRGE = $EMIRGE_DB" | tee /dev/fd/3

NAMEDB=($(ls "$EMIRGE_DB"/*.fasta*)) 
#echo "name = " | tee /dev/fd/3
#printf '%s\n' "${DBNAME[@]}" | tee /dev/fd/3

REF_NAME=$(basename ""${NAMEDB[@]}"")  

DBNAME=($(ls "$EMIRGE_DB"/*.4.ebwt* | sed 's/.4.ebwt//')) 
#echo "name = " | tee /dev/fd/3
#printf '%s\n' "${DBNAME[@]}" | tee /dev/fd/3

BWT_NAME=$(basename ""${DBNAME[@]}"") 

MAX_LENGTH=$(awk '/^MAX_LENGTH/{print $3}' "${CONFIG}")
IDENTITY=$(awk '/^IDENTITY/{print $3}' "${CONFIG}")
NUM_ITERATION=$(awk '/^NUM_ITERATION/{print $3}' "${CONFIG}")
MEAN_INSERT_SIZE=$(awk '/^MEAN_INSERT_SIZE/{print $3}' "${CONFIG}")
STD_DEV=$(awk '/^STD_DEV/{print $3}' "${CONFIG}")
MIN_COV=$(awk '/^MIN_COV/{print $3}' "${CONFIG}")


EMIRGE_FILE="$RESULTS/SSU_sequences/output_emirge/"$SHORTNAME"_renamed_16S18S_recons.fasta"

if [ -f "$EMIRGE_FILE" ];
then
	echo "Emirge file : '${EMIRGE_FILE}' is present."
else

rm -rf "$RESULTS"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/*

echo "command line :" | tee /dev/fd/3
echo "emirge_amplicon.py \
	-1 $RESULTS/output_sortmerna/"$SHORTNAME"_16S18S.fastq \
	-f $EMIRGE_DB/$REF_NAME \
	-b $EMIRGE_DB/$BWT_NAME \
	-l $MAX_LENGTH \
	-j $IDENTITY \
	-n $NUM_ITERATION \
	-i $MEAN_INSERT_SIZE \
	-s $STD_DEV \
	-a $THREAD \
	-c "$MIN_COV" \
	--phred33 $RESULTS/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons" | tee /dev/fd/3


emirge_amplicon.py \
	-1 "$RESULTS"/output_sortmerna/"$SHORTNAME"_16S18S.fastq \
	-f "$EMIRGE_DB"/"$REF_NAME" \
	-b "$EMIRGE_DB"/$BWT_NAME \
	-l "$MAX_LENGTH" \
	-j "$IDENTITY" \
	-n "$NUM_ITERATION" \
	-i "$MEAN_INSERT_SIZE" \
	-s "$STD_DEV" \
	-a "$THREAD" \
	-c "$MIN_COV" \
	--phred33 "$RESULTS"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons

#-j 1 : 100 % identity
# -i : mean insert size of 300 bp
#-s : standard deviation of 100 bp 
# -l : max. read length of 300 bp (careful, if not given correctly, command will stop!!!!!)
# -n : number of iterations
#-a : number of cores

echo "Merging iterations into fasta sequence... " | tee /dev/fd/3
echo "command line :" | tee /dev/fd/3
echo "emirge_rename_fasta.py --no_N $RESULTS/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/iter.$NUM_ITERATION > $RESULTS/SSU_sequences/output_emirge/"$SHORTNAME"_renamed_16S18S_recons.fasta" | tee /dev/fd/3

emirge_rename_fasta.py --no_N "$RESULTS"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/iter.$NUM_ITERATION > "$RESULTS"/SSU_sequences/output_emirge/"$SHORTNAME"_renamed_16S18S_recons.fasta

#rm -rv "$RESULTS"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/!("initial_mapping"|"iter.$NUM_ITERATION")

#zip -jr  "$SHORTNAME"_amplicon_16S18S_recons.zip "$RESULTS"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/*

#cd "$RESULTS"/SSU_sequences/output_emirge && zip -qrm "$SHORTNAME"_amplicon_16S18S_recons.zip "$SHORTNAME"_amplicon_16S18S_recons/ && cd -

#rm -r "$RESULTS"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons
if [ -d ""$RESULTS"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons" ]; 
	then 
		mv "$RESULTS"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/"iter.$NUM_ITERATION" ""$RESULTS"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons"/"last.$NUM_ITERATION"
		rm -rf ""$RESULTS"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons"/*iter*; 
		mv "$RESULTS"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/"last.$NUM_ITERATION" ""$RESULTS"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons"/"iter.$NUM_ITERATION"
fi


fi


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#			Reconstructing 16S/18S full length sequences using MetaRib
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MetaRib_FILE="$RESULTS/SSU_sequences/output_MetaRib/MetaRib_SSU.fasta"

if [ -f "$MetaRib_FILE" ];
then
	echo "MetaRib file : '${MetaRib_FILE}' is present."
else
echo ">Running MetaRib to reconstruct 16S/18S full length sequences..."`date` | tee /dev/fd/3

#SAMPLE=$(awk '{s++}END{print s/4}' "$RESULTS"/output_sortmerna/"$SHORTNAME"_R1_16S18Sreads.fastq)

mkdir -p  "$RESULTS"/SSU_sequences/output_MetaRib/"$SHORTNAME"

echo $SHORTNAME > ""$RESULTS"/SSU_sequences/output_MetaRib/samples.list.txt"

echo "command line :" | tee /dev/fd/3
echo "python2 $RiboTaxa_DIR/scripts/run_MetaRib_SE.py -cfg $CONFIG_PATH -p $RESULTS/quality_control -b $EMIRGE_DB/$BWT_NAME -l $EMIRGE_DB/$REF_NAME" -o "$RESULTS"/SSU_sequences/output_MetaRib | tee /dev/fd/3

python2 "$RiboTaxa_DIR"/scripts/run_MetaRib_SE.py -cfg "$CONFIG_PATH" -p "$RESULTS"/quality_control -b "$EMIRGE_DB"/$BWT_NAME -l "$EMIRGE_DB"/"$REF_NAME" -o "$RESULTS"/SSU_sequences/output_MetaRib

#cd "$RESULTS"/output_MetaRib && zip -qrm Iteration.zip Iteration/ && cd -

if [ -d ""$RESULTS"/SSU_sequences/output_MetaRib/Iteration/iter_1/emirge_amp/"iter.$NUM_ITERATION"" ]; 
	then 
		mv "$RESULTS"/SSU_sequences/output_MetaRib/Iteration/iter_1/emirge_amp/"iter.$NUM_ITERATION" "$RESULTS"/SSU_sequences/output_MetaRib/Iteration/iter_1/emirge_amp/"last.$NUM_ITERATION"
		rm -rf "$RESULTS"/SSU_sequences/output_MetaRib/Iteration/iter_1/emirge_amp/*iter*; 
		mv "$RESULTS"/SSU_sequences/output_MetaRib/Iteration/iter_1/emirge_amp/"last.$NUM_ITERATION" "$RESULTS"/SSU_sequences/output_MetaRib/Iteration/iter_1/emirge_amp/"iter.$NUM_ITERATION"
fi

if [ -d "$RESULTS/SSU_sequences/output_MetaRib/Iteration/iter_2_L/emirge_amp/"iter.$NUM_ITERATION"" ]; 
	then 
		mv "$RESULTS"/SSU_sequences/output_MetaRib/Iteration/iter_2_L/emirge_amp/"iter.$NUM_ITERATION" "$RESULTS"/SSU_sequences/output_MetaRib/Iteration/iter_2_L/emirge_amp/"last.$NUM_ITERATION"
		rm -rf "$RESULTS"/SSU_sequences/output_MetaRib/Iteration/iter_2_L/emirge_amp/*iter*; 
		mv "$RESULTS"/SSU_sequences/output_MetaRib/Iteration/iter_2_L/emirge_amp/"last.$NUM_ITERATION" "$RESULTS"/SSU_sequences/output_MetaRib/Iteration/iter_2_L/emirge_amp/"iter.$NUM_ITERATION"
fi


#if [ -d "$RESULTS/output_MetaRib/Iteration/iter_1/emirge_amp/" ]; then rm -Rf $RESULTS/output_MetaRib/Iteration/iter_1/emirge_amp/!("initial_mapping"|"iter.$NUM_ITERATION"); fi
#if [ -d "$RESULTS/output_MetaRib/Iteration/iter_2_L/emirge_amp/" ]; then rm -Rf $RESULTS/output_MetaRib/Iteration/iter_2_L/emirge_amp/!("initial_mapping"|"iter.$NUM_ITERATION"); fi

#mkdir -p "$RESULTS/SSU_sequences/output_MetaRib"
#mkdir -p "$RESULTS/SSU_sequences/output_MetaRib/$SHORTNAME"
#mkdir -p "$RESULTS/SSU_sequences/Abundance_calc"

#mv "$RESULTS"/output_MetaRib/Abundance "$RESULTS"/SSU_sequences/output_MetaRib/"$SHORTNAME"
#mv "$RESULTS"/output_MetaRib/Iteration "$RESULTS"/SSU_sequences/output_MetaRib/"$SHORTNAME"
mv "$RESULTS"/SSU_sequences/output_MetaRib/Abundance/all.dedup.fasta "$RESULTS"/SSU_sequences/output_MetaRib
mv "$RESULTS"/SSU_sequences/output_MetaRib/all.dedup.fasta "$RESULTS"/SSU_sequences/output_MetaRib/MetaRib_SSU.fasta

#rm -r "$RESULTS"/output_MetaRib
#rm "$RESULTS"/SSU_sequences/output_MetaRib/"$SHORTNAME"/samples.list.txt


#cat "$RESULTS"/SSU_sequences/output_emirge/"$SHORTNAME"_renamed_16S18S_recons.fasta "$RESULTS"/SSU_sequences/output_MetaRib/"$SHORTNAME"/MetaRib_SSU.fasta > "$RESULTS"/SSU_sequences/emirge_metarib_SSU_sequences.fasta

#clustering at 97%
#vsearch --cluster_fast "$RESULTS"/SSU_sequences/emirge_metarib_SSU_sequences.fasta --centroids "$RESULTS"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.fasta --id 0.97

#covert small letters into capital letters
#awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' "$RESULTS"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.fasta > "$RESULTS"/SSU_sequences/"$SHORTNAME"_SSU_sequences.fasta

#rm -rv "$RESULTS"/output_MetaRib/"$SHORTNAME"/Iteration/iter_1/emirge_amp/!("initial_mapping"|"iter.$NUM_ITERATION")


#zip -r Iteration.zip "$RESULTS"/output_MetaRib/"$SHORTNAME"/Iteration/

#rm -r "$RESULTS"/output_MetaRib/"$SHORTNAME"/Iteration

#mv Iteration.zip "$RESULTS"/output_MetaRib/"$SHORTNAME"

echo "Saving results..." | tee /dev/fd/3

fi

echo ">Reconstructing 16S/18S sequences ends successfully on : "`date` | tee /dev/fd/3

echo "Finalising reconstructed sequences..." | tee /dev/fd/3

cat "$RESULTS"/SSU_sequences/output_emirge/"$SHORTNAME"_renamed_16S18S_recons.fasta "$RESULTS"/SSU_sequences/output_MetaRib/MetaRib_SSU.fasta |sed '/>/s/[ ].*//' > "$RESULTS"/SSU_sequences/"$SHORTNAME"_emirge_metarib_SSU_sequences.fasta


echo ">Reconstructing 16S/18S sequences ends successfully on : "`date` | tee /dev/fd/3


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#			Abundance calculation
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

echo "" | tee /dev/fd/3
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo ">Relative abundance calculation starting on : "`date` | tee /dev/fd/3
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" | tee /dev/fd/3
echo "" | tee /dev/fd/3

ABUNDANCE_FILE=""$RESULTS"/SSU_sequences/"$SHORTNAME"_Abundance.tsv"

if [ -f "$ABUNDANCE_FILE" ];
then
	echo "Abundance file : '${ABUNDANCE_FILE}' is present."

else

mkdir -p ""$RESULTS"/SSU_sequences/Abundance_calc"

echo ">Calculating relative abundances of reconstructed sequences..." | tee /dev/fd/3

echo "command line :" | tee /dev/fd/3
echo "bbmap.sh -Xmx${RAM}g in=$RESULTS/quality_control/"$SHORTNAME"_trimmed.fastq \
	ref=$RESULTS/SSU_sequences/"$SHORTNAME"_emirge_metarib_SSU_sequences.fasta \
	path=$RESULTS/SSU_sequences/ \
	covstats=$RESULTS/SSU_sequences/"$SHORTNAME"_covstats.txt \
	threads=$THREAD \
	32bit=t \
	scafstats=$RESULTS/SSU_sequences/"$SHORTNAME"_scafstats.txt " | tee /dev/fd/3

bbmap.sh -Xmx${RAM}g in="$RESULTS"/quality_control/"$SHORTNAME"_trimmed.fastq \
	ref="$RESULTS"/SSU_sequences/"$SHORTNAME"_emirge_metarib_SSU_sequences.fasta \
	path="$RESULTS"/SSU_sequences/ \
	covstats="$RESULTS"/SSU_sequences/"$SHORTNAME"_covstats.txt \
	scafstats="$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats.txt \
	threads=$THREAD \
	32bit=t


vsearch --cluster_fast "$RESULTS"/SSU_sequences/"$SHORTNAME"_emirge_metarib_SSU_sequences.fasta --centroids "$RESULTS"/SSU_sequences/"$SHORTNAME"_emirge_metarib_clustered_SSU_sequences.fasta --id 0.97 --uc "$RESULTS"/SSU_sequences/"$SHORTNAME"_emirge_metarib_clustered_SSU_sequences.tsv

awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' "$RESULTS"/SSU_sequences/"$SHORTNAME"_emirge_metarib_clustered_SSU_sequences.fasta > "$RESULTS"/SSU_sequences/"$SHORTNAME"_SSU_sequences.fasta

cat "$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats.txt  | sed 1d | awk '{ print $1,$8 }' | sort -k1 | tr ' ' \\t > "$RESULTS"/SSU_sequences/"$SHORTNAME"_sorted_emirge_metarib_scafstats.txt
cat "$RESULTS"/SSU_sequences/"$SHORTNAME"_emirge_metarib_clustered_SSU_sequences.tsv |awk '!($1=="C"){print}'|awk '{ print $2,$9 }' |sort -k2 |awk '{ print $2,$1 }' > "$RESULTS"/SSU_sequences/"$SHORTNAME"_sorted_emirge_metarib_clustered_SSU_sequences.tsv
join "$RESULTS"/SSU_sequences/"$SHORTNAME"_sorted_emirge_metarib_scafstats.txt "$RESULTS"/SSU_sequences/"$SHORTNAME"_sorted_emirge_metarib_clustered_SSU_sequences.tsv |tr  ' ' \\t > "$RESULTS"/SSU_sequences/"$SHORTNAME"_reads_count.tsv

cat "$RESULTS"/SSU_sequences/"$SHORTNAME"_reads_count.tsv | awk '{ print $3,$2 }' |sort -k1 |awk '{ sum[$1] += $2; count[$1] += 1 } END { for ( key in count ) { print key, count[key], sum[key] } }' |tr  ' ' \\t |sort -k1 -n > "$RESULTS"/SSU_sequences/"$SHORTNAME"_reads_count_byCluster.tsv
cat "$RESULTS"/SSU_sequences/"$SHORTNAME"_emirge_metarib_clustered_SSU_sequences.tsv|awk '($1=="C"){print}'|awk '{ print $2,$9 }'|tr  ' ' \\t | sort -k1 -n > "$RESULTS"/SSU_sequences/"$SHORTNAME"_Cluster_refseq.tsv
join "$RESULTS"/SSU_sequences/"$SHORTNAME"_reads_count_byCluster.tsv "$RESULTS"/SSU_sequences/"$SHORTNAME"_Cluster_refseq.tsv  |awk '{ print $4,$1,$2,$3 }'|tr  ' ' \\t |sort -k1 > "$RESULTS"/SSU_sequences/"$SHORTNAME"_reads_count_bySeq.tsv


cat "$RESULTS"/SSU_sequences/"$SHORTNAME"_reads_count_bySeq.tsv | awk '{ print $1 }' |sort -k1 > "$RESULTS"/SSU_sequences/"$SHORTNAME"_id.tsv
cat "$RESULTS"/SSU_sequences/"$SHORTNAME"_covstats.txt |sed 1d |awk '{ print $1,$3 }' |tr  ' ' \\t |sort -k1 > "$RESULTS"/SSU_sequences/"$SHORTNAME"_table.tsv

awk -F '\t' 'NR==FNR {id[$1]; next} $1 in id' "$RESULTS"/SSU_sequences/"$SHORTNAME"_id.tsv "$RESULTS"/SSU_sequences/"$SHORTNAME"_table.tsv | tr  ' ' \\t |sort -k1 > "$RESULTS"/SSU_sequences/"$SHORTNAME"_length_bySeq.tsv
join "$RESULTS"/SSU_sequences/"$SHORTNAME"_reads_count_bySeq.tsv "$RESULTS"/SSU_sequences/"$SHORTNAME"_length_bySeq.tsv |awk '{ print $1,$4,$5,$2,$3 }' |awk '{sub("0","1", $2)}1' |tr  ' ' \\t > "$RESULTS"/SSU_sequences/"$SHORTNAME"_readsCount_length.tsv

total=$(awk '{s+=$2}END{print s}' "$RESULTS"/SSU_sequences/"$SHORTNAME"_readsCount_length.tsv)

awk -v total=$total '{ printf ("%s\t%s\t%.6f\n", $1, $2, ($2/total)*100)}' "$RESULTS"/SSU_sequences/"$SHORTNAME"_readsCount_length.tsv | tr  ' ' \\t |sort -k1  > "$RESULTS"/SSU_sequences/"$SHORTNAME"_RA_length.tsv

join "$RESULTS"/SSU_sequences/"$SHORTNAME"_RA_length.tsv "$RESULTS"/SSU_sequences/"$SHORTNAME"_length_bySeq.tsv |awk '{ print $1,$4,$2,$3 }' | awk '{sub("0","1", $2)}1' |tr  ' ' \\t |sort -k1 |awk 'BEGIN{print "Seq_ID\tlength\treads_count\tRelative_abundance"}1' > "$RESULTS"/SSU_sequences/"$SHORTNAME"_Abundance.tsv

mv "$RESULTS"/SSU_sequences/"$SHORTNAME"_emirge_metarib_clustered_SSU_sequences.tsv "$RESULTS"/SSU_sequences/Abundance_calc/
mv "$RESULTS"/SSU_sequences/Abundance_calc/"$SHORTNAME"_emirge_metarib_clustered_SSU_sequences.tsv "$RESULTS"/SSU_sequences/Abundance_calc/"$SHORTNAME"_allclusters.tsv
mv "$RESULTS"/SSU_sequences/"$SHORTNAME"_readsCount_length.tsv "$RESULTS"/SSU_sequences/Abundance_calc/
mv "$RESULTS"/SSU_sequences/Abundance_calc/"$SHORTNAME"_readsCount_length.tsv "$RESULTS"/SSU_sequences/Abundance_calc/"$SHORTNAME"_RCountsLen_byCluster.tsv
mv "$RESULTS"/SSU_sequences/"$SHORTNAME"_reads_count.tsv "$RESULTS"/SSU_sequences/Abundance_calc/
mv "$RESULTS"/SSU_sequences/Abundance_calc/"$SHORTNAME"_reads_count.tsv "$RESULTS"/SSU_sequences/Abundance_calc/"$SHORTNAME"_Rcounts_allClusters.tsv
mv "$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats.txt "$RESULTS"/SSU_sequences/Abundance_calc/
mv "$RESULTS"/SSU_sequences/"$SHORTNAME"_covstats.txt "$RESULTS"/SSU_sequences/Abundance_calc/

#cleaning
rm "$RESULTS"/SSU_sequences/"$SHORTNAME"_sorted_emirge_metarib_scafstats.txt 
rm "$RESULTS"/SSU_sequences/"$SHORTNAME"_sorted_emirge_metarib_clustered_SSU_sequences.tsv
rm "$RESULTS"/SSU_sequences/"$SHORTNAME"_Cluster_refseq.tsv
rm "$RESULTS"/SSU_sequences/"$SHORTNAME"_reads_count_byCluster.tsv
rm "$RESULTS"/SSU_sequences/"$SHORTNAME"_table.tsv 
rm "$RESULTS"/SSU_sequences/"$SHORTNAME"_id.tsv
rm "$RESULTS"/SSU_sequences/"$SHORTNAME"_reads_count_bySeq.tsv
rm "$RESULTS"/SSU_sequences/"$SHORTNAME"_length_bySeq.tsv
rm "$RESULTS"/SSU_sequences/"$SHORTNAME"_RA_length.tsv
rm "$RESULTS"/SSU_sequences/"$SHORTNAME"_emirge_metarib_SSU_sequences.fasta
rm "$RESULTS"/SSU_sequences/"$SHORTNAME"_emirge_metarib_clustered_SSU_sequences.fasta

rm -rf "$RESULTS"/SSU_sequences/ref

#bbmap.sh -Xmx3g in1="$RESULTS"/quality_control/"$SHORTNAME"_1_trimmed.fastq \
#	in2="$RESULTS"/quality_control/"$SHORTNAME"_2_trimmed.fastq \
#	ref="$RESULTS"/SSU_sequences/"$SHORTNAME"_SSU_sequences.fasta \
#	scafstats="$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats.txt \
#	covstats="$RESULTS"/SSU_sequences/"$SHORTNAME"_covstats.txt


#cat "$RESULTS"/SSU_sequences/all_scafstats.txt | sed 1d | tr ',' \\t | awk '{ print $1,$8 }' | sort -k1 -k2 | awk '!($2==0){print}' | awk '{print $1}' > "$RESULTS"/SSU_sequences/id_file.txt

#cat "$RESULTS"/SSU_sequences/all_scafstats.txt | sed 1d | tr ',' \\t | awk '{ print $1,$8 }' | sort -k1 -k2 | awk '!($2==0){print}' > "$RESULTS"/SSU_sequences/"$SHORTNAME"_scafstats.txt

#awk 'NR==FNR{ids[$0]; next} ($1 in ids){ printf ">" $0 }' "$RESULTS"/SSU_sequences/id_file.txt RS='>' "$RESULTS"/SSU_sequences/all_SSU_sequences.fasta > "$RESULTS"/SSU_sequences/"$SHORTNAME"_SSU_sequences.fasta

#awk -F '\t' 'NR==FNR {id[$1]; next} $1 in id' "$RESULTS"/SSU_sequences/id_file.txt "$RESULTS"/SSU_sequences/all_covstats.txt > "$RESULTS"/SSU_sequences/"$SHORTNAME"_covstats.txt

#rm "$RESULTS"/SSU_sequences/all_scafstats.txt
#rm "$RESULTS"/SSU_sequences/all_covstats.txt
#rm "$RESULTS"/SSU_sequences/id_file.txt
#rm "$RESULTS"/SSU_sequences/all_SSU_sequences.fasta


#rm "$RESULTS"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.fasta
#rm -r "$RESULTS"/SSU_sequences/output_MetaRib/"$SHORTNAME"/Abundance

echo "Saving results..." | tee /dev/fd/3

fi

echo ">Relative abundance calculation ends successfully on : "`date` | tee /dev/fd/3

conda deactivate
#echo "RiboTaxa virtual environment has been deactivated successfully..." | tee /dev/fd/3
