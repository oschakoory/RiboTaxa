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

#fomrat of files (can either be fastq or fastq.gz)
FORMAT=$(awk '/^FORMAT/{print $3}' "${CONFIG}")

#Set up results sub-directories
mkdir -p "$OUTPUT/quality_control"
mkdir -p "$OUTPUT/quality_control/before_fastqc"
mkdir -p "$OUTPUT/quality_control/before_fastqc/multiqc"
mkdir -p "$OUTPUT/quality_control/after_fastqc"
mkdir -p "$OUTPUT/quality_control/after_fastqc/multiqc"

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

echo "Run the First quality control on raw data... " | tee /dev/fd/3

echo "command line :" | tee /dev/fd/3
echo "fastqc $DATA_DIR/"$SHORTNAME"_R1.$FORMAT  $DATA_DIR/"$SHORTNAME"_R2.$FORMAT  -dir $OUTPUT -o $OUTPUT/quality_control/before_fastqc" | tee /dev/fd/3

fastqc "$DATA_DIR"/"$SHORTNAME"_R1.$FORMAT  "$DATA_DIR"/"$SHORTNAME"_R2.$FORMAT  -dir $OUTPUT -o "$OUTPUT"/quality_control/before_fastqc

echo "Removing adapters from sequences..." | tee /dev/fd/3
echo "command line :" | tee /dev/fd/3
echo "bbduk.sh -Xmx${RAM}g \
	in1=$DATA_DIR/"$SHORTNAME"_R1.$FORMAT  \
	in2=$DATA_DIR/"$SHORTNAME"_R2.$FORMAT  \
	out1=$OUTPUT/quality_control/"$SHORTNAME"_1_noadapt.$FORMAT  \
	out2=$OUTPUT/quality_control/"$SHORTNAME"_2_noadapt.$FORMAT  \
	ref=$RiboTaxa_DIR/adapters/TruSeq3-PE.fa \
	ktrim=$KTRIM \
	k=$KMER \
	mink=11 \
	tpe \
	tbo" | tee /dev/fd/3

bbduk.sh -Xmx${RAM}g \
	in1="$DATA_DIR"/"$SHORTNAME"_R1.$FORMAT  \
	in2="$DATA_DIR"/"$SHORTNAME"_R2.$FORMAT  \
	out1="$OUTPUT"/quality_control/"$SHORTNAME"_1_noadapt.$FORMAT  \
	out2="$OUTPUT"/quality_control/"$SHORTNAME"_2_noadapt.$FORMAT  \
	ref="$RiboTaxa_DIR"/adapters/TruSeq3-PE.fa \
	ktrim=$KTRIM \
	k=$KMER \
	mink=11 \
	tpe \
	tbo


echo "Trimming sequences..." | tee /dev/fd/3
echo "command line :" | tee /dev/fd/3
echo "bbduk.sh -Xmx${RAM}g \
	in1=$OUTPUT/quality_control/"$SHORTNAME"_1_noadapt.$FORMAT  \
	in2=$OUTPUT/quality_control/"$SHORTNAME"_2_noadapt.$FORMAT  \
	out1=$OUTPUT/quality_control/"$SHORTNAME"_1_trimmed.$FORMAT  \
	out2=$OUTPUT/quality_control/"$SHORTNAME"_2_trimmed.$FORMAT  \
	minlen=$MINLENGTH \
	qtrim=$QTRIM \
	trimq=$TRIMQ \
	maxns=$MAXNS" | tee /dev/fd/3

bbduk.sh -Xmx${RAM}g \
	in1="$OUTPUT"/quality_control/"$SHORTNAME"_1_noadapt.$FORMAT  \
	in2="$OUTPUT"/quality_control/"$SHORTNAME"_2_noadapt.$FORMAT  \
	out1="$OUTPUT"/quality_control/"$SHORTNAME"_1_trimmed.$FORMAT  \
	out2="$OUTPUT"/quality_control/"$SHORTNAME"_2_trimmed.$FORMAT  \
	minlen=$MINLENGTH \
	qtrim=$QTRIM \
	trimq=$TRIMQ \
	maxns=$MAXNS


echo "Run the second quality control on trimmed data..." | tee /dev/fd/3

echo "command line :" | tee /dev/fd/3
echo "fastqc $OUTPUT/quality_control/"$SHORTNAME"_1_trimmed.$FORMAT  $OUTPUT/quality_control/"$SHORTNAME"_2_trimmed.$FORMAT  -dir $OUTPUT -o $OUTPUT/quality_control/after_fastqc" | tee /dev/fd/3

fastqc "$OUTPUT"/quality_control/"$SHORTNAME"_1_trimmed.$FORMAT  "$OUTPUT"/quality_control/"$SHORTNAME"_2_trimmed.$FORMAT  -dir $OUTPUT -o "$OUTPUT"/quality_control/after_fastqc

echo "Running multiQC..." | tee /dev/fd/3

multiqc -f "$OUTPUT"/quality_control/before_fastqc/* --outdir "$OUTPUT"/quality_control/before_fastqc/multiqc/

multiqc -f "$OUTPUT"/quality_control/after_fastqc/* --outdir "$OUTPUT"/quality_control/after_fastqc/multiqc/

echo "Saving results for quality control..." | tee /dev/fd/3

echo ">Quality control ends successfully on : "`date` | tee /dev/fd/3
 
rm "$OUTPUT"/quality_control/*_noadapt.$FORMAT 

#unzip fastq files if compressed as sortmerna takes uncompressed files only
for NAME in `ls "$OUTPUT"/quality_control/*"$SHORTNAME"*`; do
	FILE=($(basename ""${NAME[@]}""))
	#echo $FILE
	if [[ ${NAME[@]} == *.gz* ]]; then
		#echo "$FILE is gzipped"
		gunzip -f $NAME
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
echo "command line :" | tee /dev/fd/3
echo "reformat.sh in1=$OUTPUT/quality_control/"$SHORTNAME"_1_trimmed.fastq in2=$OUTPUT/quality_control/"$SHORTNAME"_2_trimmed.fastq out=$OUTPUT/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq overwrite=t" | tee /dev/fd/3

#bash merge-paired-reads.sh "$OUTPUT"/quality_control/"$SHORTNAME"_1_trimmed.fastq "$OUTPUT"/quality_control/"$SHORTNAME"_2_trimmed.fastq "$OUTPUT"/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq
reformat.sh in1="$OUTPUT"/quality_control/"$SHORTNAME"_1_trimmed.fastq in2="$OUTPUT"/quality_control/"$SHORTNAME"_2_trimmed.fastq out="$OUTPUT"/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq overwrite=t


echo "Filtering 16S/18S reads...." | tee /dev/fd/3
echo "command line :" | tee /dev/fd/3
echo "sortmerna --ref $SORTMERNA_DB/$SORTME_NAME.fasta,$SORTMERNA_DB/$SORTME_NAME \
	--reads $OUTPUT/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq \
	--fastx \
	--aligned $OUTPUT/output_sortmerna/"$SHORTNAME"_16S18S \
	--paired_in \
	-a $THREAD \
	--log \
	-v" | tee /dev/fd/3

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
echo "command line :" | tee /dev/fd/3
echo "reformat.sh in=$OUTPUT/output_sortmerna/"$SHORTNAME"_16S18S.fastq out1=$OUTPUT/output_sortmerna/"$SHORTNAME"_R1_16S18Sreads.fastq out2=$OUTPUT/output_sortmerna/"$SHORTNAME"_R2_16S18Sreads.fastq overwrite=t" | tee /dev/fd/3

reformat.sh in="$OUTPUT"/output_sortmerna/"$SHORTNAME"_16S18S.fastq out1="$OUTPUT"/output_sortmerna/"$SHORTNAME"_R1_16S18Sreads.fastq out2="$OUTPUT"/output_sortmerna/"$SHORTNAME"_R2_16S18Sreads.fastq overwrite=t


echo "Saving results..." | tee /dev/fd/3

echo ">Filtering 16S/18S using sortmerna ends successfully on : "`date` | tee /dev/fd/3

rm "$OUTPUT"/output_sortmerna/"$SHORTNAME"_16S18S.fastq
rm "$OUTPUT"/output_sortmerna/"$SHORTNAME"_mergedpaired.fastq

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
mkdir -p "$OUTPUT/SSU_sequences/output_emirge"
rm -rf "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/*

#echo "Setting up indexed database for EMIRGE..." | tee /dev/fd/3
EMIRGE_DB=$(awk '/^EMIRGE_DB/{print $3}' "${CONFIG}")
echo "Database directory for EMIRGE = $EMIRGE_DB" | tee /dev/fd/3

NAME=($(ls "$EMIRGE_DB"/*.fasta*)) 
#echo "name = " | tee /dev/fd/3
#printf '%s\n' "${NAME[@]}" | tee /dev/fd/3

REF_NAME=$(basename ""${NAME[@]}"")  

NAME=($(ls "$EMIRGE_DB"/*.4.ebwt* | sed 's/.4.ebwt//')) 
#echo "name = " | tee /dev/fd/3
#printf '%s\n' "${NAME[@]}" | tee /dev/fd/3

BWT_NAME=$(basename ""${NAME[@]}"")  

MAX_LENGTH=$(awk '/^MAX_LENGTH/{print $3}' "${CONFIG}")
IDENTITY=$(awk '/^IDENTITY/{print $3}' "${CONFIG}")
NUM_ITERATION=$(awk '/^NUM_ITERATION/{print $3}' "${CONFIG}")
MEAN_INSERT_SIZE=$(awk '/^MEAN_INSERT_SIZE/{print $3}' "${CONFIG}")
STD_DEV=$(awk '/^STD_DEV/{print $3}' "${CONFIG}")

echo "command line :" | tee /dev/fd/3
echo "emirge_amplicon.py \
	-1 $OUTPUT/output_sortmerna/"$SHORTNAME"_R1_16S18Sreads.fastq \
	-f $EMIRGE_DB/$REF_NAME \
	-b $EMIRGE_DB/$BWT_NAME \
	-l $MAX_LENGTH \
	-2 $OUTPUT/output_sortmerna/"$SHORTNAME"_R2_16S18Sreads.fastq \
	-j $IDENTITY \
	-n $NUM_ITERATION \
	-i $MEAN_INSERT_SIZE \
	-s $STD_DEV \
	-a $THREAD \
	--phred33 $OUTPUT/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons" | tee /dev/fd/3


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
echo "command line :" | tee /dev/fd/3
echo "emirge_rename_fasta.py --no_N $OUTPUT/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/iter.$NUM_ITERATION > $OUTPUT/SSU_sequences/output_emirge/"$SHORTNAME"_renamed_16S18S_recons.fasta" | tee /dev/fd/3

emirge_rename_fasta.py --no_N "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/iter.$NUM_ITERATION > "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_renamed_16S18S_recons.fasta

#rm -rv "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/!("initial_mapping"|"iter.$NUM_ITERATION")

#zip -jr  "$SHORTNAME"_amplicon_16S18S_recons.zip "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons/*

cd "$OUTPUT"/SSU_sequences/output_emirge && zip -qrm "$SHORTNAME"_amplicon_16S18S_recons.zip "$SHORTNAME"_amplicon_16S18S_recons/ && cd -

#rm -r "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_amplicon_16S18S_recons

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#			Reconstructing 16S/18S full length sequences using MetaRib
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


echo ">Running MetaRib to reconstruct 16S/18S full length sequences..."`date` | tee /dev/fd/3

#SAMPLE=$(awk '{s++}END{print s/4}' "$OUTPUT"/output_sortmerna/"$SHORTNAME"_R1_16S18Sreads.fastq)


echo $SHORTNAME > "$OUTPUT"/quality_control/samples.list.txt

echo "command line :" | tee /dev/fd/3
echo "python2 $RiboTaxa_DIR/run_MetaRib.py -cfg $CONFIG_PATH -p $OUTPUT/quality_control -b $EMIRGE_DB/$BWT_NAME -l $EMIRGE_DB/$REF_NAME" | tee /dev/fd/3

python2 "$RiboTaxa_DIR"/run_MetaRib.py -cfg "$CONFIG_PATH" -p "$OUTPUT"/quality_control -b "$EMIRGE_DB"/$BWT_NAME -l "$EMIRGE_DB"/"$REF_NAME"

cd "$OUTPUT"/output_MetaRib && zip -qrm Iteration.zip Iteration/ && cd -

mkdir -p "$OUTPUT/SSU_sequences/output_MetaRib"
mkdir -p "$OUTPUT/SSU_sequences/output_MetaRib/$SHORTNAME"
mkdir -p "$OUTPUT/SSU_sequences/Abundance_calc"

mv "$OUTPUT"/output_MetaRib/Abundance "$OUTPUT"/SSU_sequences/output_MetaRib/"$SHORTNAME"
mv "$OUTPUT"/output_MetaRib/Iteration.zip "$OUTPUT"/SSU_sequences/output_MetaRib/"$SHORTNAME"
mv "$OUTPUT"/SSU_sequences/output_MetaRib/"$SHORTNAME"/Abundance/all.dedup.fasta "$OUTPUT"/SSU_sequences/output_MetaRib/"$SHORTNAME"
mv "$OUTPUT"/SSU_sequences/output_MetaRib/"$SHORTNAME"/all.dedup.fasta "$OUTPUT"/SSU_sequences/output_MetaRib/"$SHORTNAME"/MetaRib_SSU.fasta

rm -r "$OUTPUT"/output_MetaRib
rm "$OUTPUT"/quality_control/samples.list.txt

echo "Finalising reconstructed sequences..." | tee /dev/fd/3

cat "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_renamed_16S18S_recons.fasta "$OUTPUT"/SSU_sequences/output_MetaRib/"$SHORTNAME"/MetaRib_SSU.fasta |sed '/>/s/[ ].*//' > "$OUTPUT"/SSU_sequences/emirge_metarib_SSU_sequences.fasta

#cat "$OUTPUT"/SSU_sequences/output_emirge/"$SHORTNAME"_renamed_16S18S_recons.fasta "$OUTPUT"/SSU_sequences/output_MetaRib/"$SHORTNAME"/MetaRib_SSU.fasta > "$OUTPUT"/SSU_sequences/emirge_metarib_SSU_sequences.fasta

#clustering at 97%
#vsearch --cluster_fast "$OUTPUT"/SSU_sequences/emirge_metarib_SSU_sequences.fasta --centroids "$OUTPUT"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.fasta --id 0.97

#covert small letters into capital letters
#awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' "$OUTPUT"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.fasta > "$OUTPUT"/SSU_sequences/"$SHORTNAME"_SSU_sequences.fasta

#rm -rv "$OUTPUT"/output_MetaRib/"$SHORTNAME"/Iteration/iter_1/emirge_amp/!("initial_mapping"|"iter.$NUM_ITERATION")


#zip -r Iteration.zip "$OUTPUT"/output_MetaRib/"$SHORTNAME"/Iteration/

#rm -r "$OUTPUT"/output_MetaRib/"$SHORTNAME"/Iteration

#mv Iteration.zip "$OUTPUT"/output_MetaRib/"$SHORTNAME"

echo "Saving results..." | tee /dev/fd/3

echo ">Reconstructing 16S/18S sequences ends successfully on : "`date` | tee /dev/fd/3


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#			Abundance calculation
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

echo ">Calculating relative abundances of reconstructed sequences..." | tee /dev/fd/3

echo "command line :" | tee /dev/fd/3
echo "bbmap.sh -Xmx${RAM}g in1=$OUTPUT/quality_control/"$SHORTNAME"_1_trimmed.fastq \
	in2=$OUTPUT/quality_control/"$SHORTNAME"_2_trimmed.fastq \
	ref=$OUTPUT/SSU_sequences/emirge_metarib_SSU_sequences.fasta \
	path=$OUTPUT/SSU_sequences/ \
	covstats=$OUTPUT/SSU_sequences/"$SHORTNAME"_covstats.txt \
	scafstats=$OUTPUT/SSU_sequences/"$SHORTNAME"_scafstats.txt " | tee /dev/fd/3

bbmap.sh -Xmx${RAM}g in1="$OUTPUT"/quality_control/"$SHORTNAME"_1_trimmed.fastq \
	in2="$OUTPUT"/quality_control/"$SHORTNAME"_2_trimmed.fastq \
	ref="$OUTPUT"/SSU_sequences/emirge_metarib_SSU_sequences.fasta \
	path="$OUTPUT"/SSU_sequences/ \
	covstats="$OUTPUT"/SSU_sequences/"$SHORTNAME"_covstats.txt \
	scafstats="$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats.txt \
	32bit=t


vsearch --cluster_fast "$OUTPUT"/SSU_sequences/emirge_metarib_SSU_sequences.fasta --centroids "$OUTPUT"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.fasta --id 0.97 --uc "$OUTPUT"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.tsv

awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' "$OUTPUT"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.fasta > "$OUTPUT"/SSU_sequences/"$SHORTNAME"_SSU_sequences.fasta

cat "$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats.txt  | sed 1d | awk '{ print $1,$8 }' | sort -k1 | tr ' ' \\t > "$OUTPUT"/SSU_sequences/sorted_emirge_metarib_scafstats.txt
cat "$OUTPUT"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.tsv |awk '!($1=="C"){print}'|awk '{ print $2,$9 }' |sort -k2 |awk '{ print $2,$1 }' > "$OUTPUT"/SSU_sequences/sorted_emirge_metarib_clustered_SSU_sequences.tsv
join "$OUTPUT"/SSU_sequences/sorted_emirge_metarib_scafstats.txt "$OUTPUT"/SSU_sequences/sorted_emirge_metarib_clustered_SSU_sequences.tsv |tr  ' ' \\t > "$OUTPUT"/SSU_sequences/reads_count.tsv

cat "$OUTPUT"/SSU_sequences/reads_count.tsv | awk '{ print $3,$2 }' |sort -k1 |awk '{ sum[$1] += $2; count[$1] += 1 } END { for ( key in count ) { print key, count[key], sum[key] } }' |tr  ' ' \\t |sort -k1 -n > "$OUTPUT"/SSU_sequences/reads_count_byCluster.tsv
cat "$OUTPUT"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.tsv|awk '($1=="C"){print}'|awk '{ print $2,$9 }'|tr  ' ' \\t | sort -k1 -n > "$OUTPUT"/SSU_sequences/Cluster_refseq.tsv
join "$OUTPUT"/SSU_sequences/reads_count_byCluster.tsv "$OUTPUT"/SSU_sequences/Cluster_refseq.tsv  |awk '{ print $4,$1,$2,$3 }'|tr  ' ' \\t |sort -k1 > "$OUTPUT"/SSU_sequences/reads_count_bySeq.tsv


cat "$OUTPUT"/SSU_sequences/reads_count_bySeq.tsv | awk '{ print $1 }' |sort -k1 > "$OUTPUT"/SSU_sequences/id.tsv
cat "$OUTPUT"/SSU_sequences/"$SHORTNAME"_covstats.txt |sed 1d |awk '{ print $1,$3 }' |tr  ' ' \\t |sort -k1 > "$OUTPUT"/SSU_sequences/table.tsv

awk -F '\t' 'NR==FNR {id[$1]; next} $1 in id' "$OUTPUT"/SSU_sequences/id.tsv "$OUTPUT"/SSU_sequences/table.tsv | tr  ' ' \\t |sort -k1 > "$OUTPUT"/SSU_sequences/length_bySeq.tsv
join "$OUTPUT"/SSU_sequences/reads_count_bySeq.tsv "$OUTPUT"/SSU_sequences/length_bySeq.tsv |awk '{ print $1,$4,$5,$2,$3 }' |tr  ' ' \\t > "$OUTPUT"/SSU_sequences/readsCount_length.tsv

total=$(awk '{s+=$2}END{print s}' "$OUTPUT"/SSU_sequences/readsCount_length.tsv)

awk -v total=$total '{ printf ("%s\t%s\t%.6f\n", $1, $2, ($2/total)*100)}' "$OUTPUT"/SSU_sequences/readsCount_length.tsv | tr  ' ' \\t |sort -k1  > "$OUTPUT"/SSU_sequences/RA_length.tsv

join "$OUTPUT"/SSU_sequences/RA_length.tsv "$OUTPUT"/SSU_sequences/length_bySeq.tsv |awk '{ print $1,$4,$2,$3 }'|tr  ' ' \\t |sort -k1 |awk 'BEGIN{print "Seq_ID\tlength\treads_count\tRelative_abundance"}1' > "$OUTPUT"/SSU_sequences/"$SHORTNAME"_Abundance.tsv

mv "$OUTPUT"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.tsv "$OUTPUT"/SSU_sequences/Abundance_calc
mv "$OUTPUT"/SSU_sequences/Abundance_calc/emirge_metarib_clustered_SSU_sequences.tsv "$OUTPUT"/SSU_sequences/Abundance_calc/"$SHORTNAME"_allclusters.tsv
mv "$OUTPUT"/SSU_sequences/readsCount_length.tsv "$OUTPUT"/SSU_sequences/Abundance_calc
mv "$OUTPUT"/SSU_sequences/Abundance_calc/readsCount_length.tsv "$OUTPUT"/SSU_sequences/Abundance_calc/"$SHORTNAME"_RCountsLen_byCluster.tsv
mv "$OUTPUT"/SSU_sequences/reads_count.tsv "$OUTPUT"/SSU_sequences/Abundance_calc
mv "$OUTPUT"/SSU_sequences/Abundance_calc/reads_count.tsv "$OUTPUT"/SSU_sequences/Abundance_calc/"$SHORTNAME"_Rcounts_allClusters.tsv
mv "$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats.txt "$OUTPUT"/SSU_sequences/Abundance_calc
mv "$OUTPUT"/SSU_sequences/"$SHORTNAME"_covstats.txt "$OUTPUT"/SSU_sequences/Abundance_calc

#cleaning
rm "$OUTPUT"/SSU_sequences/sorted_emirge_metarib_scafstats.txt 
rm "$OUTPUT"/SSU_sequences/sorted_emirge_metarib_clustered_SSU_sequences.tsv
rm "$OUTPUT"/SSU_sequences/Cluster_refseq.tsv
rm "$OUTPUT"/SSU_sequences/reads_count_byCluster.tsv
rm "$OUTPUT"/SSU_sequences/table.tsv 
rm "$OUTPUT"/SSU_sequences/id.tsv
rm "$OUTPUT"/SSU_sequences/reads_count_bySeq.tsv
rm "$OUTPUT"/SSU_sequences/length_bySeq.tsv
rm "$OUTPUT"/SSU_sequences/RA_length.tsv
rm "$OUTPUT"/SSU_sequences/emirge_metarib_SSU_sequences.fasta
rm "$OUTPUT"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.fasta


#bbmap.sh -Xmx3g in1="$OUTPUT"/quality_control/"$SHORTNAME"_1_trimmed.fastq \
#	in2="$OUTPUT"/quality_control/"$SHORTNAME"_2_trimmed.fastq \
#	ref="$OUTPUT"/SSU_sequences/"$SHORTNAME"_SSU_sequences.fasta \
#	scafstats="$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats.txt \
#	covstats="$OUTPUT"/SSU_sequences/"$SHORTNAME"_covstats.txt


#cat "$OUTPUT"/SSU_sequences/all_scafstats.txt | sed 1d | tr ',' \\t | awk '{ print $1,$8 }' | sort -k1 -k2 | awk '!($2==0){print}' | awk '{print $1}' > "$OUTPUT"/SSU_sequences/id_file.txt

#cat "$OUTPUT"/SSU_sequences/all_scafstats.txt | sed 1d | tr ',' \\t | awk '{ print $1,$8 }' | sort -k1 -k2 | awk '!($2==0){print}' > "$OUTPUT"/SSU_sequences/"$SHORTNAME"_scafstats.txt

#awk 'NR==FNR{ids[$0]; next} ($1 in ids){ printf ">" $0 }' "$OUTPUT"/SSU_sequences/id_file.txt RS='>' "$OUTPUT"/SSU_sequences/all_SSU_sequences.fasta > "$OUTPUT"/SSU_sequences/"$SHORTNAME"_SSU_sequences.fasta

#awk -F '\t' 'NR==FNR {id[$1]; next} $1 in id' "$OUTPUT"/SSU_sequences/id_file.txt "$OUTPUT"/SSU_sequences/all_covstats.txt > "$OUTPUT"/SSU_sequences/"$SHORTNAME"_covstats.txt

#rm "$OUTPUT"/SSU_sequences/all_scafstats.txt
#rm "$OUTPUT"/SSU_sequences/all_covstats.txt
#rm "$OUTPUT"/SSU_sequences/id_file.txt
#rm "$OUTPUT"/SSU_sequences/all_SSU_sequences.fasta


#rm "$OUTPUT"/SSU_sequences/emirge_metarib_clustered_SSU_sequences.fasta
#rm -r "$OUTPUT"/SSU_sequences/output_MetaRib/"$SHORTNAME"/Abundance

echo "Saving results..." | tee /dev/fd/3

echo ">Relative abundance calculation ends successfully on : "`date` | tee /dev/fd/3

conda deactivate
#echo "RiboTaxa virtual environment has been deactivated successfully..." | tee /dev/fd/3
