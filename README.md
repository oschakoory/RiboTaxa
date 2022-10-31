<img align="right" src="docs/RiboaTaxa_Icon.png" width="200" alt="RiboTaxa logo"/>

# RiboTaxa v1.4

by Oshma Chakoory, Sophie Marre, and Pierre Peyret.

RiboTaxa is a complete pipeline to rapidly filter and reconstruct the full length SSU rRNA gene from Illumina (meta)genomic dataset and perform taxonomic classification on the reconstructed sequences.

RiboTaxa takes as input singled-end or paired-end files which can be in compressed format (fastq.gz) or uncompressed format (.fastq).

Tools used in RiboTaxa pipeline:
- For quality control :<a class="reference external" href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/" target="_blank" rel="noopener noreferrer">FastQC</a>, <a class="reference external" href="https://multiqc.info/" target="_blank" rel="noopener noreferrer">MultiQC</a>
- For adapters removal and trimming: <a class="reference external" href="https://jgi.doe.gov/data-and-tools/bbtools/" target="_blank" rel="noopener noreferrer">BBTOOLS</a>
- To filter 16S/18S reads: <a class="reference external" href="https://academic.oup.com/bioinformatics/article/28/24/3211/246053" target="_blank" rel="noopener noreferrer">SortMeRNA</a> 
- To reconstruct full-length SSU rRNA sequences: <a class="reference external" href="https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-5-r44" target="_blank" rel="noopener noreferrer">EMIRGE</a>, <a class="reference external" href="https://pubmed.ncbi.nlm.nih.gov/32167532/" target="_blank" rel="noopener noreferrer">MetaRib</a> 
- Classify the full-length reconstructed SSU sequences: <a class="reference external" href="https://docs.qiime2.org/2020.8/plugins/available/feature-classifier/classify-sklearn/">sklearn classifier</a> of <a class="reference external" href="https://docs.qiime2.org/2020.8/" target="_blank" rel="noopener noreferrer">QIIME2</a>

<img align="center" src="docs/RiboTaxa.png" width="700" alt="RiboTaxa Pipeline"/>

## Quick-start

### Install Miniconda

Miniconda provides the conda environment and package manager, and is the recommended way to install RiboTaxa. Follow the <a class="reference external" href="https://conda.io/projects/conda/en/latest/user-guide/install/index.html" target="_blank" rel="noopener noreferrer">Miniconda instructions</a> for downloading and installing Miniconda. It is important to follow all of the directions provided in the <a class="reference external" href="https://conda.io/projects/conda/en/latest/user-guide/install/index.html" target="_blank" rel="noopener noreferrer">Miniconda instructions</a>, particularly ensuring that you run the following commands at the end of the installation process, to ensure that your Miniconda installation is fully installed and configured correctly:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

```

### Updating Miniconda

After installing Miniconda and opening a new terminal, make sure you’re running the latest version of conda:

```bash
conda update conda
```

### Install RiboTaxa within a conda environment

Once you have Miniconda installed, run the following command to create the required conda environments for RiboTaxa. Creating a virtual environment is recommended to avoid conflicts between required dependencies and those in existing environment. 

```bash
git clone https://github.com/oschakoory/RiboTaxa.git
cd RiboTaxa
bash conda_virt_env.sh
```
The first virtual environment ```RiboTaxa_py27``` uses python 2.7 and contains all the necessary tools for quality control, to filter 16S18S reads and reconstruct full length rRNAs sequences.

The second virtual environment ```RiboTaxa_py36``` uses python 3.6 and runs sklearn classifier of QIIME2 to taxonomically classify full length reconstruted rRNAs sequences.

### Activate the conda environment

Activation of virtual environment is done automatically in scripts. No need to activate environment manually, unless you want to scrutinize each tool used in RiboTaxa.

```bash
conda activate RiboTaxa_py27
conda list

emirge_amplicon.py --help
```
```bash
conda activate RiboTaxa_py36
conda list
```

### Install usearch (if you do not have usearch in your $path)

USEARCH is a unique sequence analysis tool which offers search and clustering algorithms that are often orders of magnitude faster than BLAST. To reconstruct full length rRNA sequences, ```EMIRGE``` uses ```USEARCH```. 

Please ensure that you have ```USEARCH``` installed in your ```$PATH``` using the following command:

```bash
usearch --version
echo $PATH
```

If not, please follow the <a class="reference external" href="https://pbertinblog.wordpress.com/usearch-installation/" target="_blank" rel="noopener noreferrer">usearch-installation instructions</a> to download ```USEARCH``` and install it in your ```$PATH```.

To make sure that usearch is successfully installed, run

```bash
usearch --version
```

Now, you are ready to use RiboTaxa !!!

## Let the power of RiboTaxa begin...

### Indexing databases

RiboTaxa pipeline includes tools like sortmerna and emirge, both of which need indexed databases of their own. The latest database SILVA SSU 138.1 can be downloaded <a class="reference external" href="https://www.arb-silva.de/" target="_blank" rel="noopener noreferrer">here</a> and the indexed database SILVA SSU 138.1 for RiboTaxa are available <a class="reference external" href="https://ucafr-my.sharepoint.com/:f:/g/personal/oshma_chakoory_uca_fr/EkqL_9bA2UpApbEHiB8cfJ4B4BBsMB25O9Xtl29wi2QSIw?e=VqqWM9" target="_blank" rel="noopener noreferrer">here</a>. 

To index your own database, you will need to fill the config file ```indexDB_arguments.conf```.If you are not sure of certains parameters, leave as defined except for directories and input files.

```bash
The configuration file is very important and each parameter needs to be filled to avoid errors.

script used: indexDB_RiboTaxa.sh

[Setting up directories...]

####set up RiboTaxa directory path **
RiboTaxa_DIR = /home/user/Documents/RiboTaxa

[Setting up database path...]

####set up database path+database file in fasta format **
DB_DIR = /home/user/Documents/Databases/SILVA_138_SSURef_Nr99_tax_silva.fasta

#### Set up output directory **
OUTPUT = /home/user/Documents/Databases

#### set up the number of threads/CPUs to be used for indexing
THREAD = 8

#### set up the clustering id threshold (between 0.7-1)
CLUSTER_ID = 0.97
```

Once it is filled with all the necessary information, you can run the following command and index your database.

```bash
bash -i RiboTaxa_DIR/indexDB_RiboTaxa.sh PATH_TO/indexDB_arguments.conf
```

Indexing database takes a while. Using the maximum number of available threads/CPUs will save time. This step will produce two directories in your ```OUTPUT``` path:
- sortmerna_indexed_DB : containing indexed files for sortmeRNA
- bowtie_indexed_DB : containing indexed files by <a class="reference external" href="http://bowtie-bio.sourceforge.net/manual.shtml" target="_blank" rel="noopener noreferrer">Bowtie</a> to be used for EMIRGE and MetaRib

For the taxonomic classification by sklearn classifier, the database used is the trained classifier ```SILVA 138 reference sequence``` downloaded from the <a class="reference external" href="https://docs.qiime2.org/2020.8/data-resources/" target="_blank" rel="noopener noreferrer">Data resources</a> of <a class="reference external" href="https://docs.qiime2.org/2020.8/" target="_blank" rel="noopener noreferrer">Qiime2</a>. To use other databases such as Greengenes or UNITE, you can download already trained classifer or train your own database by following the <a class="reference external" href="https://docs.qiime2.org/2020.8/data-resources/" target="_blank" rel="noopener noreferrer">Data resources</a> instructions.

### Running RiboTaxa pipeline

RiboTaxa pipeline will 
- Remove adapters and trim (meta)genomics data using <a class="reference external" href="https://jgi.doe.gov/data-and-tools/bbtools/" target="_blank" rel="noopener noreferrer">BBTOOLS</a>
- Filter 16S/18S reads using <a class="reference external" href="https://academic.oup.com/bioinformatics/article/28/24/3211/246053" target="_blank" rel="noopener noreferrer">SortMeRNA</a> 
- Reconstruct full-length SSU rRNA sequences using <a class="reference external" href="https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-5-r44" target="_blank" rel="noopener noreferrer">EMIRGE</a> and <a class="reference external" href="https://pubmed.ncbi.nlm.nih.gov/32167532/" target="_blank" rel="noopener noreferrer">MetaRib</a>
- Classify the full-length reconstructed SSU sequences using <a class="reference external" href="https://docs.qiime2.org/2020.8/plugins/available/feature-classifier/classify-sklearn/">sklearn classifier</a> of <a class="reference external" href="https://docs.qiime2.org/2020.8/" target="_blank" rel="noopener noreferrer">QIIME2</a>


RiboTaxa can be used for one singled-end/paired-end dataset or multiple singled-end/paired-end datasets in the same folder.

To run RiboTaxa, you will need to fill the config file ```RiboTaxa_arguments.conf```. If you are not sure of certains parameters, leave as defined except for fields denoted by **.

```bash
##The configuration file is very important and each parameter needs to be filled to avoid errors.
##Mandatory parameters are denoted by **. Most of them depend on the configuration of your computer capacities/configuration and of your sequencing data.
Other tool parameters have been optimized and can be left at default. 

#-----------------------
#Setting up directories
#-----------------------

[BASE]
####set up RiboTaxa directory **
RiboTaxa_DIR = /home/user/RiboTaxa

####set up data directory containing only raw reads in fastq/fastq.gz format **
###paired-end files should be metagenome_R1.fastq/fastq.gz and metagenome_R2.fastq/fastq.gz
###singled-end file should be metagenome.fastq/fastq.gz
DATA_DIR = /home/user/Documents/raw_reads

####format of your paired end files **
## fastq : if files are not compressed
## fastq.gz: if files are compressed in gz format
FORMAT = fastq

####set up output directory **
OUTPUT = /home/user/Documents/RiboTaxa_results

####number of threads/CPUS to be used through the pipeline **
THREAD = 8

[BBMAP]
####RAM limit to be used by BBTOOLS/BBMAP during quality control and mapping **
##depends of the computer RAM. Use approx 80% of available RAM
##example: Available RAM = 16GB, therefore RAM = 80/100*16 = 12GB
RAM = 12

#-----------------------------
#Quality control using BBTOOLS
#-----------------------------

####Trim reads to remove bases matching adapter sequences (Default value = r)
## f (don't trim)
## r (trim to the right)
## l (trim to the left)
##In ktrim=r mode, once a reference kmer is matched in a read, 
##that kmer and all the bases to the right will be trimmed, leaving 
##only the bases to the left; this is the normal mode for adapter trimming.

ktrim = r

####Kmer length used for finding adapters (Default value = 21)
kmer = 21

####Reads shorter than this length (bases) after trimming will be discarded (Default value = 60)
minlength = 60

#Regions with average quality BELOW this will be trimmed (Default value = 20)
trimq = 20

####Trim read ends to remove bases with quality below trimq (Default value = rl)
# rl (trim both ends), 
# f (neither end), 
# r (right end only), 
# l (left end only),
# w (sliding window)

qtrim = rl

####reads with more Ns than this (after trimming) will be discarded (Default value = 1)
maxns = 1

#------------------------------------
#Filter 16S/18S reads using SortmeRNA
#------------------------------------

####indexed database directory for sortmerna **
##This directory should contain one .clustered.fasta and several .clustered. files
SORTMERNA_DB = /home/user/Documents/Databases/sortmerna_indexed_DB


#----------------------------------------------------------
#Reconstructing 16S/18S sequences using EMIRGE and MetaRIB
#----------------------------------------------------------

[EMIRGE]

####set up database directory containing indexed files for emirge **
##This directory should contain one fasta file and several .ebwt files
EMIRGE_DB = /home/user/Documents/Databases/bowtie_indexed_DB

####length of longest reads **
MAX_LENGTH = 300

####identity threshold (Default value = 1)
##This the JOIN_TRESHOLD parameter of EMIRGE. If two candidate sequences 
##share >= this fractional identity over their bases with mapped reads, then
##merge the two sequences into one for the next iteration. Fixed to 1, sequence
##reconstruction gave the best results on controlled samples (mock or synthetic communities).
IDENTITY = 1 

####number of iterations (Default value = 40)
##Number of iterations to perform by EMIRGE for sequence reconstruction.
##The default value can fit to most of metagenomic data. EMIRGE authors recommended
to increase this value for complex communities. 
NUM_ITERATION = 40

####mean insert size **
##Insert size distribution mean
MEAN_INSERT_SIZE = 300 

####standard deviation **
##Insert size distribution standard deviation
STD_DEV = 100 

####minimum fraction of the length of a candidate
##reference sequence that must be covered by mapped
##reads (Default=0.3) Range [0.0,1.0]
MIN_COV = 0.3

[METARIB]
####Subsampling reads number in each iteration (Default=1000000)
SAMPLING_NUM = 1000000

#-----------------------------------------------------------
#Taxonomic classfication using sklearn_classifier of qiime2
#-----------------------------------------------------------

####Set up path+database name for sklearn classifier **
SKLEARN_DB = /home/user/Documents/Databases/qiime2020.8_silva138/silva-138-99-nb-classifier.qza

####Confidence threshold for limiting taxonomic depth (default = 0.7)
##This threshold ensures qualititative affiliation of reconstructed sequences. 
##We do not recommend to lower this value under 0.7. 
CONFIDENCE = 0.7

####Number of reads to process in each batch (default = 0)
##Use BATCH = 1 if you have less than 16GB RAM to avoid errors
BATCH = 0
```

Once it is filled with all the necessary information, you can use the following command to run the RiboTaxa pipeline for

##### Singled-end dataset(s)

```bash
bash -i RiboTaxa_DIR/Pipeline_RiboTaxa_SE.sh PATH_TO/RiboTaxa_arguments.conf
```

##### Paired-end dataset(s)

```bash
bash -i RiboTaxa_DIR/Pipeline_RiboTaxa_PE.sh PATH_TO/RiboTaxa_arguments.conf
```

For each singled-end/paired-end sample, RiboTaxa will create one directory using the sample name in your OUTPUT path of your ```RiboTaxa_arguments.conf``` file. Each sample directory will contain the 4 following sub-directories:
- ```quality_control``` : This directory contains your (meta)genomics files after adpaters removal and trimming. It also has two sub_directories ```before_fastqc``` and ```after_fastqc``` containing quality reports of your sequence files before and after trimming. You may look at the ```.html``` files in each sub-directory to have an overview of each (meta)genomics file or look into ```multiqc``` folder to have an overview of all the (meta)genomics files given to this pipeline.

- ```output_sortmerna``` : This folder contains filtered 16S/18S sequences from your trimmed (meta)genomics sequence files: ```metagenome_R1_16S18Sreads.fastq```,```metagenome_R2_16S18Sreads.fastq```(for paired-end) or ```metagenome_16S18Sreads.fastq``` (for singled-end) and a ```metagenome.log``` indicating the % of 16S/18S reads filtered from your (meta)genomics dataset.

-  ```SSU_sequences``` : This folder contains two subfolders: ```output_emirge``` and ```output_metarib```. These subfolders contain iterations perfomed by EMIRGE and MetaRib to reconstruct full-length/nearly full-length rRNA 16S/18S gene sequences. 

- The file ```metagenome_SSU_sequences.fasta``` contains the final full-length/nearly full-length rRNA 16S/18S gene sequences recontructed by EMIRGE and MetaRib.
	
- Taxonomy : This folder contains the taxonomic classification of the full-length SSU rRNA sequences reconstructed by EMIRGE. It takes ```metagenome_SSU_sequences.fasta``` as input, converts it to ```metagenome_SSU_sequences_qiime2.qza``` and returns the classification of each sequence as ```metagenome_renamed_16S18S_recons_qiime2_taxonomy.qza```. To calculate relative abundace of each reconstructed SSU sequence, BBmap is used to map short reads onto SSU sequence and the number of assigned reads per sequence is divided by the sum of all the assigned reads (and multiplied by 100). Relative abundance is thus expressed in %. To view the classification and relative abundance of each sequence, open ```metagenome_SSU_taxonomy_abundance.tsv```.

. The ```metagenome_SSU_taxonomy_abundance.tsv``` contains the following column names:

```bash
Sequence_ID		Domain		Phylum			Class			Order		Family		Genus			Species		Confidence 	Length(bp)	Assigned reads		Relative_Abundance(%)
3|EU334524.1.1558	Bacteria	Desulfobacterota	Desulfuromonadia	Geobacterales	Geobacteraceae	Trichlorobacter	Geobacter_lovleyi	0.999874566	1425			68			4.6687
```


### Running RiboTaxa pipeline on test data

To run RiboTaxa pipeline on test data, you need to download the indexed database of SILVA SSU 138.1 available <a class="reference external" href="https://ucafr-my.sharepoint.com/:f:/g/personal/oshma_chakoory_uca_fr/EkqL_9bA2UpApbEHiB8cfJ4B4BBsMB25O9Xtl29wi2QSIw?e=VqqWM9" target="_blank" rel="noopener noreferrer">here</a>

For the taxonomic classification by sklearn classifier, the database used is the trained classifier ```SILVA 138 reference sequence``` downloaded from the <a class="reference external" href="https://docs.qiime2.org/2020.8/data-resources/" target="_blank" rel="noopener noreferrer">Data resources</a> of <a class="reference external" href="https://docs.qiime2.org/2020.8/" target="_blank" rel="noopener noreferrer">Qiime2</a>.

Once all the indexed databases are successfully downloaded into ```~/database```, you need to update the directory path of each in your ```RiboTaxa_arguments.conf```. The remaining parameters can be left as default.

Then run: 

```bash
bash -i RiboTaxa_DIR/Pipeline_RiboTaxa_PE.sh RiboTaxa_DIR/test_data/RiboTaxa_arguments.conf
```

The test data directory also contains the results of RiboTaxa of the test sample.


### Group multiple samples taxonomy files into one

To group taxonomic files of multiple samples into one, you will need to copy all taxonomic files (```*_SSU_taxonomy_abundance.tsv```) into the same folder along with a ```sample.csv``` (see in ```RiboTaxa_DIR/test_data/multiple_samples_Taxonomy```):

```bash
Sample01,Control
Sample02,Control
Sample03,Control
Sample04,Treated
Sample05,Treated
Sample06,Treated
```

Then run:

```bash
cd RiboTaxa_DIR/scripts

./RiboTaxa_group_taxonomy.sh $INPUT_PATH $OUTPUT_PATH
```
whereby:
```$INPUT_PATH``` is the input path containing the files
```$OUTPUT_PATH``` is the desired output path

An exemple test has been conducted with the following parameters:

```$INPUT_PATH``` = ```RiboTaxa_DIR/test_data/multiple_samples_Taxonomy```

```$OUTPUT_PATH``` = ```RiboTaxa_DIR/test_data/multiple_samples_Taxonomy/Output```

The results are in ```RiboTaxa_DIR/test_data/multiple_samples_Taxonomy/Output```.



In the output folder, there are 4 files:
1. ```Complete_taxonomy_abundance.csv```: containing the abundance of the complete taxonomy in all the samples
2. ```Family_abundance.csv```: containing the abundance of the families in all the samples
3. ```Genus_abundance.csv```: containing the abundance of the genera in all the samples
4. ```Species_abundance.csv```: containing the abundance of the species in all the samples

Example of the ```Species_abundance.csv```:

| Study                     | Control  | Control     | Control     | Treated     | Treated     | Treated     |
|---------------------------|----------|-------------|-------------|-------------|-------------|-------------|
| Sample                    | Sample01 | Sample02    | Sample03    | Sample04    | Sample05    | Sample06    |
| Bifidobacterium_bifidum   | 0        |  0.2050580  |  0.4260770  | 12.6970194  |  7.3577775  |  0.5264176  |
| Bifidobacterium_sp.       | 0        | 0.008041001 | 0.098325006 | 0.007881000 | 0.045245990 | 0.642682495 |
| Cosenzaea_myxofaciens     | 0        | 0.04422800  | 0           | 0           | 0           | 0.03552537  |
| Enterobacter_sp.          | 0        | 1.881710    | 0.732972    | 0.165496    | 0.010442    | 0           |
| Lactobacillus_acidophilus | 0        |  2.850710   |  7.460818   |  0.058543   |  4.284497   | 12.540370   |
| Lactobacillus_kitasatonis | 0        | 0.0241240   | 0.5095060   | 0           | 0.8405428   | 1.8279293   |
| Proteus_mirabilis         | 0        |  0.05629001 | 18.44943811 |  6.40029719 | 14.55544994 | 42.99186000 |
| Proteus_sp.               | 0        | 0.036187    | 0           | 0           | 0           | 0           |
| Streptomyces_cinnamoneus  | 0        | 0.06433201  | 0           | 0           | 0.01044200  | 0           |
| uncultured_Syntrophomonas | 0        | 0.016083    | 0           | 0           | 0           | 0           |
