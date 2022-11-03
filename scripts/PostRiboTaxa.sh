#!/usr/bin/env Rscript

###########################################################################

# Post-RiboTaxa
# 31 October
# Author : Oshma Chakoory

###########################################################################

#activate virtual environment for RiboTaxa
#source activate RiboTaxa_py36
#echo "Qiime2 virtual environment has been activated successfully..." | tee /dev/fd/3

args <- commandArgs() # get arguments

#load library
library(plyr)
library(dplyr)
library(stringr)

file_name <- as.character(args[6])

file <- read.table(file = paste0(file_name,"_SSU_taxonomy_abundance.tsv"), sep = '\t', header=TRUE, fill=TRUE)
file$X <- NULL
file1 <- arrange(file,Species,Genus) #arrange by species and genus
no_vertebrate_file <- file1 %>% 
  filter(!str_detect(Phylum, 'Vertebrata')) #discard all vertebrate
no_vertebrate_file <- as_tibble(no_vertebrate_file)
no_vertebrate_file[12] <- apply(no_vertebrate_file[,12],2,function(x){x/sum(x)*100}) #recalculate relative abundance


no_vertebrate_file <- no_vertebrate_file %>% 
  mutate(ID = str_replace(ID, "^[0-9|]+", ""))

no_vertebrate_file$Domain <- as.factor(no_vertebrate_file$Domain)
no_vertebrate_file$Phylum <- as.factor(no_vertebrate_file$Phylum)
no_vertebrate_file$Class <- as.factor(no_vertebrate_file$Class)
no_vertebrate_file$Order <- as.factor(no_vertebrate_file$Order)
no_vertebrate_file$Family <- as.factor(no_vertebrate_file$Family)
no_vertebrate_file$Genus <- as.factor(no_vertebrate_file$Genus)
no_vertebrate_file$Species <- as.factor(no_vertebrate_file$Species)

#------------------------------
#working on completing taxonomy
#------------------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#working on uncultured species
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#replace all uncultured by uncultured previous rank
no_vertebrate_file$Phylum[no_vertebrate_file$Phylum=='p__uncultured'] <- NA #replace all uncultured_bacterium by NA
no_vertebrate_file <- no_vertebrate_file %>%
  mutate(Phylum = coalesce(Phylum,Domain)) #copy Domain to replace NA in Phylum column

#class
no_vertebrate_file$Class[no_vertebrate_file$Class=='c__uncultured'] <- NA #replace all uncultured_bacterium by NA
no_vertebrate_file <- no_vertebrate_file %>%
  mutate(Class = coalesce(Class,Phylum)) #copy Phylum to replace NA in Class column

#order
no_vertebrate_file$Order[no_vertebrate_file$Order=='o__uncultured'] <- NA #replace all uncultured_bacterium by NA
no_vertebrate_file <- no_vertebrate_file %>%
  mutate(Order = coalesce(Order,Class)) #copy Class to replace NA in Order column

#family
no_vertebrate_file$Family[no_vertebrate_file$Family=='f__uncultured'] <- NA #replace all uncultured_bacterium by NA
no_vertebrate_file <- no_vertebrate_file %>%
  mutate(Family = coalesce(Family,Order)) #copy Order to replace NA in Family column

#genus
no_vertebrate_file$Genus[no_vertebrate_file$Genus=='g__uncultured'] <- NA #replace all uncultured_bacterium by NA
no_vertebrate_file <- no_vertebrate_file %>%
  mutate(Genus = coalesce(Genus,Family)) #copy Family to replace NA in genus column

#levels(no_empty_species$Species) <- sub('g__', 's__unclassified_', levels(no_empty_species$Species)) #replace g__ by s__uncultured_
no_vertebrate_file$Species[no_vertebrate_file$Species=='s__uncultured_bacterium'] <- NA #replace all uncultured_bacterium by NA
no_vertebrate_file$Species[no_vertebrate_file$Species=='s__bacterium_New'] <- NA #replace all bacterium_New by NA
no_vertebrate_file$Species[no_vertebrate_file$Species=='s__gut_metagenome'] <- NA
no_vertebrate_file$Species[no_vertebrate_file$Species=='s__human_gut'] <- NA #replace all gut_metagenome by NA
no_vertebrate_file$Species[no_vertebrate_file$Species=='s__uncultured_organism'] <- NA #replace all uncultured_organism by NA
no_vertebrate_file$Species[no_vertebrate_file$Species=='s__metagenome'] <- NA #replace all metagenome by NA
no_vertebrate_file$Species[no_vertebrate_file$Species=='s__unidentified'] <- NA #replace all unidentified by NA
no_vertebrate_file$Species[no_vertebrate_file$Species=='s__uncultured_delta'] <- NA #replace all unidentified by NA
no_vertebrate_file$Species[no_vertebrate_file$Species=='s__marine_metagenome'] <- NA #replace all unidentified by NA
no_vertebrate_file$Species[no_vertebrate_file$Species=='s__uncultured_gamma'] <- NA #replace all unidentified by NA
no_vertebrate_file$Species[no_vertebrate_file$Species=='s__uncultured_soil'] <- NA #replace all unidentified by NA
no_vertebrate_file$Species[no_vertebrate_file$Species=='s__uncultured_marine'] <- NA #replace all unidentified by NA
no_vertebrate_file$Species[no_vertebrate_file$Species=='s__marine_gamma'] <- NA #replace all unidentified by NA
no_vertebrate_file$Species[no_vertebrate_file$Species=='s__uncultured_rumen'] <- NA #replace all unidentified by NA
no_vertebrate_file$Species[no_vertebrate_file$Species=='s__swine_fecal'] <- NA #replace all unidentified by NA


no_vertebrate_file <- no_vertebrate_file %>%
  mutate(Species = coalesce(Species,Genus)) #copy genus to empty spaces in species column
#levels(no_vertebrate_file$Species) <- sub('g__', 's__uncultured_', levels(no_vertebrate_file$Species)) #replace g__ by s__uncultured_

levels(no_vertebrate_file$Phylum) <- sub('d__', 'p__uncultured_', levels(no_vertebrate_file$Phylum)) #replace f__ by g__uncultured_

levels(no_vertebrate_file$Class) <- sub('d__', 'c__uncultured_', levels(no_vertebrate_file$Class)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Class) <- sub('p__', 'c__uncultured_', levels(no_vertebrate_file$Class)) #replace f__ by g__uncultured_

levels(no_vertebrate_file$Order) <- sub('d__', 'o__uncultured_', levels(no_vertebrate_file$Order)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Order) <- sub('p__', 'o__uncultured_', levels(no_vertebrate_file$Order)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Order) <- sub('c__', 'o__uncultured_', levels(no_vertebrate_file$Order)) #replace f__ by g__uncultured_

levels(no_vertebrate_file$Family) <- sub('d__', 'f__uncultured_', levels(no_vertebrate_file$Family)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Family) <- sub('p__', 'f__uncultured_', levels(no_vertebrate_file$Family)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Family) <- sub('c__', 'f__uncultured_', levels(no_vertebrate_file$Family)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Family) <- sub('o__', 'f__uncultured_', levels(no_vertebrate_file$Family)) #replace f__ by g__uncultured_

levels(no_vertebrate_file$Genus) <- sub('d__', 'g__uncultured_', levels(no_vertebrate_file$Genus)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Genus) <- sub('p__', 'g__uncultured_', levels(no_vertebrate_file$Genus)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Genus) <- sub('c__', 'g__uncultured_', levels(no_vertebrate_file$Genus)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Genus) <- sub('o__', 'g__uncultured_', levels(no_vertebrate_file$Genus)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Genus) <- sub('f__', 'g__uncultured_', levels(no_vertebrate_file$Genus)) #replace f__ by g__uncultured_

levels(no_vertebrate_file$Species) <- sub('d__', 's__uncultured_', levels(no_vertebrate_file$Species)) #replace g__ by s__uncultured_
levels(no_vertebrate_file$Species) <- sub('p__', 's__uncultured_', levels(no_vertebrate_file$Species)) #replace g__ by s__uncultured_
levels(no_vertebrate_file$Species) <- sub('c__', 's__uncultured_', levels(no_vertebrate_file$Species)) #replace g__ by s__uncultured_
levels(no_vertebrate_file$Species) <- sub('o__', 's__uncultured_', levels(no_vertebrate_file$Species)) #replace g__ by s__uncultured_
levels(no_vertebrate_file$Species) <- sub('f__', 's__uncultured_', levels(no_vertebrate_file$Species)) #replace g__ by s__uncultured_
levels(no_vertebrate_file$Species) <- sub('g__', 's__uncultured_', levels(no_vertebrate_file$Species)) #replace g__ by s__uncultured_

#replace all empty spaces by NA
no_vertebrate_file$Phylum[no_vertebrate_file$Phylum==''] <- NA
no_vertebrate_file <- no_vertebrate_file %>%
  mutate(Phylum = coalesce(Phylum,Domain)) #copy Domain to replace NA in Phylum column

no_vertebrate_file$Class[no_vertebrate_file$Class==''] <- NA
no_vertebrate_file <- no_vertebrate_file %>%
  mutate(Class = coalesce(Class,Phylum)) #copy Phylum to replace NA in Class column

no_vertebrate_file$Order[no_vertebrate_file$Order==''] <- NA
no_vertebrate_file <- no_vertebrate_file %>%
  mutate(Order = coalesce(Order,Class)) #copy Class to replace NA in Order column

no_vertebrate_file$Family[no_vertebrate_file$Family==''] <- NA
no_vertebrate_file <- no_vertebrate_file %>%
  mutate(Family = coalesce(Family,Order)) #copy Order to replace NA in Family column

no_vertebrate_file$Genus[no_vertebrate_file$Genus==''] <- NA
no_vertebrate_file <- no_vertebrate_file %>%
  mutate(Genus = coalesce(Genus,Family)) #copy Family to replace NA in genus column

no_vertebrate_file$Species[no_vertebrate_file$Species==''] <- NA #replace all empty spaces by NA
no_vertebrate_file <- no_vertebrate_file %>%
  mutate(Species = coalesce(Species,Genus)) #copy genus to replace NA in species column

levels(no_vertebrate_file$Phylum) <- sub('d__', 'p__unclassified_', levels(no_vertebrate_file$Phylum)) #replace f__ by g__uncultured_

levels(no_vertebrate_file$Class) <- sub('d__', 'c__unclassified_', levels(no_vertebrate_file$Class)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Class) <- sub('p__', 'c__unclassified_', levels(no_vertebrate_file$Class)) #replace f__ by g__uncultured_

levels(no_vertebrate_file$Order) <- sub('d__', 'o__unclassified_', levels(no_vertebrate_file$Order)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Order) <- sub('p__', 'o__unclassified_', levels(no_vertebrate_file$Order)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Order) <- sub('c__', 'o__unclassified_', levels(no_vertebrate_file$Order)) #replace f__ by g__uncultured_

levels(no_vertebrate_file$Family) <- sub('d__', 'f__unclassified_', levels(no_vertebrate_file$Family)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Family) <- sub('p__', 'f__unclassified_', levels(no_vertebrate_file$Family)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Family) <- sub('c__', 'f__unclassified_', levels(no_vertebrate_file$Family)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Family) <- sub('o__', 'f__unclassified_', levels(no_vertebrate_file$Family)) #replace f__ by g__uncultured_

levels(no_vertebrate_file$Genus) <- sub('d__', 'g__unclassified_', levels(no_vertebrate_file$Genus)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Genus) <- sub('p__', 'g__unclassified_', levels(no_vertebrate_file$Genus)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Genus) <- sub('c__', 'g__unclassified_', levels(no_vertebrate_file$Genus)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Genus) <- sub('o__', 'g__unclassified_', levels(no_vertebrate_file$Genus)) #replace f__ by g__uncultured_
levels(no_vertebrate_file$Genus) <- sub('f__', 'g__unclassified_', levels(no_vertebrate_file$Genus)) #replace f__ by g__uncultured_

levels(no_vertebrate_file$Species) <- sub('d__', 's__unclassified_', levels(no_vertebrate_file$Species)) #replace g__ by s__uncultured_
levels(no_vertebrate_file$Species) <- sub('p__', 's__unclassified_', levels(no_vertebrate_file$Species)) #replace g__ by s__uncultured_
levels(no_vertebrate_file$Species) <- sub('c__', 's__unclassified_', levels(no_vertebrate_file$Species)) #replace g__ by s__uncultured_
levels(no_vertebrate_file$Species) <- sub('o__', 's__unclassified_', levels(no_vertebrate_file$Species)) #replace g__ by s__uncultured_
levels(no_vertebrate_file$Species) <- sub('f__', 's__unclassified_', levels(no_vertebrate_file$Species)) #replace g__ by s__uncultured_
levels(no_vertebrate_file$Species) <- sub('g__', 's__unclassified_', levels(no_vertebrate_file$Species)) #replace g__ by s__uncultured_
levels(no_vertebrate_file$Species) <- sub('s__unclassified_uncultured_', 's__unclassified_', levels(no_vertebrate_file$Species)) #replace g__ by s__uncultured_

#save results in csv files
write.table(no_vertebrate_file, file = paste0(file_name,"_SSU_taxonomy_abundance.tsv"), row.names = FALSE, sep = "\t")
