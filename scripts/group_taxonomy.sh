#!/usr/bin/env Rscript

#!/usr/bin/env Rscript

###########################################################################

# Post-RiboTaxa
# 31 October
# Author : Oshma Chakoory

###########################################################################
args <- commandArgs() # get arguments
PATH <- as.character(args[6])
OUTPUT <- as.character(args[7])
setwd(PATH)

#load library
library(plyr)
library(dplyr)
library(reshape2)
library(stringr)
#create empty dataframe
#Taxo_abundance <- data.frame(Taxonomy=character(),Abundance=numeric(), Study=character(), Sample=factor())
ComID_abundance <- data.frame(Taxonomy=character(),Abundance=numeric(), Study=character(), Sample=factor())
Species_abundance <- data.frame(Species=character(),Abundance=numeric(), Study=character(), Sample=factor())
Genus_abundance <- data.frame(Genus=character(),Abundance=numeric(), Study=character(), Sample=factor())
Family_abundance <- data.frame(Family=character(),Abundance=numeric(), Study=character(), Sample=factor())

#file in the format csv without header (sample,condition)
#example
#file.csv
#   ERR1776143,Control
#   ERR1776148,CD
#   ERR1776156,CD
df=read.csv(file='sample.csv', header = FALSE)
df$V1 <- as.character(df$V1)
df$V2 <- as.character(df$V2)

a <- list(strsplit(df$V1, ","))
b <- list(strsplit(df$V2, ","))

i=1
j=1
for (filename in a[[1]]) {
  file_name <- a[[1]][[i]]
  study <- b[[1]][[j]]
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
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  com_ID <- mutate(Taxo = paste(no_vertebrate_file$ID,no_vertebrate_file$Domain,no_vertebrate_file$Phylum,no_vertebrate_file$Class, no_vertebrate_file$Order, no_vertebrate_file$Family,no_vertebrate_file$Genus, no_vertebrate_file$Species, sep=";"),
                   no_vertebrate_file)
  com_ID <- com_ID[,c("Taxo","Relative_abundance")]
  com_ID_grouped <- com_ID %>% group_by(Taxo = as.factor(Taxo)) %>% summarise(Abundance = sum(Relative_abundance))
  Study_taxo <- as.data.frame(replicate(nrow(com_ID_grouped), study ))
  Sample_taxo <- as.data.frame(replicate(nrow(com_ID_grouped), file_name))
  
  final_ID <- cbind(com_ID_grouped, Study_taxo, Sample_taxo)
  colnames(final_ID) <-c("Taxonomy","Abundance","Study","Sample")
  ComID_abundance <- rbind(ComID_abundance,final_ID)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Getting Family, Genus, species table
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  family <- no_vertebrate_file[,c("Family","Relative_abundance")] #extract genus and %abundance
  family_grouped <- family %>% group_by(Family = as.factor(Family)) %>% summarise(Abundance = sum(Relative_abundance)) #group by species and sum up abundance of same species
  family_final <- subset(family_grouped, Family != "")  
  #family_final <- family_grouped %>% 
  #  filter(str_detect(Family, ''))
  Study_family <- as.data.frame(replicate(nrow(family_final), study))
  Sample_family <- as.data.frame(replicate(nrow(family_final), file_name))
  final_family <- cbind(family_final, Study_family, Sample_family)
  colnames(final_family) <-c("Family","Abundance","Study","Sample")
  Family_abundance <- rbind(Family_abundance,final_family)
  
  genus <- no_vertebrate_file[,c("Genus","Relative_abundance")] #extract genus and %abundance
  genus_grouped <- genus %>% group_by(Genus = as.factor(Genus)) %>% summarise(Abundance = sum(Relative_abundance)) #group by species and sum up abundance of same species
  genus_final <- subset(genus_grouped, Genus != "")  
  #genus_final <- genus_grouped %>% 
  #  filter(str_detect(Genus, ''))
  Study_Ge <- as.data.frame(replicate(nrow(genus_final), study))
  Sample_Ge <- as.data.frame(replicate(nrow(genus_final), file_name))
  final_Ge <- cbind(genus_final, Study_Ge, Sample_Ge)
  colnames(final_Ge) <-c("Genus","Abundance","Study","Sample")
  Genus_abundance <- rbind(Genus_abundance,final_Ge)
  
  species <- no_vertebrate_file[,c("Species","Relative_abundance")] #extract species and %abundance
  species_grouped <- species %>% group_by(Species = as.factor(Species)) %>% summarise(Abundance = sum(Relative_abundance)) #group by species and sum up abundance of same species
  species_final <- subset(species_grouped, Species != "")  
  #species_final <- species_grouped %>% 
  #  filter(str_detect(Species, ''))
  Study_Sp <- as.data.frame(replicate(nrow(species_final), study))
  Sample_Sp <- as.data.frame(replicate(nrow(species_final), file_name))
  final_Sp <- cbind(species_final, Study_Sp, Sample_Sp)
  colnames(final_Sp) <-c("Species","Abundance","Study","Sample")
  Species_abundance <- rbind(Species_abundance,final_Sp)
  
  i=i+1
  j=j+1
}

Genus_abundance$Genus <- str_replace_all(Genus_abundance$Genus, "g__", "")
Species_abundance$Species <- str_replace_all(Species_abundance$Species, "s__", "")
Family_abundance$Family <- str_replace_all(Family_abundance$Family, "f__", "")

setwd(OUTPUT)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Convert long table to wide table 
# row represent microorganism
# column represent sample
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#convert long table to wide table for whole Taxonomy
#long_taxo <- read.table(file = 'ComID_Taxonomy_abundance.csv', sep = ',', header = TRUE)
long_taxo <- arrange(ComID_abundance,Study) #arrange by study
wide_table <- reshape(long_taxo, idvar = c("Study","Sample"), timevar = c("Taxonomy"), direction = "wide") 
colnames(wide_table) = gsub("Abundance.", "", colnames(wide_table))
wide_table <- as.data.frame(t(wide_table),stringsAsFactors=FALSE)
wide_table[is.na(wide_table)] <- 0
write.table(wide_table,"Complete_Taxonomy_abundance.csv", row.names = TRUE, col.names = FALSE, sep = "\t")

#convert long table to wide table for genus
#long_taxo <- read.table(file = 'Family_Taxonomy_abundance.csv', sep = ',', header = TRUE)
long_taxo <- arrange(Family_abundance,Study) #arrange by study
wide_table <- reshape(long_taxo, idvar = c("Study","Sample"), timevar = c("Family"), direction = "wide") 
colnames(wide_table) = gsub("Abundance.", "", colnames(wide_table))
wide_table <- as.data.frame(t(wide_table),stringsAsFactors=FALSE)
wide_table[is.na(wide_table)] <- 0
write.table(wide_table,"Family_Taxonomy_abundance.csv", row.names = TRUE, col.names = FALSE, sep = "\t")

#convert long table to wide table for genus
#long_taxo <- read.table(file = 'Genus_Taxonomy_abundance.csv', sep = ',', header = TRUE)
long_taxo <- arrange(Genus_abundance,Study) #arrange by study
wide_table <- reshape(long_taxo, idvar = c("Study","Sample"), timevar = c("Genus"), direction = "wide") 
colnames(wide_table) = gsub("Abundance.", "", colnames(wide_table))
wide_table <- as.data.frame(t(wide_table),stringsAsFactors=FALSE)
wide_table[is.na(wide_table)] <- 0
write.table(wide_table,"Genus_Taxonomy_abundance.csv", row.names = TRUE, col.names = FALSE, sep = "\t")

#convert long table to wide table for Species
#long_taxo <- read.table(file = 'Species_Taxonomy_abundance.csv', sep = ',', header = TRUE)
long_taxo <- arrange(Species_abundance,Study) #arrange by study
wide_table <- reshape(long_taxo, idvar = c("Study","Sample"), timevar = c("Species"), direction = "wide") 
colnames(wide_table) = gsub("Abundance.", "", colnames(wide_table))
wide_table <- as.data.frame(t(wide_table),stringsAsFactors=FALSE)
wide_table[is.na(wide_table)] <- 0
write.table(wide_table,"Species_Taxonomy_abundance.csv", row.names = TRUE, col.names = FALSE, sep = "\t")
