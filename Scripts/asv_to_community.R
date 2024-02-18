#!/usr/bin/env Rscript

###############################################################################
# script name: asv_to_community.R
# developed by: Savvas Paragkamian
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is provide a streamline of files with all the information
# required for further analyses. That is to use ASVs and the taxonomy to filter
# and normalise the reads of each sample biodiversity and to combine sample 
# metadata.
#
###############################################################################
# OUTPUT:
#
# the main output of this script is 
# Results/crete_biodiversity_asv.tsv
# which contains the sample - asv occurrences with abundance along with 
# taxonomic information.
# 
###############################################################################
# RUNNING TIME: 9 minutes
###############################################################################
# usage:./Scripts/isd_crete_biodiversity.R
###############################################################################
source("Scripts/functions.R")
library(magrittr)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)

#### set directory
create_dir("Results")
create_dir("Figures")

################################## Load data ##################################
# Load data
print("loading data")

abundance_asv <- readRDS("Results/dada2_output/taxonomy/seqtab_nochim.RDS") # use this one if data are RDS
asv_fasta <- read_delim("Results/dada2_output/taxonomy/asv_fasta_ids.tsv", delim="\t")
taxa_asv <- readRDS("Results/dada2_output/taxonomy/dada2_taxonomy.RDS")
species_asv <- readRDS("Results/dada2_output/taxonomy/dada2_taxa_species.RDS")

########################## transform the matrices to long tables ##########################

print("transforming data")
## this is for the species exact matching results from DADA2 addSpecies function
species <- data.frame(genus = species_asv[,6], species = species_asv[,8]) |> 
    rownames_to_column("taxonomic_unit") |>
    mutate(speciesname = paste(genus, species,sep=" ")) |>
    group_by(speciesname) |>
    summarise(n=n())

## the abundance matrix has all the biodiverstity information
## sampleID, ASV and abundance. NOTE that contains many zeros!!!
## UNCOMMENT if data are from RDS!
abundance_asv_long <- abundance_asv |> 
    data.frame() |>
    rownames_to_column("file") |>
    pivot_longer(!file, names_to = "sequence", values_to = "abundance") |>
    mutate(ENA_RUN=gsub("_1_filt.fastq.gz", "", file, perl=T))

## the taxonomy matrix has all the taxonomic information of each ASV
taxa_asv_l <- taxa_asv |> 
    data.frame() |>
    rownames_to_column("sequence") |>
    mutate(Species = ifelse(is.na(Species),
                                 NA,
                                 paste(Genus, Species,sep=" "))) |>
    pivot_longer(!sequence,
                 names_to = "higherClassification",
                 values_to = "scientificName") |>
    na.omit("scientificName") |>
    ungroup() |>
    left_join(asv_fasta, by=c("sequence"="asv")) |>
    rename("taxonomic_unit"="asv_id")


taxa_asv_s <- taxa_asv_l |>
    group_by(taxonomic_unit) |>
    summarise(higherClassification = paste0(higherClassification,
                                            collapse = "; "), 
              taxonomy = paste0(scientificName,
                                       collapse = "; ")) |>
    ungroup()

# keep the lowest taxonomic inforfation per ASV
taxa_asv_s$scientificName <- sapply(taxa_asv_s$taxonomy,
                                           keep_last,
                                           simplify = T, USE.NAMES=F)

taxa_asv_s$classification <- sapply(taxa_asv_s$higherClassification,
                                           keep_last,
                                           simplify = T, USE.NAMES=F)

taxa_asv_all <- taxa_asv_l |> 
    pivot_wider(names_from=higherClassification, 
                values_from=scientificName) |> 
    left_join(taxa_asv_s, by=c("taxonomic_unit"="taxonomic_unit"))

####################### merge abundance and taxonomy ##########################

crete_biodiversity_all <- abundance_asv_long |>
    left_join(taxa_asv_all,by=c("sequence"="sequence")) 

## total ASVs and taxonomy
asv_no_taxonomy <- crete_biodiversity_all |> 
    distinct(taxonomic_unit,classification) |>
    filter(is.na(classification))

print(paste0("ASVs without taxonomy: ", nrow(asv_no_taxonomy)))

########################## Master file ####################################
crete_biodiversity_asv <- crete_biodiversity_all |>
    filter(abundance>0) |>
    filter(!is.na(taxonomy)) 

write_delim(crete_biodiversity_asv, "Results/crete_biodiversity_asv.tsv", delim="\t")
print("master file created")
