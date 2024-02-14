#!/usr/bin/env Rscript

###############################################################################
# script name: otu_to_community.R
# developed by: Savvas Paragkamian
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is provide a streamline of files with all the information
# required for further analyses. That is to use OTUs and the taxonomy to filter
# and normalise the reads of each sample biodiversity and to combine sample 
# metadata.
#
###############################################################################
# OUTPUT:
#
# the main output of this script is 
# 0. community_matrix_l.tsv
# 1. crete_biodiversity_asv.tsv
# which contains the sample - asv occurrences with abundance along with 
# taxonomic information.
# 2. sample_metadata.tsv
# 3. asv_metadata.tsv
# 
# Other 4 files are produced
# crete_biodiversity_matrix.RDS, a matrix of abundances
# tax_tab.RDS, taxonomy table with the remaining asvs
# sample_stats.tsv
# phyla_samples_summary.tsv
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
abundance_table <- read_delim("Results/pema215_vsearch/my_taxon_assign/finalTable.tsv", delim="\t") |>
    rename(taxonomy = Classification) |>
    mutate(taxonomy=gsub("Main genome;","",taxonomy))

fasta_q <- readLines("Results/pema215_vsearch/all.otus.fasta")
ids <- grepl(">", fasta_q)

fasta <- data.frame(
  OTU = sub(">", "", fasta_q[ids]),
  sequence = tapply(fasta_q[!ids], cumsum(ids)[!ids], function(x) {
    paste(x, collapse = "")
  }))

#silva_132 <- read_delim("Data/Metadata/tax_slv_ssu_132.txt", delim="\t", col_names=F)

#colnames(silva_132) <- c("taxonomy", "n_taxa", "classification", "5")

#silva_132 <- silva_132 |>
#    dplyr::select(taxonomy,classification) |>
#    mutate(taxonomy=gsub(";$","",taxonomy)) 

#silva_132$scientificName <- sapply(silva_132$taxonomy,
#                                           keep_last_2,
#                                           simplify = T, USE.NAMES=F)

colnames <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_tab_tmp <- abundance_table |>
    dplyr::select(OTU,taxonomy)

tax_tab <- separate(tax_tab_tmp,
                     taxonomy,
                     colnames,
                     sep = ";",
                     remove = TRUE,
                     convert = FALSE,
                     extra = "warn",
                     fill = "warn")
################################### Metadata ####################################

metadata_long <- read_delim("Data/Metadata/ena_isd_2016_attributes.tsv", delim="\t") %>%
    mutate(VALUE=gsub("\\r(?!\\n)","", VALUE, perl=T))
# metadata to wide format

metadata_wide <- metadata_long %>% 
    dplyr::select(-c(UNITS)) %>%
    filter(!(TAG %in% c("ENA-FIRST-PUBLIC", "ENA-LAST-UPDATE"))) %>%
    mutate(TAG=gsub(" ","_", TAG, perl=T)) %>%
    mutate(TAG=gsub("-","_", TAG, perl=T)) %>%
    pivot_wider(names_from=TAG, 
                values_from=VALUE)

metadata_wide$total_nitrogen <- as.numeric(metadata_wide$total_nitrogen)
metadata_wide$water_content <- as.numeric(metadata_wide$water_content)
metadata_wide$total_organic_carbon <- as.numeric(metadata_wide$total_organic_carbon)
metadata_wide$sample_volume_or_weight_for_DNA_extraction <- as.numeric(metadata_wide$sample_volume_or_weight_for_DNA_extraction)
metadata_wide$DNA_concentration <- as.numeric(metadata_wide$dna_concentration)
metadata_wide$latitude <- as.numeric(metadata_wide$`geographic_location_(latitude)`)
metadata_wide$longitude <- as.numeric(metadata_wide$`geographic_location_(longitude)`)
metadata_wide$elevation <- as.numeric(metadata_wide$`geographic_location_(elevation)`)
metadata_wide$amount_or_size_of_sample_collected <- as.numeric(metadata_wide$amount_or_size_of_sample_collected)

metadata <- metadata_wide %>%
    dplyr::select(ENA_RUN, source_material_identifiers, total_nitrogen, water_content,total_organic_carbon,sample_volume_or_weight_for_DNA_extraction,DNA_concentration,latitude, longitude,elevation, amount_or_size_of_sample_collected, vegetation_zone) %>%
    arrange(ENA_RUN) %>%
    mutate(route = sub("isd_(.*.)_site.*" ,"\\1" , source_material_identifiers))

## location pairs of each site
metadata$sites <- gsub("_loc_*.","", metadata$source_material_identifiers)
metadata$location <- gsub(".*(loc_.*)","\\1", metadata$source_material_identifiers)

## C:N ratio
## Six samples have 0 total_nitrogen. So to avoid the inf 
## use the if else statement.
metadata$carbon_nitrogen_ratio <- ifelse(metadata$total_nitrogen==0,
                                         metadata$total_organic_carbon,
                                         metadata$total_organic_carbon/metadata$total_nitrogen)

########################## transform the matrices to long tables ##########################

print("transforming data")

## the abundance matrix has all the biodiverstity information
## sampleID, ASV and abundance. NOTE that contains many zeros!!!
## UNCOMMENT if data are from RDS!
abundance_long <- abundance_table |>
    pivot_longer(!c(OTU,taxonomy),
                 names_to = "ENA_RUN",
                 values_to = "abundance") |>
    left_join(tax_tab, by=c("OTU"="OTU")) 


abundance_long_all <- abundance_long |>
    left_join(fasta, by=c("OTU"="OTU"))


# keep the lowest taxonomic inforfation per ASV
taxa_asv_s$scientificName <- sapply(taxa_asv_s$taxonomy,
                                           keep_last_2,
                                           simplify = T, USE.NAMES=F)

taxa_asv_s$classification <- sapply(taxa_asv_s$higherClassification,
                                           keep_last,
                                           simplify = T, USE.NAMES=F)

taxa_asv_all <- taxa_asv_l %>% 
    pivot_wider(names_from=higherClassification, 
                values_from=scientificName) %>% 
    left_join(taxa_asv_s, by=c("asv"="asv"))

####################### merge abundance and taxonomy ##########################

crete_biodiversity_all <- abundance_asv_long %>%
    left_join(taxa_asv_all,by=c("asv"="asv")) %>%
    filter(abundance>0)

########################### Master file ######################################
## total ASVs and taxonomy
asv_no_taxonomy <- crete_biodiversity_all %>% 
    distinct(asv,classification) %>%
    filter(is.na(classification))

print(paste0("ASVs without taxonomy: ", nrow(asv_no_taxonomy)))
## create a abundance matrix
crete_biodiversity_m <- crete_biodiversity_all %>%
    filter(!is.na(classification)) %>%
    dplyr::select(ENA_RUN, asv_id, abundance) %>%
    pivot_wider(names_from=ENA_RUN, values_from=abundance, values_fill = 0) %>%
    as.matrix()

crete_biodiversity_matrix <- crete_biodiversity_m[,-1]
crete_biodiversity_matrix <- apply(crete_biodiversity_matrix, 2, as.numeric)
rownames(crete_biodiversity_matrix) <- crete_biodiversity_m[,1]
saveRDS(crete_biodiversity_matrix, "results/crete_biodiversity_matrix.RDS")
print("master file created")
