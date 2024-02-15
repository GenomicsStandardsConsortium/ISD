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
### abundance table
abundance_table <- read_delim("Results/pema215_vsearch/my_taxon_assign/finalTable.tsv", delim="\t") |>
    rename(taxonomy = Classification) |>
    mutate(taxonomy=gsub("Main genome;","",taxonomy)) |>
    mutate(taxonomy=gsub(";","; ",taxonomy))

### sequences
fasta_q <- readLines("Results/pema215_vsearch/all.otus.fasta")
ids <- grepl(">", fasta_q)

fasta <- data.frame(
  OTU = sub(">", "", fasta_q[ids]),
  sequence = tapply(fasta_q[!ids], cumsum(ids)[!ids], function(x) {
    paste(x, collapse = "")
  }))

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

##################################### taxonomy ###################################
colnames <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_tab_tmp <- abundance_table |>
    dplyr::select(OTU,taxonomy) |>
    separate(taxonomy,
                     colnames,
                     sep = "; ",
                     remove = FALSE,
                     convert = FALSE,
                     extra = "warn",
                     fill = "warn") 

tax_tab_tmp$keep_last <- sapply(tax_tab_tmp$taxonomy,
                                           keep_last,
                                           simplify = T, USE.NAMES=F)
tax_tab <- tax_tab_tmp |> 
    pivot_longer(-c(OTU,taxonomy,keep_last),
                 names_to="classification",
                 values_to="scientificName") |>
    filter(scientificName==keep_last) |>
    dplyr::select(-keep_last) |>
    separate(taxonomy,
                 colnames,
                 sep = "; ",
                 remove = FALSE,
                 convert = FALSE,
                 extra = "warn",
                 fill = "warn")

tax_tab_tmp_no <- tax_tab_tmp[which(!c(tax_tab_tmp$OTU %in% tax_tab$OTU)),]
print("there are some OTUs that are excluded, Chroroplasts")
nrow(tax_tab_tmp_no)
########################## transform the matrices to long tables ##########################

print("transforming data")

####################### merge abundance and taxonomy ##########################
## the abundance matrix has all the biodiverstity information
abundance_table_long <- abundance_table |>
    pivot_longer(!c(OTU,taxonomy),
                 names_to = "ENA_RUN",
                 values_to = "abundance") |>
    dplyr::select(-taxonomy) |> 
    left_join(tax_tab, by=c("OTU"="OTU"))


crete_biodiversity_all <- abundance_table_long |>
    left_join(fasta, by=c("OTU"="OTU")) |>
    mutate(taxonomy=if_else(taxonomy=="No hits", NA, taxonomy))

########################### Master file ######################################
## total ASVs and taxonomy
print("there are some OTUs that are excluded, Chroroplasts and No hits")
no_taxonomy <- crete_biodiversity_all %>% 
    distinct(OTU,taxonomy) %>%
    filter(is.na(taxonomy))

print(paste0("OTUs without taxonomy: ", nrow(no_taxonomy)))

crete_biodiversity_otu <- crete_biodiversity_all |>
    filter(abundance>0) |>
    filter(!is.na(taxonomy))

write_delim(crete_biodiversity_otu, "Results/crete_biodiversity_otu.tsv", delim="\t")

### move from here the rest

## create a abundance matrix
crete_biodiversity_m <- crete_biodiversity_all %>%
    filter(!is.na(taxonomy)) %>%
    dplyr::select(ENA_RUN, OTU, abundance) %>%
    pivot_wider(names_from=ENA_RUN, values_from=abundance, values_fill = 0) %>%
    as.matrix()

crete_biodiversity_matrix <- crete_biodiversity_m[,-1]
crete_biodiversity_matrix <- apply(crete_biodiversity_matrix, 2, as.numeric)
rownames(crete_biodiversity_matrix) <- crete_biodiversity_m[,1]
saveRDS(crete_biodiversity_matrix, "Results/crete_otu_matrix.RDS")
print("master file created")
