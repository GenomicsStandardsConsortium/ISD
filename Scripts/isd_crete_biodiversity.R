#!/usr/bin/env Rscript

###############################################################################
# script name: isd_crete_biodiversity.R
# developed by: Savvas Paragkamian, Johanna Holms
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of the script is to perform biodiversity analysis of the community matrix
# and provide sample metadata.
# It is the intermediate step before the numerical ecology analysis.
#
###############################################################################
# OUTPUT:
# 
# Other 4 files are produced
# crete_biodiversity_matrix.RDS, a matrix of abundances
# tax_tab.RDS, taxonomy table with the remaining asvs
# sample_stats.tsv
# phyla_samples_summary.tsv
###############################################################################
# RUNNING TIME: 9 minutes
###############################################################################
# usage:
# ./Scripts/isd_crete_biodiversity.R Results/crete_biodiversity_otu.tsv "OTU"
###############################################################################
source("Scripts/functions.R")
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(SRS) # normalisation
library(vegan)

################################## User input #################################
# test if there is at least one argument: if not, return an error
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least two arguments must be supplied (input file) and a prefix.n", call.=FALSE)
} else if (length(args)==1) {
  stop("At least two arguments must be supplied (input file) and a prefix.n", call.=FALSE)
} 

abundance_table_path <- args[1]
prefix <- args[2]
################################## Load data ##################################

# Load data
print("loading data")

crete_biodiversity_all <- read_delim(abundance_table_path, delim="\t")


#### set directory
create_dir("Results")
create_dir("Figures")
################################### Metadata ##################################

metadata_long <- read_delim("Data/Metadata/ena_isd_2016_attributes.tsv", delim="\t") |>
    mutate(VALUE=gsub("\\r(?!\\n)","", VALUE, perl=T))
# metadata to wide format

metadata_wide <- metadata_long |> 
    dplyr::select(-c(UNITS)) |>
    filter(!(TAG %in% c("ENA-FIRST-PUBLIC", "ENA-LAST-UPDATE"))) |>
    mutate(TAG=gsub(" ","_", TAG, perl=T)) |>
    mutate(TAG=gsub("-","_", TAG, perl=T)) |>
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

metadata <- metadata_wide |>
    dplyr::select(ENA_RUN, source_material_identifiers, total_nitrogen, water_content,total_organic_carbon,sample_volume_or_weight_for_DNA_extraction,DNA_concentration,latitude, longitude,elevation, amount_or_size_of_sample_collected, vegetation_zone) |>
    arrange(ENA_RUN) |>
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


######################## Decontamination #####################################
## very important step, we need the samples controls.

# This is the MASTER dataset of Crete biodiversity and taxonomy

### move from here the rest

############################ create a abundance matrix ########################
crete_biodiversity_m <- crete_biodiversity_all |>
    filter(!is.na(taxonomy)) |>
    dplyr::select(ENA_RUN, taxonomic_unit, abundance) |>
    pivot_wider(names_from=ENA_RUN, values_from=abundance, values_fill = 0) |>
    as.matrix()

crete_biodiversity_matrix <- crete_biodiversity_m[,-1]
crete_biodiversity_matrix <- apply(crete_biodiversity_matrix, 2, as.numeric)
rownames(crete_biodiversity_matrix) <- crete_biodiversity_m[,1]
saveRDS(crete_biodiversity_matrix, paste0("Results/",prefix,"_crete_matrix.RDS", sep=""))
print("master file created")

########################## Normalisation Compositionality #####################
# SRS compositional preparation
print("compositionanlity with SRS method")

Cmin <- min(colSums(crete_biodiversity_matrix))
Cmax <- max(colSums(crete_biodiversity_matrix))
summary(colSums(crete_biodiversity_matrix))
table(colSums(crete_biodiversity_matrix) > 10000)

## SRS curve
png(file=paste0("Figures/",prefix,"_isd_crete_srs_curve.jpeg", sep=""),
    width = 50,
    height = 30,
    res=300,
    units = "cm",
    bg="white")

SRScurve(crete_biodiversity_matrix,
         metric = "richness",
         step = 500,
         xlim=c(0,Cmax),
         xlab = "sample size",
         ylab = "richness",
         label = F,
         col = c("#5A4A6F", "#E47250",  "#EBB261"))

dev.off()

crete_biodiversity_matrix_df <- as.data.frame(crete_biodiversity_matrix)
rownames(crete_biodiversity_matrix_df) <- rownames(crete_biodiversity_matrix)

biodiversity_srs <- SRS(crete_biodiversity_matrix_df,
                        10000,
                        set_seed = TRUE,
                        seed = 1)

rownames(biodiversity_srs) <- rownames(crete_biodiversity_matrix_df)

# how many samples don't have taxonomic units
print("how many samples don't have ASVs/OTUs")
length(which(colSums(biodiversity_srs)==0))
# how many ASVs have 0 abundance after the SRS
print("how many samples have 0 abundance after the SRS")
length(which(rowSums(biodiversity_srs)==0))
dim(biodiversity_srs)
## filter empty
biodiversity_srs <- biodiversity_srs[-which(rowSums(biodiversity_srs)==0) ,]

#### Save ######
#saveRDS(biodiversity_srs, "results/biodiversity_srs.RDS")
#biodiversity_srs <- readRDS("results/biodiversity_srs.RDS")

### remove the asvs that don't have a taxonomy
############################# Merge SRS ################################
biodiversity_srs_l <- dist_long(biodiversity_srs,"srs_abundance")

crete_biodiversity <- crete_biodiversity_all |>
    filter(!is.na(classification)) |> 
    left_join(biodiversity_srs_l, by=c("taxonomic_unit"="rowname", "ENA_RUN"="colname"))

write_delim(crete_biodiversity,paste0("Results/",prefix,"_crete_biodiversity.tsv", sep=""),delim="\t")

#crete_biodiversity <- read_delim("results/crete_biodiversity_asv.tsv",delim="\t")

print("how many occurrences don't have srs abundace")
nrow(crete_biodiversity[is.na(crete_biodiversity$srs_abundance),])

print("how many occurrences with srs abundace and Phylum")
nrow(crete_biodiversity[!is.na(crete_biodiversity$srs_abundance) & !is.na(crete_biodiversity$Phylum),])

print("how many occurrences with srs abundace and Family")
nrow(crete_biodiversity[!is.na(crete_biodiversity$srs_abundance) & !is.na(crete_biodiversity$Family),])

print("how many occurrences with srs abundace and Genus")
nrow(crete_biodiversity[!is.na(crete_biodiversity$srs_abundance) & !is.na(crete_biodiversity$Genus),])

print("how many occurrences with srs abundace and Species")
nrow(crete_biodiversity[!is.na(crete_biodiversity$srs_abundance) & !is.na(crete_biodiversity$Species),])

################################ Community Matrix ###########################################
## here is the last filtering

community_matrix_l <- crete_biodiversity |>
    filter(srs_abundance>0,
           !is.na(srs_abundance),
           !is.na(Phylum)) |>
    group_by(ENA_RUN,Kingdom,Phylum,Class,Order,Family,Genus,Species,scientificName,classification,taxonomy) |>
    summarise(n_taxonomic_units=n(),
              reads_srs_mean=mean(srs_abundance),
              reads_srs_sum=sum(srs_abundance), .groups="keep") |>
    group_by(ENA_RUN) |>
    mutate(relative_srs=reads_srs_sum/sum(reads_srs_sum)) |>
    ungroup()

write_delim(community_matrix_l,paste0("Results/",prefix,"_community_matrix_l.tsv", sep=""),delim="\t")
#community_matrix_l <- read_delim("results/community_matrix_l.tsv",delim="\t") 

community_matrix <- community_matrix_l |>
    dplyr::select(ENA_RUN,reads_srs_sum,scientificName) |>
    pivot_wider(names_from=scientificName,
                values_from=reads_srs_sum,
                values_fill=0) |>
    as.data.frame()

write_delim(community_matrix,paste0("Results/",prefix,"_community_matrix.tsv", sep=""), delim="\t")

rownames(community_matrix) <- community_matrix[,1]
community_matrix <- community_matrix[,-1]
saveRDS(community_matrix, paste0("Results/",prefix,"_community_matrix.RDS", sep=""))

################################ save Faprotax format ############################
#### faprotax community matrix

faprotax_community_matrix <- community_matrix_l |>
    pivot_wider(id_cols=c(scientificName,taxonomy),
                names_from=ENA_RUN,
                values_from=reads_srs_sum, values_fill=0) |>
    relocate(taxonomy, .after = last_col())

write_delim(faprotax_community_matrix,paste0("Results/",prefix,"_faprotax_community_matrix.tsv"),delim="\t")

######################## clean environment ######################## 

#variables_keep <- c(biodiversity_srs_l,
#                    biodiversity_srs,
#                    crete_biodiversity,
#                    metadata)

#rm(list=setdiff(ls(), variables_keep))

########################## Samples diversity #########################
print("sample summary")
################################# Indices ###################################
# Explore alpha-metrics

shannon <- data.frame(shannon = as.matrix(diversity(community_matrix, index = "shannon")))
observed <- data.frame(t(estimateR(community_matrix)))

biodiversity_index <- cbind(shannon,observed)

biodiversity_index$ENA_RUN <- rownames(biodiversity_index)
rownames(biodiversity_index) <- NULL
biodiversity_index <- as.data.frame(biodiversity_index)
## taxonomic, asv and read diversity per sample

sample_taxonomy_stats <- community_matrix_l |> 
    group_by(ENA_RUN,classification) |>
    summarise(taxa=n(),
              n_taxonomic_units=sum(n_taxonomic_units),
              reads_srs=sum(reads_srs_sum),
              .groups="keep") |>
#    pivot_wider(names_from=classification,values_from=n_taxa) |>
    dplyr::ungroup()

write_delim(sample_taxonomy_stats,
            paste0("Results/",prefix,"_sample_taxonomy_stats.tsv"),
            delim="\t")

sample_stats_total <- sample_taxonomy_stats |>
    group_by(ENA_RUN) |>
    summarise(taxa=sum(taxa),
              n_taxonomic_units=sum(n_taxonomic_units),  # total ASV per sample
              reads_srs=sum(reads_srs))

### highest species biodiversity sample

crete_biodiversity_s <- community_matrix_l |>
    distinct(ENA_RUN, Species) |> 
    group_by(ENA_RUN) |>
    summarise(Species=n())

print("sample with the highest microbial species diversity")
crete_biodiversity_s[which(crete_biodiversity_s$Species==max(crete_biodiversity_s$Species)),]

## keep only the sample metadata after filtering
## filter also metadata and taxonomy
metadata_all <-  metadata |>
    left_join(biodiversity_index, by=c("ENA_RUN"="ENA_RUN")) |>
    left_join(sample_stats_total, by=c("ENA_RUN"="ENA_RUN"))

write_delim(metadata_all,
            paste0("Results/",prefix,"_sample_metadata.tsv"),
            delim="\t")

########################## taxonomic units summary ###############################
print("Taxonomic summary")
# Create taxonomy table of the remaining taxonomic units

taxonomic_units_stats <- crete_biodiversity |> 
    group_by(taxonomic_unit, classification, scientificName) |> 
    summarise(n_samples=n(),
              reads=sum(abundance),
              reads_srs=sum(srs_abundance, na.rm=T),
              reads_srs_mean=mean(srs_abundance, na.rm=T),
              reads_srs_sd=sd(srs_abundance, na.rm=T),
              .groups="keep") |>
    dplyr::ungroup()

write_delim(taxonomic_units_stats,
            paste0("Results/",prefix,"_taxonomic_units_metadata.tsv"),
            delim="\t")

## singletons
singletons <- taxonomic_units_stats |> 
    filter(reads==1) |>
    nrow()

print(paste0("there are ",singletons," singletons"))

## taxonomic units and samples distribution
taxonomic_units_sample_dist <- taxonomic_units_stats |>
    group_by(n_samples) |>
    summarise(n_taxonomic_units=n(),
              reads=sum(reads),
              reads_srs=sum(reads_srs, na.rm=T))


################################# Taxonomy #####################################
tax_tab1 <- community_matrix_l |>
    distinct(scientificName, Kingdom, Phylum, Class, Order, Family, Genus, Species) |>
    as.matrix()

tax_tab <- tax_tab1[,-1]
rownames(tax_tab) <- tax_tab1[,1]
saveRDS(tax_tab, paste0("Results/",prefix,"_tax_tab.RDS"))

print("taxonomic summary")
## how the information of communities of Cretan soils is distributed across 
## the taxonomic levels

taxonomy_taxa <- community_matrix_l |>
    distinct(classification,scientificName) |>
    group_by(classification) |> summarise(n_taxa=n())

taxonomy_levels_occurrences <- community_matrix_l |>
    group_by(classification) |>
    summarise(n_occurrences=n(),
              reads_srs=sum(reads_srs_sum, na.rm=T),
              n_taxonomic_units=sum(n_taxonomic_units))

write_delim(taxonomy_levels_occurrences,
            paste0("Results/",prefix,"_taxonomy_levels_occurrences.tsv"),
            delim="\t")
## Phyla distribution, average relative abundance and ubiquity

phyla_samples_summary <- community_matrix_l |>
    group_by(ENA_RUN,Phylum) |>
    summarise(n_taxonomic_units=sum(n_taxonomic_units),
              abundance_mean=mean(reads_srs_sum),
              abundance_sum=sum(reads_srs_sum), .groups="keep") |>
    group_by(ENA_RUN) |>
    mutate(relative_abundance=abundance_sum/sum(abundance_sum))
#    na.omit(Phylum)

write_delim(phyla_samples_summary,
            paste0("Results/",prefix,"_phyla_samples_summary.tsv"),
            delim="\t")

## phyla stats
total_samples <- length(unique(community_matrix_l$ENA_RUN))
phyla_stats <- phyla_samples_summary |> 
    group_by(Phylum) |>
    summarise(n_samples=n(),
              total_taxonomic_units=sum(n_taxonomic_units),
              total_abundance=sum(abundance_sum),
              proportion_sample=n_samples/total_samples,
              average_relative=mean(relative_abundance)) |>
    arrange(desc(average_relative))

write_delim(phyla_stats,
            paste0("Results/",prefix,"_phyla_stats.tsv"),
            delim="\t")

## Genera stats

genera_phyla_samples <- community_matrix_l |>
    filter(!is.na(Genus)) |>
    group_by(Phylum,Genus,ENA_RUN) |>
    summarise(n_taxonomic_units=sum(n_taxonomic_units),
              reads_srs_mean=mean(reads_srs_sum),
              reads_srs_sum=sum(reads_srs_sum), .groups="keep") |>
    group_by(ENA_RUN) |>
    mutate(relative_srs=reads_srs_sum/sum(reads_srs_sum)) |>
    ungroup()

write_delim(genera_phyla_samples,
            paste0("Results/",prefix,"_genera_phyla_samples.tsv"),
            delim="\t")

genera_phyla_stats <- genera_phyla_samples |>
    group_by(Phylum,Genus) |>
    summarise(n_samples=n(),
              relative_abundance_mean=mean(relative_srs), 
              relative_abundance_sd=sd(relative_srs), 
              total_taxonomic_units=sum(n_taxonomic_units),
              reads_srs_mean=mean(reads_srs_sum),
              reads_srs_sd=sd(reads_srs_sum),
              total_reads_srs=sum(reads_srs_sum),
              proportion_sample=n_samples/total_samples, .groups="keep") |>
    group_by(Phylum) |>
    mutate(n_genera=n()) |>
    arrange(desc(n_samples), desc(relative_abundance_mean))

write_delim(genera_phyla_stats,
            paste0("Results/",prefix,"_genera_phyla_stats.tsv"),
            delim="\t")

################################ save network format ############################

network_genera_community_matrix <- community_matrix_l |>
    pivot_wider(id_cols=c(ENA_RUN),
                names_from=scientificName,
                values_from=reads_srs_sum,
                values_fill=0) |>
    column_to_rownames("ENA_RUN")

write_delim(network_genera_community_matrix,
            paste0("Results/",prefix,"_gnetwork_genera_community_matrix.tsv"),
            delim="\t")

##################################################################################

print("biodiversity script finished")
