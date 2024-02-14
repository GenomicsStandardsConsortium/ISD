#!/usr/bin/env Rscript

###############################################################################
# script name: isd_crete_biodiversity.R
# developed by: Savvas Paragkamian, Johanna Holms
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of the script is to perform biodiversity and numerical ecology analysis
# of the community matrix
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
# usage:./Scripts/isd_crete_biodiversity.R
###############################################################################
source("Scripts/functions.R")
library(magrittr)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(vegan)
library(SRS) # normalisation

#### set directory
create_dir("Results")
create_dir("Figures")

################################## Load data ##################################

# Load data
print("loading data")

######################## Decontamination #####################################
## very important step, we need the samples controls.

# This is the MASTER dataset of Crete biodiversity and taxonomy


########################## Normalisation Compositionality #####################
# SRS compositional preparation
print("compositionanlity with SRS method")

Cmin <- min(colSums(crete_biodiversity_matrix))
Cmax <- max(colSums(crete_biodiversity_matrix))
summary(colSums(crete_biodiversity_matrix))
table(colSums(crete_biodiversity_matrix) > 10000)

## SRS curve
png(file="figures/isd_crete_srs_curve.jpeg",
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

# how many samples don't have ASVs
print("how many samples don't have ASVs")
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

crete_biodiversity <- crete_biodiversity_all %>%
    filter(!is.na(classification)) %>% 
    left_join(biodiversity_srs_l, by=c("asv_id"="rowname", "ENA_RUN"="colname"))

write_delim(crete_biodiversity,"results/crete_biodiversity_asv.tsv",delim="\t")

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

community_matrix_l <- crete_biodiversity %>%
    filter(srs_abundance>0,
           !is.na(srs_abundance),
           !is.na(Phylum)) %>%
    group_by(ENA_RUN,Kingdom,Phylum,Class,Order,Family,Genus,Species,scientificName,classification,taxonomy) %>%
    summarise(asvs=n(),
              reads_srs_mean=mean(srs_abundance),
              reads_srs_sum=sum(srs_abundance), .groups="keep") %>%
    group_by(ENA_RUN) %>%
    mutate(relative_srs=reads_srs_sum/sum(reads_srs_sum)) %>%
    ungroup()

write_delim(community_matrix_l,"results/community_matrix_l.tsv",delim="\t")
#community_matrix_l <- read_delim("results/community_matrix_l.tsv",delim="\t") 

community_matrix <- community_matrix_l %>%
    dplyr::select(ENA_RUN,reads_srs_sum,scientificName) %>%
    pivot_wider(names_from=scientificName,
                values_from=reads_srs_sum,
                values_fill=0) %>%
    as.data.frame()

write_delim(community_matrix,"results/community_matrix.tsv", delim="\t")

rownames(community_matrix) <- community_matrix[,1]
community_matrix <- community_matrix[,-1]
saveRDS(community_matrix, "results/community_matrix.RDS")

################################ save Faprotax format ############################
#### faprotax community matrix

faprotax_community_matrix <- community_matrix_l %>%
    pivot_wider(id_cols=c(scientificName,taxonomy),
                names_from=ENA_RUN,
                values_from=reads_srs_sum, values_fill=0) %>%
    relocate(taxonomy, .after = last_col())

write_delim(faprotax_community_matrix,"results/faprotax_community_matrix.tsv",delim="\t")


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

sample_taxonomy_stats <- community_matrix_l %>% 
    group_by(ENA_RUN,classification) %>%
    summarise(taxa=n(),
              asvs=sum(asvs),
              reads_srs=sum(reads_srs_sum),
              .groups="keep") %>%
#    pivot_wider(names_from=classification,values_from=n_taxa) %>%
    dplyr::ungroup()

write_delim(sample_taxonomy_stats,
            "results/sample_taxonomy_stats.tsv",
            delim="\t")

sample_stats_total <- sample_taxonomy_stats %>%
    group_by(ENA_RUN) %>%
    summarise(taxa=sum(taxa),
              asvs=sum(asvs),  # total ASV per sample
              reads_srs=sum(reads_srs))

### highest species biodiversity sample

crete_biodiversity_s <- community_matrix_l %>%
    distinct(ENA_RUN, Species) %>% 
    group_by(ENA_RUN) %>%
    summarise(Species=n())

print("sample with the highest microbial species diversity")
crete_biodiversity_s[which(crete_biodiversity_s$Species==max(crete_biodiversity_s$Species)),]

## keep only the sample metadata after filtering
## filter also metadata and taxonomy
metadata_all <-  metadata %>%
    left_join(biodiversity_index, by=c("ENA_RUN"="ENA_RUN")) %>%
    left_join(sample_stats_total, by=c("ENA_RUN"="ENA_RUN"))

write_delim(metadata_all,
            "results/sample_metadata.tsv",
            delim="\t")

########################## ASV summary ###############################
print("ASV summary")
# Create taxonomy table of the remaining asvs

asv_stats <- crete_biodiversity %>% 
    group_by(asv_id, classification, scientificName) %>% 
    summarise(n_samples=n(),
              reads=sum(abundance),
              reads_srs=sum(srs_abundance, na.rm=T),
              reads_srs_mean=mean(srs_abundance, na.rm=T),
              reads_srs_sd=sd(srs_abundance, na.rm=T),
              .groups="keep") %>%
    dplyr::ungroup()

write_delim(asv_stats,"results/asv_metadata.tsv",delim="\t")

## singletons
singletons <- asv_stats %>% 
    filter(reads==1) %>%
    nrow()

print(paste0("there are ",singletons," singletons asvs"))

## asv and samples distribution
asv_sample_dist <- asv_stats %>%
    group_by(n_samples) %>%
    summarise(n_asv=n(),
              reads=sum(reads),
              reads_srs=sum(reads_srs, na.rm=T))


################################# Taxonomy #####################################
tax_tab1 <- community_matrix_l %>%
    distinct(scientificName, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
    as.matrix()

tax_tab <- tax_tab1[,-1]
rownames(tax_tab) <- tax_tab1[,1]
saveRDS(tax_tab, "results/tax_tab.RDS")

print("taxonomic summary")
## how the information of communities of Cretan soils is distributed across 
## the taxonomic levels

taxonomy_taxa <- community_matrix_l %>%
    distinct(classification,scientificName) %>%
    group_by(classification) %>% summarise(n_taxa=n())

taxonomy_levels_occurrences <- community_matrix_l %>%
    group_by(classification) %>%
    summarise(n_occurrences=n(),
              reads_srs=sum(reads_srs_sum, na.rm=T),
              asvs=sum(asvs))

write_delim(taxonomy_levels_occurrences,"results/taxonomy_levels_occurrences.tsv",delim="\t")
## Phyla distribution, average relative abundance and ubiquity

phyla_samples_summary <- community_matrix_l %>%
    group_by(ENA_RUN,Phylum) %>%
    summarise(asvs=sum(asvs),
              abundance_mean=mean(reads_srs_sum),
              abundance_sum=sum(reads_srs_sum), .groups="keep") %>%
    group_by(ENA_RUN) %>%
    mutate(relative_abundance=abundance_sum/sum(abundance_sum))
#    na.omit(Phylum)

write_delim(phyla_samples_summary,"results/phyla_samples_summary.tsv",delim="\t")

## phyla stats
total_samples <- length(unique(community_matrix_l$ENA_RUN))
phyla_stats <- phyla_samples_summary %>% 
    group_by(Phylum) %>%
    summarise(n_samples=n(),
              total_asvs=sum(asvs),
              total_abundance=sum(abundance_sum),
              proportion_sample=n_samples/total_samples,
              average_relative=mean(relative_abundance)) %>%
    arrange(desc(average_relative))

write_delim(phyla_stats,"results/phyla_stats.tsv",delim="\t")

## Genera stats

genera_phyla_samples <- community_matrix_l %>%
    filter(!is.na(Genus)) %>%
    group_by(Phylum,Genus,ENA_RUN) %>%
    summarise(asvs=sum(asvs),
              reads_srs_mean=mean(reads_srs_sum),
              reads_srs_sum=sum(reads_srs_sum), .groups="keep") %>%
    group_by(ENA_RUN) %>%
    mutate(relative_srs=reads_srs_sum/sum(reads_srs_sum)) %>%
    ungroup()

write_delim(genera_phyla_samples,"results/genera_phyla_samples.tsv",delim="\t")

genera_phyla_stats <- genera_phyla_samples %>%
    group_by(Phylum,Genus) %>%
    summarise(n_samples=n(),
              relative_abundance_mean=mean(relative_srs), 
              relative_abundance_sd=sd(relative_srs), 
              total_asvs=sum(asvs),
              reads_srs_mean=mean(reads_srs_sum),
              reads_srs_sd=sd(reads_srs_sum),
              total_reads_srs=sum(reads_srs_sum),
              proportion_sample=n_samples/total_samples, .groups="keep") %>%
    group_by(Phylum) %>%
    mutate(n_genera=n()) %>%
    arrange(desc(n_samples), desc(relative_abundance_mean))

write_delim(genera_phyla_stats,"results/genera_phyla_stats.tsv",delim="\t")

################################ save network format ############################

network_genera_community_matrix <- community_matrix_l %>%
    pivot_wider(id_cols=c(ENA_RUN),
                names_from=scientificName,
                values_from=reads_srs_sum,
                values_fill=0) %>%
    column_to_rownames("ENA_RUN")

write_delim(network_genera_community_matrix,"results/network_genera_community_matrix.tsv",delim="\t")

##################################################################################

print("biodiversity script finished")
