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
# Results/crete_biodiversity_otu.tsv
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
    filter(!is.na(taxonomy)) |>
    rename(taxonomic_unit=OTU)

write_delim(crete_biodiversity_otu, "Results/crete_biodiversity_otu.tsv", delim="\t")

