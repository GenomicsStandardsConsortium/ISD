#!/usr/bin/env Rscript

###############################################################################
# script name: isd_crete_numerical_ecology.R
# developed by: Savvas Paragkamian, Johanna Holms
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is to use Crete biodiversity data of taxonomy and 
# the sample metadata to perform ecological analyses on biodiversity, ordination
# and multivariate comparison.
#
###############################################################################
# OUTPUT:
#
###############################################################################
# usage:./Scripts/isd_crete_numerical_ecology.R Results/otu_community_matrix_l.tsv Results/otu_sample_metadata.tsv "otu"
###############################################################################
# for interactive execution run first
# setwd("../")

source("scripts/functions.R")
library(vegan)
library(ape)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)
library(ggplot2)
library(dendextend) 
################################## User input #################################
#args = c("Results/asv_community_matrix_l.tsv", "Results/asv_sample_metadata.tsv","asv")
# test if there is at least one argument: if not, return an error
args = commandArgs(trailingOnly=TRUE)

if (length(args)<3) {
  stop("Three arguments must be supplied, community matrix long format, sample metadata and a prefix.n", call.=FALSE)
} else if (length(args)>3) {
  stop("Three arguments must be supplied, community matrix long format, sample metadata and a prefix.n", call.=FALSE)
} 

community_matrix_l <- read_delim(args[1], delim="\t")

metadata <- read_delim(args[2], delim="\t")
prefix <- args[3]
################################## Load data ##################################
print("community matrix")

#community_matrix_l <- read_delim("Results/otu_community_matrix_l.tsv", delim="\t")
#metadata <- read_delim("Results/asv_sample_metadata.tsv", delim="\t")
#prefix <- "asv"


######################################community matrix##########################
community_matrix <- community_matrix_l |>
    dplyr::select(ENA_RUN,reads_srs_sum,scientificName) |>
    pivot_wider(names_from=scientificName,
                values_from=reads_srs_sum,
                values_fill=0) |>
    as.data.frame()

write_delim(community_matrix,paste0("Results/",prefix,"_community_matrix.tsv", sep=""), delim="\t")

rownames(community_matrix) <- community_matrix[,1]
community_matrix <- community_matrix[,-1]

# taxonomy
taxa <- community_matrix_l |> distinct(Kingdom,Phylum,Class,Order,Family,Genus,Species,scientificName,classification)
# Metadata

metadata <- metadata |> filter(ENA_RUN %in% rownames(community_matrix))

### 
print("samples with highest values of physicochemical properties")
metadata |> arrange(desc(total_nitrogen)) |> head(n=2) # ERR3697708 , ERR3697732
metadata |>
    dplyr::select(ENA_RUN, source_material_identifiers, total_organic_carbon) |>
    arrange(desc(total_organic_carbon)) |> head(n=10) # ERR3697655, ERR3697675
metadata |> arrange(desc(water_content)) |> head(n=2) ## ERR3697703, ERR3697702 
metadata |>
    dplyr::select(ENA_RUN, source_material_identifiers,shannon, taxa,n_taxonomic_units) |>
    arrange(desc(shannon)) |>
    head(n=2) ## ERR3697693 isd_3_site_4_loc_2, ERR3697675 isd_2_site_5_loc_2 4.57  4.56
metadata |>
    dplyr::select(ENA_RUN, source_material_identifiers,shannon, taxa,n_taxonomic_units) |>
    arrange(desc(n_taxonomic_units)) |> head(n=2) #ERR3697759 isd_7_site_8_loc_2 1075, ERR3697765 isd_7_site_11_loc_2, 1051
metadata |>
    dplyr::select(ENA_RUN, source_material_identifiers,shannon, taxa,n_taxonomic_units) |>
    arrange(desc(taxa)) |>
    head(n=2) # ERR3697703 isd_4_site_5_loc_2 871, ERR3697702 isd_4_site_5_loc_1 869
################################# Metadata correlations #############################
# correlations of diversity and other numerical metadata
with(metadata, cor(shannon, water_content))
with(metadata, cor(taxa, water_content))
with(metadata, cor(n_taxonomic_units, water_content))
with(metadata, cor(shannon, total_organic_carbon))
with(metadata, cor(taxa, total_organic_carbon))
with(metadata, cor(n_taxonomic_units, total_organic_carbon))
with(metadata, cor(shannon, total_nitrogen))
with(metadata, cor(taxa, total_nitrogen))
with(metadata, cor(n_taxonomic_units, total_nitrogen))
with(metadata, cor(shannon, elevation))
with(metadata, cor(taxa, elevation))
with(metadata, cor(n_taxonomic_units, elevation))

with(metadata, cor.test(n_taxonomic_units, total_organic_carbon))
with(metadata, cor.test(taxa, water_content))
with(metadata, cor.test(taxa, elevation)) 

metadata_n <- metadata
rownames(metadata_n) <- metadata$ENA_RUN
nums <- unlist(lapply(metadata_n, is.numeric), use.names = FALSE)
metadata_n <- metadata_n[,c(nums)]
cc <- cor(metadata_n)

cc_sp <- cor(metadata_n, method="spearman")

write.table(cc_sp,
            paste0("Results/",prefix,"_metadata_spearman.tsv"),
            sep="\t",
            row.names=T,
            col.names=NA)
############################## Dissimilarity ###########################
# use the vegan package, the matrix must be transposed
print("(dis)similarities")

########################## Phylum level ########################
## Phyla distribution, average relative abundance and ubiquity
## Biogeography of soil bacteria and archaea across France

#### Community matrix

bray <- vegdist(community_matrix,
                method="bray")
hc <- hclust(bray)
dend <- as.dendrogram(hc) |>
    color_branches(k = 6) |>
    color_labels(k = 6)

png(file=paste0("Figures/", prefix, "_clustering_bray_hclust_samples.png"),
    width = 55,
    height = 20,
    res=300,
    units = "cm",
    bg="white")
plot(dend)
dev.off()

samples_clusters <- as.tibble(as.ggdend(dend)$labels) |>
    dplyr::select(-c(x,y,cex)) 

colnames(samples_clusters) <- c("ENA_RUN", "dend_cluster_color")

write_delim(samples_clusters,paste0("Results/",prefix,"_samples_dend_clusters.tsv"),delim="\t")

## for taxa

bray_tax <- vegdist(t(community_matrix),method="bray")

png(file=paste0("Figures/", prefix, "_clustering_hclust_taxa.png"),
    width = 50,
    height = 50,
    res=300,
    units = "cm",
    bg="white")
plot(hclust(bray_tax))
dev.off()

bray_samples <- vegdist(community_matrix,method="bray")
#homoscedasticity_s <- betadisper(bray_samples, metadata$LABEL1, type = c("median","centroid"), bias.adjust = FALSE)

bray_l <- dist_long(bray, "bray")

jaccard <- vegdist(community_matrix,
                method="jaccard",
                binary=TRUE)

aitchison <- vegdist(community_matrix,
                method="robust.aitchison")

jaccard_l <- dist_long(jaccard, "jaccard")
aitchison_l <- dist_long(aitchison, "robust.aitchison")

# beta diversity
z <- betadiver(community_matrix, "z")
#mod <- with(metadata, betadisper(z, LABEL1))
#sac <- specaccum(biodiversity_srs_t)

######################### Ordination ############################
####################### PCoA #########################
print("starting PCoA")

#### sites
pcoa_bray <- ape::pcoa(bray)
pcoa_bray_m <- pcoa_bray$vectors |> as.data.frame() |> rownames_to_column("ENA_RUN")

write_delim(pcoa_bray_m,
            paste0("Results/",prefix,"_ordination_pcoa_bray_sites.tsv"),
            delim="\t")

####################### nMDS #########################
print("starting nMDS")
nmds_isd <- vegan::metaMDS(community_matrix,
                       k=2,
                       distance = "bray",
                       trymax=100)

# fit environmental numerical vectors
env_isd <- metadata |>
    filter(ENA_RUN %in% rownames(community_matrix)) |> 
    column_to_rownames(var="ENA_RUN")# |>

print("starting envfit")
envfit_isd <- envfit(nmds_isd, env_isd, permutations = 999, na.rm=T) 
env_scores_isd <- as.data.frame(scores(envfit_isd, display = "vectors"))

write_delim(env_scores_isd,
            paste0("Results/",prefix,"_env_scores_isd.tsv"),
            delim="\t")

# plotting
png(file=paste0("Figures/", prefix, "_ordination_nmds_stressplot.png"),
    width = 30,
    height = 30,
    res=300,
    units = "cm",
    bg="white")
stressplot(nmds_isd)
dev.off()

png(file=paste0("Figures/", prefix, "_ordination_nmds_sites_lat.png"),
    width = 30,
    height = 30,
    res=300,
    units = "cm",
    bg="white")
ordiplot(nmds_isd,display="sites", cex=1.25)
ordisurf(nmds_isd,env_isd$latitude,main="",col="firebrick") ## interesting
#ordisurf(nmds,metadata$dem,main="",col="orange")
dev.off()

png(file=paste0("Figures/", prefix, "_ordination_nmds_sites_dem.png"),
    width = 30,
    height = 30,
    res=300,
    units = "cm",
    bg="white")
ordiplot(nmds_isd,display="sites", cex=1.25)
ordisurf(nmds_isd,env_isd$elevation,main="",col="firebrick") ## interesting
dev.off()


nmds_isd_taxa <- as.data.frame(scores(nmds_isd, "species")) |>
    rownames_to_column("scientificName") |>
    left_join(taxa, by=c("scientificName"="scientificName"))

write_delim(nmds_isd_taxa,
            paste0("Results/",prefix,"_nmds_isd_taxa.tsv"),
            delim="\t")

nmds_isd_sites <- as.data.frame(scores(nmds_isd,"sites")) |>
    rownames_to_column("ENA_RUN") |>
    left_join(metadata[,c("ENA_RUN","elevation_bin", "vegetation_zone")],
              by=c("ENA_RUN"="ENA_RUN"))

write_delim(nmds_isd_sites,
            paste0("Results/",prefix,"_nmds_isd_sites.tsv"),
            delim="\t")

############################ nmds k3 ###########################

#nmds_isd_k3 <- vegan::metaMDS(community_matrix,
#                       k=3,
#                       distance = "bray",
#                       trymax=100)
#nmds_isd_taxa_k3 <- as.data.frame(scores(nmds_isd_k3,"species"))
#nmds_isd_sites_k3 <- as.data.frame(scores(nmds_isd_k3,"sites"))
############################# dbRDA ############################

#dbrda_isd <- dbrda(community_matrix ~ elevation + latitude + longitude + total_organic_carbon + total_nitrogen + water_content,env_isd, dist="bray")

############################# UMAP ############################
# the python script isd_crete_umap.py
# performs the UMAP algorithm
print("starting UMAP")
system(paste0("python Scripts/isd_crete_umap.py",
              " ",
              paste0("Results/",prefix,"_community_matrix.tsv", sep=""),
              " ",
              paste0(prefix)),
       wait=TRUE)

print("ordination methods join")
umap_isd_sites <- read_delim(paste0("Results/", prefix,"_umap_samples_2.tsv"), delim="\t")
#umap_isd_sites_k1 <- read_delim("results/umap_samples_1.tsv", delim="\t")
#colnames(umap_isd_sites_k1) <- c("id", "UCIE")

metadata <- metadata |>
    left_join(umap_isd_sites, by=c("ENA_RUN"="id")) |>
    left_join(pcoa_bray_m) |>
    left_join(nmds_isd_sites)
#    left_join(umap_isd_sites_k1 ,by=c("ENA_RUN"="id"))
metadata$elevation_bin <- factor(metadata$elevation_bin,
                        levels=unique(metadata$elevation_bin)[order(sort(unique(metadata$elevation_bin)))])
################################# statistics ##########################

print("starting correlations")

##### regression
## diversity
cor.test(metadata$shannon, metadata$total_nitrogen)
cor.test(metadata$shannon, metadata$carbon_nitrogen_ratio) 

cor.test(metadata$shannon, metadata$elevation)
cor.test(metadata$shannon, metadata$water_content) 

gradient_scatterplot(metadata, "total_organic_carbon","shannon", "elevation_bin", prefix) 

####### Drivers numerical
cor.test(metadata$shannon, metadata$total_nitrogen)
cor.test(metadata$shannon, metadata$total_organic_carbon)
cor.test(metadata$shannon, metadata$carbon_nitrogen_ratio)
cor.test(metadata$shannon, metadata$water_content)
cor.test(metadata$shannon, metadata$elevation)

cor.test(metadata$shannon, metadata$UMAP1)

lm_s <- lm(metadata$shannon ~ metadata$elevation_bin + metadata$vegetation_zone + metadata$total_organic_carbon)
summary(lm_s)
anova(lm_s)

### drivers of major axis of ordination
# first axis
lm_o <- lm(metadata$UMAP1 ~ metadata$elevation_bin + metadata$total_organic_carbon)
summary(lm_o)
anova(lm_o)
gradient_scatterplot(metadata, "total_organic_carbon","UMAP1", "none",prefix) 
gradient_scatterplot(metadata, "total_nitrogen","UMAP1", "none", prefix) 
cor.test(metadata$UMAP1, metadata$total_organic_carbon)
cor.test(metadata$UMAP1, metadata$total_nitrogen)

boxplot_single(metadata, "UMAP1", "elevation_bin", "total_nitrogen", prefix)
# second axis
lm_o2 <- lm(metadata$UMAP2 ~ metadata$total_organic_carbon + metadata$water_content)
summary(lm_o2)
anova(lm_o2)
kruskal.test(UMAP2 ~ elevation_bin, data = metadata)
cor.test(metadata$UMAP2, metadata$total_organic_carbon)
cor.test(metadata$UMAP2, metadata$total_nitrogen)
cor.test(metadata$UMAP2, metadata$water_content)
gradient_scatterplot(metadata, "water_content","UMAP2", "none", prefix) 
gradient_scatterplot(metadata, "total_nitrogen","UMAP2", "none", prefix) 
gradient_scatterplot(metadata, "total_organic_carbon","UMAP2", "none", prefix) 
boxplot_single(metadata, "UMAP2", "total_organic_carbon", "elevation_bin", prefix)
####### Drivers categorical
kruskal.test(shannon ~ elevation_bin, data = metadata)
pairwise.wilcox.test(metadata$shannon, metadata$elevation_bin, p.adjust.method="BH")

########### community dissimilarity tests #############
# calculate the bray dissimilatiry
print("starting community dis-simmilaties")
bray <- vegdist(community_matrix)

# multivariate dispersion (variance) for a group of samples is to calculate
# the average distance of group members to the group centroid or spatial
# median (both referred to as 'centroid' from now on unless stated otherwise)
# in multivariate space. 
## test to see if there are any significant differences 
### Pairwise comparisons of group mean dispersions can also be performed using
### permutest.betadisper. An alternative to the classical comparison of group
### dispersions, is to calculate Tukey's Honest Significant Differences between
### groups, via TukeyHSD.betadisper. This is a simple wrapper to TukeyHSD. The
### user is directed to read the help file for TukeyHSD before using this
### function. In particular, note the statement about using the function with unbalanced designs.

# total nitrogen
mod <- betadisper(bray, metadata$total_nitrogen,type="centroid")
png(file=paste0("Figures/", prefix, "_community_betadisper_nitrogen_box.png"),
    res=300,
    width=60,
    height=40,
    unit="cm")
plot(mod)
dev.off()

anova(mod)

permutest(mod, pairwise = TRUE, permutations = 99)
# elevation
mod <- betadisper(bray, metadata$elevation_bin,type="centroid")

png(file=paste0("Figures/", prefix, "_community_betadisper_elevation_box.png"),
    res=300,
    width=60,
    height=40,
    unit="cm")
boxplot(mod)
dev.off()

anova(mod)

permutest(mod, pairwise = TRUE, permutations = 99)
mod.HSD <- TukeyHSD(mod)
png(file=paste0("Figures/", prefix, "_community_betadisper_elevation_bin.png"),
    res=300,
    width=60,
    height=40,
    unit="cm")
plot(mod.HSD)

dev.off()


#### permanova

print("starting permanova")
#adonis_elevation <- adonis2(community_matrix ~ elevation_bin, data=metadata_f, permutations=99)

adonis_multiple <- adonis2(community_matrix ~ elevation_bin*total_nitrogen*carbon_nitrogen_ratio,
                           data=metadata,
                           permutations=999)

############################## Community analysis ###########################
###################### Co-occurrence of samples and ASVs ####################
print("starting co-occurrence")

#biodiversity_m <- biodiversity_srs
#biodiversity_m[biodiversity_m > 0 ] <- 1
#biodiversity_m <- as.matrix(biodiversity_m)

## matrix multiplication takes up a lot of memory and CPU, I had an error
## Error: vector memory exhausted (limit reached?)
## cd ~ ; touch .Renviron 
## echo R_MAX_VSIZE=200Gb >> .Renviron

#asv_cooccur <- biodiversity_m %*% t(biodiversity_m)
community_matrix_m <- community_matrix
community_matrix_m[community_matrix_m > 0] <- 1
community_matrix_m <- as.matrix(community_matrix_m)

sample_cooccur <- community_matrix_m %*% t(community_matrix_m)
taxa_cooccur <- t(community_matrix_m) %*% community_matrix_m

isSymmetric(taxa_cooccur) # is true so we can remove the lower triangle
taxa_cooccur[lower.tri(taxa_cooccur)] <- NA

taxa_cooccur_l <- dist_long(taxa_cooccur,"cooccurrence") |>
    filter(rowname!=colname) |>
    na.omit()

write_delim(taxa_cooccur_l,
            paste0("Results/",prefix,"_taxa_cooccur_l.tsv"),
            delim="\t")

isSymmetric(sample_cooccur) # is true so we can remove the lower triangle
sample_cooccur[lower.tri(sample_cooccur)] <- NA

sample_cooccur_l <- dist_long(sample_cooccur,"cooccurrence") |>
    filter(rowname!=colname) |>
    na.omit() |> 
    left_join(bray_l,
              by=c("rowname"="rowname", "colname"="colname")) |>
    left_join(jaccard_l,
              by=c("rowname"="rowname", "colname"="colname")) |>
    left_join(aitchison_l,
              by=c("rowname"="rowname", "colname"="colname"))


write_delim(sample_cooccur_l,
            paste0("Results/",prefix,"_sample_cooccur_l.tsv"),
            delim="\t")
######################## Site locations comparison ASV #################
samples_locations <- metadata |>
    pivot_wider(id_cols=sites,
                names_from=location,
                values_from=ENA_RUN)

dissi_loc <- samples_locations |>
    left_join(sample_cooccur_l,
              by=c("loc_1"="rowname", "loc_2"="colname"))

summary(dissi_loc)

print("finish")
