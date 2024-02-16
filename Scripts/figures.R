#!/usr/bin/env Rscript
###############################################################################
# script name: figures.R
# developed by: Savvas Paragkamian
# framework: ISD Crete
###############################################################################
# GOAL:
# Aim of this script is to load the results and with some minor counts and 
# transformations make publication-quality figures. That includes maps, bar 
# plots, scatterplots etc.
###############################################################################
# OUTPUT:
# 1. Maps
#   a) blank
#   b) elevation
# 2. Taxonomy: 
#   a) distribution
#   b) ratios
#   c) representative taxa
#   d) prevalence
# 3. Numerical ecology
#   a) diversity indices
#   b) beta diversity
#   c) ordination
#
###############################################################################
# time: takes ~1 hour because of the many maps produced
###############################################################################
# usage:./Scripts/figures.R
###############################################################################

# load packages and functions
#setwd("../")
source("Scripts/functions.R")
library(vegan)
library(tidyverse)
library(ggnewscale)
library(ggpubr)
library(pheatmap)
library(sf)
library(jpeg)
library(raster)
library(scales)
library(ucie)

################################## Load data ##################################
## biodiversity
community_matrix_l <- read_delim("Results/community_matrix_l.tsv",delim="\t")
crete_biodiversity <- read_delim("Results/crete_biodiversity_asv.tsv",delim="\t")
asv_metadata <- read_delim("Results/asv_metadata.tsv",delim="\t")
genera_phyla_samples <- read_delim("Results/genera_phyla_samples.tsv",delim="\t")

# Metadata
metadata <- read_delim("Results/sample_metadata.tsv", delim="\t")

metadata$elevation_bin <- factor(metadata$elevation_bin,
                        levels=unique(metadata$elevation_bin)[order(sort(unique(metadata$elevation_bin)))])

samples_ucie_nmds_genera <- read_delim("Results/samples_ucie_nmds_genera.tsv")

## spatial
locations_spatial <- metadata %>%
    st_as_sf(coords=c("longitude", "latitude"),
             remove=F,
             crs="WGS84")

crete_shp <- sf::st_read("spatial_data/crete/crete.shp")
crete_peaks <- read_delim("spatial_data/crete_mountain_peaks.csv", delim=";", col_names=T) %>%
    st_as_sf(coords=c("X", "Y"),
             remove=F,
             crs="WGS84")

# raster DEM hangling
dem_crete <- raster("spatial_data/dem_crete/dem_crete.tif")
dem_crete_pixel <- as(dem_crete, "SpatialPixelsDataFrame")
dem_crete_df <- as.data.frame(dem_crete_pixel) %>% filter(dem_crete>0)

####################### UCIE #########################
# UCIE needs 3 axis of ordination
print("starting UCIE")

############################### ucie with umap ########################
# sites
#umap_isd_sites_k3 <- read_delim("results/umap_samples_3.tsv", delim="\t") %>% column_to_rownames("id")

#umap_isd_genera_k3 <- read_delim("results/umap_genera_3.tsv", delim="\t") %>% column_to_rownames("id")

#umap_isd_sites_ucie <- ucie::data2cielab(umap_isd_sites_k3, LAB_coordinates = F)
#colnames(umap_isd_sites_ucie) <- c("ENA_RUN","UCIE")
#write_delim(umap_isd_sites_ucie,"results/umap_isd_sites_ucie.tsv", delim="\t")
# taxa
#umap_isd_taxa_ucie <- ucie::data2cielab(umap_isd_genera_k3)
#colnames(umap_isd_taxa_ucie) <- c("scientificName","UCIE")
#write_delim(umap_isd_taxa_ucie,"results/umap_isd_taxa_ucie.tsv", delim="\t")

############################### ucie with pcoa ########################

pcoa_isd_sites <- read_delim("Results/ordination_pcoa_bray_sites.tsv", delim="\t") %>%
    column_to_rownames("ENA_RUN") %>% dplyr::select(Axis.1, Axis.2, Axis.3)

#pcoa_isd_sites_ucie <- ucie::data2cielab(pcoa_isd_sites, LAB_coordinates = F)
#colnames(pcoa_isd_sites_ucie) <- c("ENA_RUN","UCIE")
#write_delim(pcoa_isd_sites_ucie,"results/pcoa_isd_sites_ucie.tsv", delim="\t")
pcoa_isd_sites_ucie <- read_delim("Results/pcoa_isd_sites_ucie.tsv")
############################### ucie with pcoa ########################
#nmds_isd_taxa_ucie <- data2cielab(nmds_isd_taxa_k3, Wb=1.2, S=1.6)
#colnames(nmds_isd_taxa_ucie) <- c("scientificName","UCIE")
#write_delim(nmds_isd_taxa_ucie,"results/nmds_isd_taxa_ucie.tsv", delim="\t")

############################# ucie to location data ####################
locations_spatial <- locations_spatial %>%
    left_join(pcoa_isd_sites_ucie,by=c("ENA_RUN"="ENA_RUN"))

locations_spatial$UCIE[is.na(locations_spatial$UCIE)] <- "gray" 
############################### Color palettes ##################################
okabe_ito_colors <- palette.colors(palette = "Okabe-Ito")   


##################################### 1. MAPS #####################################
# Colorblind palette
palette.colors(palette = "Okabe-Ito")
# Crete figures
cols=c("chocolate1","cornflowerblue","darkgoldenrod1", "darkolivegreen4", "darkorchid1", "goldenrod3", "palevioletred3", "peachpuff4", "turquoise","skyblue")

print("1. Maps")
print("printing base maps")
crete_dem <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_raster(dem_crete_df, mapping=aes(x=x, y=y, fill=dem_crete))+
    scale_fill_gradientn(guide = guide_colourbar(barwidth = 0.5, barheight = 3.5,
                                  title="elevation",
                                  direction = "vertical",
                                  title.vjust = 0.8),
                        colours = c("snow3","#f0e442","#d55e00","#cc79a7"),
                        breaks = c(100, 800, 1500, 2400),
                        labels = c(100, 800, 1500, 2400))+
    coord_sf(crs="wgs84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          line = element_blank(),
          axis.text=element_blank(),
          legend.text=element_text(size=8),
          legend.title = element_text(size=8))

crete_blank <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_point(locations_spatial,
            mapping=aes(x=longitude, y=latitude, color=UCIE),
            size=6,
            show.legend=F) +
    geom_jitter(width = 0.25, height = 0.25)+
    scale_color_manual(values=locations_spatial$UCIE, guide="none")+
    coord_sf(crs="wgs84") +
    theme_bw()+
    theme(
#        panel.background = element_rect(fill='transparent'), #transparent panel bg
        panel.border = element_blank(),
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill='transparent'), #transparent legend bg
        legend.box.background = element_rect(fill='transparent'), #transparent legend panel
        line = element_blank(),
        axis.title=element_blank(),
        axis.text=element_blank(),
        legend.text=element_text(size=8),
        legend.title = element_text(size=8),
        legend.position = "bottom")


ggsave("Figures/map_crete_blank.png",
       plot=crete_blank,
       bg='transparent',
       height = 30,
       width = 60,
       dpi = 300,
       units="cm",
       device="png")


crete_base <- ggplot() +
    geom_sf(crete_shp, mapping=aes()) +
    geom_raster(dem_crete_df, mapping=aes(x=x, y=y, fill=dem_crete))+
#    geom_point(locations_spatial,
#            mapping=aes(x=longitude, y=latitude, color=UCIE, shape=as.character(route)),
#            size=1.7,
#            alpha=0.6,
#            show.legend=T) +
#    geom_jitter(width = 0.25, height = 0.25)+
#    geom_sf(crete_peaks,
#            mapping=aes(),
#            colour=NA,
#            size=1,
#            alpha=1,
#            show.legend=FALSE) +
    geom_label(data = crete_peaks,
               mapping=aes(x = X, y = Y, label = name),
               size = 2,
#               nudge_x = 0.07,
#               nudge_y=0.07,
               label.padding = unit(0.1, "lines"))+
    scale_fill_gradientn(guide = guide_colourbar(barwidth = 0.5, barheight = 3.5,
                                  title="elevation",
                                  direction = "vertical",
                                  title.vjust = 0.8),
                        colours = c("snow3","#f0e442","#d55e00","#cc79a7"),
                        breaks = c(100, 800, 1500, 2400),
                        labels = c(100, 800, 1500, 2400))+
#    scale_color_manual(values=locations_spatial$UCIE, guide="none")+
#    scale_shape_manual(values=c(seq(0,9,1)),name="route")+
    coord_sf(crs="wgs84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          panel.border = element_blank(),
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.background = element_rect(fill='transparent'), #transparent legend bg
          line = element_blank(),
          axis.text=element_blank(),
          legend.text=element_text(size=8),
          legend.title = element_text(size=8),
          legend.position = c(0.95,0.8),
          legend.box.background = element_blank())


ggsave("Figures/map_fig1a.tiff",
       plot=crete_base,
       height = 10,
       width = 20,
       dpi = 300,
       units="cm",
       device="tiff")

ggsave("Figures/map_fig1a.png",
       plot=crete_base,
       height = 10,
       width = 20,
       dpi = 300,
       units="cm",
       device="png")


fig1 <- ggarrange(crete_blank,crete_base,
          labels = c("A", "B"),
          align = "hv",
          heights = c(1,1),
          ncol = 1,
          nrow = 3,
          font.label=list(color="black",size=22)) + bgcolor("white")

ggsave("Figures/map_fig1.tiff", 
       plot=fig1, 
       height = 40, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="tiff")

ggsave("Figures/map_fig1.png", 
       plot=fig1, 
       height = 40, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")

ggsave("Figures/map_fig1.pdf", 
       plot=fig1, 
       height = 40, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="pdf")

ggsave("Figures/map_fig1-small.png", 
       plot=fig1, 
       height = 40, 
       width = 30,
       dpi = 300, 
       units="cm",
       device="png")


################################# 2. Taxonomy #####################################
print("2. Taxonomy")
## Phyla distribution, average relative abundance and ubiquity
## Biogeography of soil bacteria and archaea across France
#phyla_samples_summary <- read_delim("results/phyla_samples_summary.tsv",delim="\t")

total_samples <- length(unique(community_matrix_l$ENA_RUN))

phyla_samples_summary <- community_matrix_l %>%
    group_by(ENA_RUN,Phylum) %>%
    summarise(asvs=sum(asvs),
              reads_srs_mean=mean(reads_srs_mean),
              reads_srs_sum=sum(reads_srs_sum), .groups="keep") %>%
    group_by(ENA_RUN) %>%
    mutate(relative_srs=reads_srs_sum/sum(reads_srs_sum)) %>%
    mutate(z_srs=(reads_srs_sum-mean(reads_srs_sum))/sd(reads_srs_sum)) %>%
    left_join(dplyr::select(metadata, ENA_RUN, LABEL1, LABEL2,LABEL3, elevation, elevation_bin, geology_na), by=c("ENA_RUN"="ENA_RUN"))
    
phyla_stats <- phyla_samples_summary %>% 
    group_by(Phylum) %>%
    summarise(n_samples=n(),
              total_asvs=sum(asvs),
              total_reads_srs=sum(reads_srs_sum),
              proportion_sample=n_samples/total_samples,
              average_relative=mean(relative_srs)) %>%
    arrange(desc(average_relative)) |>
    mutate(Phylum=fct_reorder(Phylum,proportion_sample, .desc=T))

### distribution phyla samples
distribution_phyla_samples <- ggplot(phyla_stats,
       aes(x = proportion_sample, y = Phylum)) +
  geom_point(size = 2) +  # Use a larger dot
  scale_y_discrete(limits=rev) +
  xlab("Samples proportion")+
  ylab("Phylum")+
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed")
  )

ggsave("Figures/taxonomy_distribution_phyla_samples.png",
       plot=distribution_phyla_samples,
       device="png",
       height = 20,
       width = 23,
       units="cm")
#####
taxa_samples_summary <- community_matrix_l %>%
    group_by(ENA_RUN,scientificName, classification) %>%
    summarise(asvs=sum(asvs),
              reads_srs_mean=mean(reads_srs_mean),
              reads_srs_sum=sum(reads_srs_sum), .groups="keep") %>%
    group_by(ENA_RUN) %>%
    mutate(relative_srs=reads_srs_sum/sum(reads_srs_sum))

n_taxa <- length(unique(taxa_samples_summary$classification))
okabe_ito_colors <- palette.colors(palette = "Okabe-Ito")   
fill_colors <- colorRampPalette(okabe_ito_colors)(n_taxa)

taxa_ratios_bar <- ggplot() + 
    geom_col(taxa_samples_summary, mapping=aes(x=ENA_RUN,
                                                y=relative_srs,
                                                fill=classification)) +
    scale_fill_manual(values=fill_colors) +
    theme_bw()+
    coord_flip()+
    scale_y_continuous(expand = c(0, 0), name="")+
    theme(axis.text.x = element_text(face="bold",
                                     size = 15),
          axis.ticks.y=element_blank(),axis.title=element_blank(),
          axis.text.y=element_text(size=10, hjust=0, vjust=0)) +
    theme(legend.position="bottom",
            panel.border = element_blank(),
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank())+
    guides(fill=guide_legend(nrow=1,byrow=TRUE))

ggsave("Figures/taxonomy_ratios_taxa_samples.png",
       plot=taxa_ratios_bar,
       device="png",
       height = 55,
       width = 38,
       units="cm")

############################### phyla ratios matter ###########################
n_phyla <- length(unique(phyla_samples_summary$Phylum))
okabe_ito_colors <- palette.colors(palette = "Okabe-Ito")   
fill_colors <- colorRampPalette(okabe_ito_colors)(n_phyla)

phyla_ratios_bar <- ggplot() + 
    geom_col(phyla_samples_summary, mapping=aes(x=ENA_RUN,
                                                y=relative_srs,
                                                fill=Phylum)) +
    scale_fill_manual(values=fill_colors) +
    theme_bw()+
    coord_flip()+
    scale_y_continuous(expand = c(0, 0), name="")+
    theme(axis.text.x = element_text(face="bold",
                                     size = 15),
          axis.ticks.y=element_blank(),axis.title=element_blank(),
          axis.text.y=element_text(size=10, hjust=0, vjust=0)) +
    theme(legend.position="bottom",
            panel.border = element_blank(),
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank())+
    guides(fill=guide_legend(nrow=7,byrow=TRUE))

ggsave("Figures/taxonomy_ratios_phyla_samples.png",
       plot=phyla_ratios_bar,
       device="png",
       height = 55,
       width = 38,
       units="cm")

## phyla ratios and elevation
categories <- c("elevation_bin")

for (cat in categories){

    phyla_samples_summary$category <- phyla_samples_summary[[cat]]

    phyla_samples_metadata <- phyla_samples_summary |>
    #    left_join(metadata, by=c("ENA_RUN"="ENA_RUN")) |>
        group_by(category, Phylum) |>
        summarise(n_samples=n(),
                  reads_sum=sum(reads_srs_sum),.groups="drop") |>
        group_by(category) |>
        mutate(relative_abundance_e=reads_sum/sum(reads_sum))

    phyla_ratios_cat_bar <- ggplot() + 
        geom_col(phyla_samples_metadata, mapping=aes(x=category,
                                                    y=relative_abundance_e,
                                                    fill=Phylum)) +
        scale_fill_manual(values=fill_colors) +
        theme_bw()+
        coord_flip()+
        scale_y_continuous(expand = c(0, 0), name="")+
        theme(axis.text.x = element_text(face="bold",
                                         size = 15),
              axis.ticks.y=element_blank(),axis.title=element_blank(),
              axis.text.y=element_text(size=10, hjust=0, vjust=0)) +
        theme(legend.position="bottom",
                panel.border = element_blank(),
                panel.grid.major = element_blank(), #remove major gridlines
                panel.grid.minor = element_blank())+
        guides(fill=guide_legend(nrow=7,byrow=TRUE))
    
    ggsave(paste0("Figures/taxonomy_ratios_phyla_",cat,".png"),
           plot=phyla_ratios_cat_bar,
           device="png",
           height = 55,
           width = 38,
           units="cm")
}

# top phyla ratios
top_phyla <- phyla_samples_summary |> 
    ungroup() |>
    mutate(topPhyla=if_else(relative_srs < 0.05, "Other",Phylum)) |> 
    mutate(Phylum=fct_reorder(Phylum, relative_srs)) |>
    ungroup()

top_phyla_names <- unique(top_phyla$topPhyla)

n_top_phyla <- length(unique(top_phyla$topPhyla))
okabe_ito_colors <- palette.colors(palette = "Okabe-Ito")   
fill_colors <- colorRampPalette(okabe_ito_colors)(n_top_phyla)

phyla_ratios_top <- ggplot() + 
    geom_col(top_phyla, mapping=aes(x=ENA_RUN,
                                                y=relative_srs,
                                                fill=topPhyla)) +
    scale_fill_manual(values=fill_colors) +
    theme_bw()+
    coord_flip()+
    scale_y_continuous(expand = c(0, 0), name="")+
    theme(axis.text.x = element_text(face="bold",
                                     size = 15),
          axis.ticks.y=element_blank(),axis.title=element_blank(),
          axis.text.y=element_text(size=10, hjust=0, vjust=0)) +
    theme(legend.position="bottom",
            panel.border = element_blank(),
            panel.grid.major = element_blank(), #remove major gridlines
            panel.grid.minor = element_blank())+
    guides(fill=guide_legend(nrow=1,
                             byrow=TRUE,
                             reverse=TRUE))

ggsave("Figures/taxonomy_ratios_top_phyla_samples.png",
       plot=phyla_ratios_top,
       device="png",
       height = 55,
       width = 38,
       units="cm")


############# Representative phyla of Cretan soils ########################

### representative phyla boxplot
representative_phyla <- top_phyla |> 
    filter(Phylum %in% unique(top_phyla$topPhyla))

representative_phyla_box <- ggplot(representative_phyla,
                                   mapping=aes(x=relative_srs,
                                               y=Phylum))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(height = 0.1,width = 0.0001, stat="identity",aes(alpha=0.5), color="gray50")+
    scale_x_continuous(breaks=seq(0,0.5,0.1))+
    xlab("Relative abundance")+
    theme_bw()+
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed"))


ggsave("Figures/taxonomy_representative_top_phyla_box.png",
       plot=representative_phyla_box,
       device="png",
       height = 20,
       width = 23,
       units="cm")

### all phyla boxplot
colors_label2 <- colors_clc_label2_v[names(colors_clc_label2_v) %in% unique(metadata$LABEL2)]

phyla_box <- ggplot(top_phyla,
                    mapping=aes(x=relative_srs, y=Phylum))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(height = 0.1,width = 0.0001, stat="identity",alpha=0.9, aes(color=LABEL2))+
    scale_x_continuous(breaks=seq(0,0.5,0.1))+
    scale_color_manual(values=colors_label2)+
    xlab("Relative abundance")+
    theme_bw()+
    theme(legend.position = c(0.75,0.2),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed"))

ggsave("Figures/taxonomy_representative_phyla_box.png",
       plot=phyla_box,
       device="png",
       height = 20,
       width = 23,
       units="cm")


############################ Heatmap distribution Phyla Samples ############################
phyla_samples_w_z <- phyla_samples_summary %>%
    pivot_wider(id_cols=ENA_RUN,
                names_from=Phylum,
                values_from=relative_srs,
                values_fill=0) %>%
    as.data.frame() %>% 
    column_to_rownames("ENA_RUN")

phyla_samples_w <- phyla_samples_summary |>
    filter(reads_srs_sum>1) |>
    mutate(log_srs=log(reads_srs_sum)) |>
    pivot_wider(id_cols=ENA_RUN,
                names_from=Phylum,
                values_from=log_srs,
                values_fill=0) |>
    as.data.frame() |>
    column_to_rownames("ENA_RUN")

dcols = vegdist(phyla_samples_w, method="bray")
drows = vegdist(t(phyla_samples_w), method="robust.aitchison")

png("Figures/taxonomy_distribution_heatmap_phyla_samples.png",
    res=300,
    width=70,
    height=30,
    unit="cm")

pheatmap(t(phyla_samples_w),
         clustering_distance_rows = drows,
         clustering_distance_cols = dcols,
         color=colorRampPalette(c("white","skyblue", "palevioletred3"))(20))

dev.off()

######################### Taxa distribution vs samples #####################

taxa_sample_dist_t <- community_matrix_l |>
    group_by(scientificName,classification) |>
    summarise(n_samples=n(),n_asvs=sum(asvs),.groups="keep") |>
    group_by(n_samples) |>
    summarise(n_taxa=n()) |>
    mutate(classification="total")

taxa_sample_dist_c <- community_matrix_l %>%
    group_by(scientificName,classification) |>
    summarise(n_samples=n(),n_asvs=sum(asvs),.groups="keep") |>
    group_by(n_samples, classification) %>%
    summarise(n_taxa=n(), .groups="keep")

taxa_sample_dist <- rbind(taxa_sample_dist_t, taxa_sample_dist_c)


taxa_sample_dist_plot <- ggplot() +
    geom_point(taxa_sample_dist_t,
               mapping=aes(y=n_taxa, x=n_samples)) +
    scale_x_continuous(breaks=seq(0,150,10), name="Number of samples") +
    scale_y_continuous(trans='log10', name = "Number of Taxa",
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=13),
          axis.title.x=element_text(face="bold", size=13),
          axis.title.y=element_text(face="bold", size=13),
          legend.position = c(0.88, 0.8))

ggsave("Figures/taxonomy_distribution_samples.png",
       plot=taxa_sample_dist_plot,
       device="png",
       height = 20,
       width = 23,
       units="cm")

taxa_sample_dist_plot_f <- ggplot() +
    geom_point(taxa_sample_dist,
               mapping=aes(y=n_taxa, x=n_samples)) +
    scale_x_continuous(breaks=seq(0,150,10), name="Number of samples") +
    scale_y_continuous(trans='log10', name = "Numer of Taxa",
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=13),
          axis.title.x=element_text(face="bold", size=13),
          axis.title.y=element_text(face="bold", size=13),
          legend.position = c(0.88, 0.8)) +
    facet_wrap(vars(classification),scales="fixed")

ggsave("Figures/taxonomy_distribution_samples_facet.png",
       plot=taxa_sample_dist_plot_f,
       device="png",
       height = 20,
       width = 43,
       units="cm")

################################ Phyla and genera #########################
############################## genera prevalence ##########################
n_phyla <- length(unique(genera_phyla_samples$Phylum))
okabe_ito_colors <- palette.colors(palette = "Okabe-Ito")   
fill_colors <- colorRampPalette(okabe_ito_colors)(n_phyla)

genera_phyla_stats <- genera_phyla_samples |>
    group_by(Phylum,Genus) |>
    summarise(n_samples = n(),
              proportion_samples = n_samples/total_samples,
              asvs=sum(asvs),
              relative_abundance_mean=mean(relative_srs),
              relative_abundance_sd=sd(relative_srs), .groups="keep") |>
    arrange(desc(proportion_samples), desc(relative_abundance_mean)) |>
    ungroup()

genera_stat_sample <- ggplot() +
    geom_point(genera_phyla_stats,
               mapping=aes(x=proportion_samples,
                           y=relative_abundance_mean,
                           color=Phylum)) +
#    geom_errorbar(genera_phyla_stats,
#                  mapping=aes(x=proportion_samples,
#                              y=relative_abundance_mean,
#                              ymin=relative_abundance_mean-relative_abundance_sd,
#                              ymax=relative_abundance_mean+relative_abundance_sd,
#                              alpha=0.5, colour=Phylum))+
    scale_color_manual(values=fill_colors)+
    xlab("Samples proportion")+
    scale_y_continuous(trans='log10', name = "Mean relative abundance",
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x))) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=13),
          axis.title.x=element_text(face="bold", size=13),
          axis.title.y=element_text(face="bold", size=13),
          legend.position = "bottom")

ggsave("Figures/taxonomy_prevalence_genera.png",
       plot=genera_stat_sample,
       device="png",
       height = 20,
       width = 35,
       units="cm")

genera_stat_sample_f <- genera_stat_sample + facet_wrap(vars(Phylum))

ggsave("Figures/taxonomy_prevalence_genera_facet.png",
       plot=genera_stat_sample_f,
       device="png",
       height = 50,
       width = 80,
       units="cm")

##################### Community ecology #########################
print("3. Community ecology")
##################### Sample Diversity boxplots #########################
# Metadata to long format for diversity indices
metadata_diversity <- metadata %>%
    pivot_longer(cols=c(asvs,
                        shannon,
                        S.obs,
                        S.chao1,
                        se.chao1,
                        S.ACE,
                        se.ACE,
                        taxa),
                 names_to="diversity",
                 values_to="value") %>%
    as.data.frame()
# metadata abiotic values
metadata_abiotic <- metadata %>%
    pivot_longer(cols=c(latitude,
                        longitude,
                        elevation,
                        total_nitrogen,
                        water_content,
                        total_organic_carbon,
                        carbon_nitrogen_ratio,
                        sample_volume_or_weight_for_DNA_extraction,
                        DNA_concentration),
                 names_to="abiotic_metadata",
                 values_to="value") %>%
    as.data.frame()
# metadata bioclim
metadata_bioclim <- metadata %>%
    pivot_longer(cols = grep("bio_*", colnames(.)),
                 names_to="bioclim_metadata",
                 values_to="value") %>%
    as.data.frame()

#################### Sample diversity gradients ######################
## similar to Structure and function of the global topsoil microbiome but
## without the statistic test.

# Numerical variables to plot against diversity indices

vars <- c("latitude",
          "longitude",
          "elevation",
          "total_nitrogen",
          "water_content",
          "carbon_nitrogen_ratio",
          "total_organic_carbon",
          "sample_volume_or_weight_for_DNA_extraction",
          "DNA_concentration",
          "route")

for (var in vars){

    # bioclimatic variables boxplots
    gradient_scatterplot(metadata_bioclim, var, "value", "bioclim_metadata")
}

# add the bioclim variables to the variables
bioclim <- grep("bio.*", colnames(metadata_diversity),value=T)
vars <- c(vars, bioclim)
for (var in vars){
    
    gradient_scatterplot(metadata_diversity, var, "value", "diversity")
}

# abiotic pairs scatterplots
abiotic <- c("total_nitrogen",
             "water_content",
             "total_organic_carbon",
             "carbon_nitrogen_ratio",
             "sample_volume_or_weight_for_DNA_extraction",
             "DNA_concentration")

abiotic_comb <- as.data.frame(t(combn(abiotic,2)))

gradient_scatterplot(metadata_diversity,"total_nitrogen", "total_organic_carbon", "elevation_bin") 
gradient_scatterplot(metadata_diversity,"total_nitrogen", "total_organic_carbon", "LABEL1") 
gradient_scatterplot(metadata_diversity,"total_nitrogen", "total_organic_carbon", "LABEL2") 

# Categorical variables to plot against diversity indices
cats <- c("vegetation_zone",
          "elevation_bin",
          "location")

metadata_diversity <- metadata_diversity %>% filter(!is.na(geology_na))
for (cat in cats){
    # diversity indices boxplots
    diversity_boxplot(metadata_diversity, cat, "value", "diversity")
    # abiotic metadata boxplots
    diversity_boxplot(metadata_abiotic, cat, "value", "abiotic_metadata") 
    # bioclimatic variables boxplots
    diversity_boxplot(metadata_bioclim, cat, "value", "bioclim_metadata")
}

###################### Sample similarity ############################

sample_cooccur_l <- read_delim("results/sample_cooccur_l.tsv", delim="\t") |>
    left_join(metadata, by=c("rowname"="ENA_RUN")) |>
    left_join(metadata, by=c("colname"="ENA_RUN")) |>
    mutate(elevation_difference=abs(elevation.x - elevation.y)) |>
    mutate(site_locations=ifelse(sites.x==sites.y,"same site","differenct site")) |>
    mutate(nitrogen_difference=abs(total_nitrogen.x-total_nitrogen.y))


elevation_difference_g <- ggplot() +
    geom_point(sample_cooccur_l,
               mapping=aes(x=elevation_difference, y=bray)) +
    geom_smooth(method = lm, show.legend=T)+
#    scale_color_manual(values=c("generalists"="#1370A1",
#                                "no classification"="#999999",
#                                "specialists"="#AE6120"),
#                       name="Prevalence")+
    xlab("Elevation difference")+
    ylab("Sample (dis)similarity")+
    scale_x_continuous(breaks=seq(0,2000,250))+
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=13),
          axis.title.x=element_text(face="bold", size=13),
          axis.title.y=element_text(face="bold", size=13),
          legend.position = c(0.9, 0.08))

ggsave("Figures/community_dissimilarity_elevation_difference.png",
       plot=elevation_difference_g,
       device="png",
       height = 20,
       width = 23,
       units="cm")

nitrogen_difference_g <- ggplot() +
    geom_point(sample_cooccur_l,
               mapping=aes(x=nitrogen_difference, y=bray)) +
    geom_smooth(method = lm, show.legend=T)+
#    scale_color_manual(values=c("generalists"="#1370A1",
#                                "no classification"="#999999",
#                                "specialists"="#AE6120"),
#                       name="Prevalence")+
    xlab("Nitrogen difference")+
    ylab("Sample (dis)similarity")+
    #scale_x_continuous(breaks=seq(0,2000,250))+
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size=13),
          axis.title.x=element_text(face="bold", size=13),
          axis.title.y=element_text(face="bold", size=13),
          legend.position = c(0.9, 0.08))

ggsave("Figures/community_dissimilarity_nitrogen_difference.png",
       plot=nitrogen_difference_g,
       device="png",
       height = 20,
       width = 23,
       units="cm")


##### difference between locations 
###

site_locations_dif <- ggplot(sample_cooccur_l, 
                             mapping=aes(x=site_locations, y=bray))+
    geom_boxplot(outlier.shape = NA)+
    geom_jitter(height = 0.00001,width = 0.3, stat="identity",aes(alpha=0.5), color="gray50")+
    scale_y_continuous(breaks=seq(0,1,0.2))+
    xlab("Sample location") +
    ylab("Sample (dis)similarity")+
    theme_bw()+
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(colour = "grey60", linetype = "dashed"))

ggsave("Figures/community_site_locations_dif.png",
       plot=site_locations_dif,
       device="png",
       height = 20,
       width = 23,
       units="cm")

#### same site different locations 
same_site <- sample_cooccur_l %>% filter(site_locations=="same site")

######################### Ordination ############################
######### plots sites ###########

nmds_isd_sites <- read_delim("results/nmds_isd_sites.tsv", delim="\t")
umap_isd_sites <- read_delim("results/umap_samples_2.tsv", delim="\t")
pcoa_isd_sites <- read_delim("results/ordination_pcoa_bray_sites.tsv", delim="\t")
#umap_isd_sites_k1 <- read_delim("results/umap_samples_1.tsv", delim="\t")
#colnames(umap_isd_sites_k1) <- c("id", "UCIE")

ordination_sites <- nmds_isd_sites %>%
    left_join(umap_isd_sites, by=c("ENA_RUN"="id")) %>%
    left_join(pcoa_isd_sites_ucie) %>%
    left_join(pcoa_isd_sites) %>% 
    left_join(metadata)
#    left_join(umap_isd_sites_k1 ,by=c("ENA_RUN"="id"))
ordination_sites$elevation_bin <- factor(ordination_sites$elevation_bin,
                        levels=unique(ordination_sites$elevation_bin)[order(sort(unique(ordination_sites$elevation_bin)))])

# Categorical variables to plot against ordination
cats <- c("vegetation_zone",
          "location")

for (i in cats){
    print(i)

    if (!is.numeric(ordination_sites[[i]]) & i!="ENA_RUN"){
        print(i)
        ordination_sites_plot(ordination_sites, i,"NMDS1","NMDS2", "nmds",i)
        ordination_sites_plot(ordination_sites, i,"UMAP1","UMAP2", "umap",i)
        ordination_sites_plot(ordination_sites, i,"Axis.1","Axis.2", "pcoa",i)

    }else{
        next
    }
}
# ordination of sites and their 2 locations. Are there differences? 
ordination_sites_plot(ordination_sites,"location","NMDS1","NMDS2","nmds_site_loc","sites")
### 3d for inspection of the dimention reduction
# 
#library(plotly)
#plot_ly(x=umap_isd_sites_k3$UMAP1, y=umap_isd_sites_k3$UMAP2, z=umap_isd_sites_k3$UMAP3, type="scatter3d", mode="markers")

########################## Ordination and Boxplots ###########################

boxplot_single(ordination_sites, "UMAP1", "LABEL2", "elevation_bin")
boxplot_single(ordination_sites, "elevation_bin", "UMAP2", "elevation_bin")
boxplot_single(ordination_sites, "UMAP1", "geology_na", "elevation_bin")
boxplot_single(ordination_sites, "UMAP1", "LABEL3", "shannon")
ordination_sites_plot(ordination_sites, "elevation_bin","UMAP1","UMAP2", "umap","elevation_bin")

gradient_scatterplot(ordination_sites,"UMAP1", "shannon", "none") 
######### plots genera ###########

nmds_isd_taxa <- read_delim("results/nmds_isd_taxa.tsv", delim="\t")
nmds_isd_taxa_ucie <- read_delim("results/nmds_isd_taxa_ucie.tsv", delim="\t")
nmds_isd_taxa <- nmds_isd_taxa |>
    left_join(nmds_isd_taxa_ucie) |>
    left_join(umap_isd_sites, by=c("scientificName"="id"))

nmds_genera_plot <- ggplot() +
    geom_point(data=nmds_isd_taxa,
               mapping=aes(x=NMDS1, y=NMDS2,color=UCIE),show.legend = F) +
    scale_color_manual(values=nmds_isd_taxa$UCIE)+
    coord_equal() +
    theme_bw()

ggsave("Figures/ordination_nmds_genera_plot.png",
       plot=nmds_genera_plot,
       device="png",
       height = 20,
       width = 23,
       units="cm")

nmds_genera_plot_f <- nmds_genera_plot + facet_wrap(vars(Phylum),
                                                    scales="fixed")

ggsave("Figures/ordination_nmds_genera_plot_f.png",
       plot=nmds_genera_plot_f,
       device="png",
       height = 50,
       width = 53,
       units="cm")

print("all Figures were generated")
