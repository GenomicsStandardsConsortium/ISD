#!/usr/bin/env Rscript

###############################################################################
# script name: functions.R
# developed by: Savvas Paragkamian
# framework: ISD Crete
###############################################################################
# GOAL:
# R Functions of the isd-crete scripts
###############################################################################
# source("scripts/functions.R")
###############################################################################

library(vegan)
#library(tidyverse)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

################################## functions ##################################
# this function keeps the last occurrence of a string separated by |
keep_last <- function(x) tail(strsplit(x, split="; ")[[1]],1)
## create dir
create_dir <- function(dir_name){

    if (!dir.exists(dir_name)){
    dir.create(dir_name)
    print(paste0(dir_name, " directory created", sep=""))
    } else{
        print(paste0(dir_name, " directory exists", sep=""))
    }
}

## transform a wide to long format
dist_long <- function(x,method){
    method <- method
    df <- as.data.frame(as.matrix(x)) |>
    rownames_to_column() |>
    pivot_longer(-rowname,
                 values_to=method,
                 names_to="colname")

    return(df)
}

boxplot_single <- function(dataset, x_axis, y_axis, color){
    
    plotname <- paste0("figures/",
                       "ordination_",
                       x_axis,
                       "_",
                       y_axis,
                       "_boxplot.png")
    x_lab <- x_axis
    y_lab <- y_axis
    dataset <- as.data.frame(dataset)
    dataset$x_axis <- dataset[,x_axis]
    dataset$y_axis <- dataset[,y_axis]
    dataset$color <- dataset[,color]

    boxplot_g <- ggplot(dataset, mapping=aes(x=x_axis, y=y_axis)) +
        geom_boxplot()+
        geom_jitter(mapping=aes(color=color),width = 0.2)+
        xlab(x_lab)+
        ylab(y_lab) +
        theme_bw()+
        theme(axis.text.x = element_text(face="bold",
                                         size = 10,
                                         angle = 0,
                                         vjust = 1,
                                         hjust=1))
    ggsave(plotname, 
           plot=boxplot_g, 
           device="png", 
           height = 25, 
           width = 30, 
           units="cm")
}

## function
diversity_boxplot <- function(dataset, x_axis, y_axis, grouping_var){
    plotname <- paste0("figures/",
                       grouping_var,
                       "_",
                       x_axis,
                       "_boxplot.png")
    x_lab <- x_axis
    y_lab <- y_axis

    dataset <- as.data.frame(dataset)
    dataset$x_axis <- dataset[,x_axis]
    dataset$y_axis <- dataset[,y_axis]
    dataset$grouping_var <- dataset[,grouping_var]

    box_diversity <- ggplot(dataset, mapping=aes(x=x_axis, y=y_axis)) +
        geom_boxplot()+
        geom_jitter(width = 0.2)+
        xlab(x_lab)+
        ylab(y_lab) +
        theme_bw()+
        theme(axis.text.x = element_text(face="bold",
                                         size = 10,
                                         angle = 45,
                                         vjust = 1,
                                         hjust=1)) +
        facet_wrap(vars(grouping_var), scales = "free")
    
    ggsave(plotname, 
           plot=box_diversity, 
           device="png", 
           height = 45, 
           width = 30, 
           units="cm")
}

gradient_scatterplot <- function(dataset, x_axis, y_axis, grouping_var){
    
    ## the dataset must be a dataframe, not a tibble, the colnames
    ## must characters

    ## keep the character names of column names to pass to plot
    ##
    plotname <- paste0("figures/",
                       y_axis,
                       "_",
                       x_axis,
                       "_",
                       grouping_var,
                       "_gradient.png")
    x_lab <- x_axis
    y_lab <- y_axis
    
    dataset <- as.data.frame(dataset)
    dataset$x_axis <- dataset[,x_axis]
    dataset$y_axis <- dataset[,y_axis]
    
    gradient <- ggplot(dataset, mapping=aes(x=x_axis, y=y_axis)) +
        geom_point()+
        geom_smooth(method = lm, show.legend=T)+
        xlab(x_lab)+
        ylab(y_lab) +
        theme_bw() +
        theme(panel.grid = element_blank(),
              axis.text = element_text(size=13),
              axis.title.x=element_text(face="bold", size=13),
              axis.title.y=element_text(face="bold", size=13),
              legend.position = c(0.88, 0.8))

    if (grouping_var=="none"){

        ggsave(plotname,
               plot=gradient,
               device="png",
               height = 20,
               width = 23,
               units="cm")
    }else{

        dataset$grouping_var <- dataset[,grouping_var]
        gradient2 <- gradient + facet_wrap(vars(dataset$grouping_var), scales = "free")

        ggsave(plotname,
               plot=gradient2,
               device="png",
               height = 20,
               width = 23,
               units="cm")
    }
}


ordination_sites_plot <- function(df,col,x_axis,y_axis,method, color){
    
    shapes <- length(unique(df[[col]]))
    print(shapes)
    sites_plot <- ggplot() +
        geom_point(data=df,
                   mapping=aes(x=.data[[x_axis]],
                               y=.data[[y_axis]],
                               color=.data[[color]],
                               shape=.data[[col]])) +
        scale_color_manual(values=df$UCIE, guide = "none")+
        scale_shape_manual(values=c(seq(0,shapes,1)))+
        coord_fixed() +
        theme_bw()

    ggsave(paste0("figures/ordination_",method,"_sites_plot_",col,".png"),
        plot=sites_plot,
        device="png",
        height = 20,
        width = 23,
        units="cm")
}
