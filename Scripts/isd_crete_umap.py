#!/usr/bin/env python3

###############################################################################
# script name: isd_crete_umap.py
# developed by: Savvas Paragkamian
# framework: ISD Crete 2016
###############################################################################
# GOAL:
# Aim of this script is to apply the umap dimention reduction method 
# to the metagenomic data.
###############################################################################
# usage:./Scripts/isd_crete_umap.py <path/to/community_matrix.tsv> <prefix>
###############################################################################
import pandas as pd
import numpy as np
import os, sys
import umap

## files
if len(sys.argv) < 2:
    print("Usage: python script.py <path/to/community_matrix.tsv> <prefix>")
    sys.exit(1)

# Access the arguments
arg1_value = sys.argv[1]
arg2_value = sys.argv[2] 

if len(sys.argv) > 3: 
    print("Usage: python script.py <arg1> [optional_arg]")
    sys.exit(1)

# Your script logic here
print("arg1:", arg1_value)
print("arg2:", arg2_value)

community_matrix_path = arg1_value  #"Results/community_matrix.tsv"
prefix = arg2_value
output_path = "Results/"

community_matrix = pd.read_csv(community_matrix_path, sep="\t")
community_array = community_matrix.drop(columns=['ENA_RUN'])

# taxa, the colnames to a new dataframe
cols_list = community_array.columns.tolist()
df_taxa = pd.DataFrame(cols_list)
# samples
df_samples = community_matrix[['ENA_RUN']]
community_array = community_array.to_numpy()

# for taxa UMAP the array needs to be transposed 
community_array_t = np.transpose(community_array)

def umap_function(array, df_ids, n_components,axis, prefix):
    
    np.random.seed(123)
    # umap function
    embedding = umap.UMAP(n_neighbors=50,
                          min_dist=0.7,
                          n_components=n_components,
                          metric='braycurtis').fit_transform(array)
    
    #add NumPy matrix as new columns in DataFrame
    df = pd.concat([df_ids, pd.DataFrame(embedding)], axis=1)
    # add colnames
    colnames =['id']
    for i in range(1,n_components+1):
        
        colname = "UMAP" + str(i)
        colnames.append(colname)

    df.columns = colnames
    # create the tsv file
    output_file = output_path + prefix + "_umap_" + axis + "_" + str(n_components) + ".tsv"
    df.to_csv(output_file, header=True, index=None, sep='\t')

    #return(df)

# samples
print("UMAP samples is ongoing")
umap_function(community_array,df_samples, 3, "samples", prefix)
umap_function(community_array,df_samples, 2, "samples", prefix)
umap_function(community_array,df_samples, 1, "samples", prefix)

# taxa
print("UMAP taxa is ongoing")
umap_function(community_array_t,df_taxa,3,"taxa", prefix)
umap_function(community_array_t,df_taxa,2,"taxa", prefix)
umap_function(community_array_t,df_taxa,1,"taxa", prefix)

print("UMAP ordination finished")
