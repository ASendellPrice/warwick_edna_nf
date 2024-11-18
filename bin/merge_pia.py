#!/usr/bin/env python

"""
File:           merge_pia.py
Author:         Ashley T. Sendell-Price
Date:           23.06.2024
Description:    Takes output from multiple runs of PIA (https://github.com/Allaby-lab/PIA) 
                and merges into a single tab delimited file with taxonomic information added
                e.g. kingdom, superphylum, phylum etc. Taxonomic info sourced using taxaranks
                tool from https://github.com/linzhi2013/taxonomy_ranks
Usage:          python merge_pia.py PATH/TO/PIA/OUTPUT/DIRECTORY (requires blast_edna env)
                run from PIA output folder for primer
"""

#Import required modules
import pandas as pd
import os

#Create a list of PIA "Summary_Basic.txt" files  
keyword = '_Basic.txt'
current_dir = os.getcwd()  # Get the current working directory
files = [f for f in os.listdir(current_dir) if os.path.isfile(f) and keyword in f]

#Load first PIA file in list as pandas df 
sampleID = files[0].split(".")[0]
df = pd.read_csv(files[0], sep='\t', skiprows=11)
df = df.rename(columns={"Reads" : sampleID})

#Load subsequent PIA files and merge
for i in range(1,len(files)-1):
    sampleID = files[i].split(".")[0]
    df2 = pd.read_csv(files[i], sep='\t', skiprows=11)
    df2 = df2.rename(columns={"Reads" : sampleID})
    df = pd.merge(df, df2, how = 'outer')
    
#Replae NAs with zero and rename column 1 to "taxa_ID"
df = df.fillna(0)
df = df.rename(columns={"# ID": "taxa_ID"})

# Add taxonomic information
df["taxa_ID"].to_csv('taxaIDs.txt', sep ='\t', index=False, header=False)
os.system('taxaranks -i taxaIDs.txt -o taxa_info.txt')
taxa = pd.read_csv('taxa_info.txt', sep='\t').rename(columns={"user_taxa": "taxa_ID"})
df = pd.merge(df, taxa, how = 'outer')

# Write dataframe to file
df.to_csv("pia_merged_incl_taxa.txt", sep ='\t', index=False)
