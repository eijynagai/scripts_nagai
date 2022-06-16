#!/usr/bin/env python
## -*- coding: utf-8 -*-

"""
Description:
author:    Luis AE Nagai
data:      2020-06-15
aim:       Provide a list of important genes from several centrality measures
input:     Centrality measure list of genes (columns: CM, rows: ranked genes)
output:    List of important genes ranked

last modified: 2022-06-16

TODO
1- modify the script to parse the matrix of genes
2- modify the main script for CAR as conventional scripts
3- include arguments for the script
    *input file
    *output file name

"""


# ---------------------------------------------------------------
# Libraries
# ---------------------------------------------------------------
import pandas as pd



# ---------------------------------------------------------------
# Input
# ---------------------------------------------------------------
# For testing, lets just use a static file.  
df1 = pd.read_csv("minus_allCMs_genes.csv")



# ---------------------------------------------------------------
# Functions
# ---------------------------------------------------------------

# Transform the dataframe into list
def mat2lst(df):
    return sum(df.values.tolist(), [])

# Remove the duplicated entries from the list
def lst2unq(lst):
    return list(set(lst))

# Input a dataframe and return the gene names 
def getnames(df):
    return lst2unq(mat2lst(df))



# ---------------------------------------------------------------
# Main CAR score funciton
# ---------------------------------------------------------------
def getCAR(df_input, file_name):

    #output name
    filename="%s_CARscore.csv" % file_name
    
    
    ### Step1
    # Create one dictionary having the sum of the positions for each gene

    # CAR score dictonary
    car = {}

    # Keys for the CAR
    #Populate the dictionary with gene names (unique entries only)
    gene_names = []
    gene_names = getnames(df_input)
    car = {gene: 0 for gene in gene_names}
    
    # Clear every row of the matrix
    weight = 0
    # Iterate over the input
    for rowIndex, row in df_input.iterrows(): #iterate over rows
        # Each position in the ranking has a weight. 
        # The less the better
        weight = rowIndex + 1
    
        for columnIndex, value in row.items():
            # For the gene located in the cell, add the weight
            # It should accumulate the with
            car[value] += weight
    

    ### Step 2
    # Create a second dictionary that will contain the gene names 
    # and the counts of genes in centralities. To put more importance
    # on them, multiply it as **2 or **3
    car2 = {}
    car2 = {gene: ((df_input.values == gene).sum())**3 for gene in gene_names}
    car2
    

    ### Step3 
    # Divide each gene sum of weights by the number of centralities identified
    car3 = {}
    car3_sorted = {}
    car3 = {x:float(car[x])/car2[x] for x in gene_names}
    car3_sorted = sorted(car3.items(), key=lambda x: x[1])
    car3_sorted
    

    ### Step4
    # Create a csv file with only names. It can be used in GO term enrichment
    # for quick checking
    ele = []
    for x in car3_sorted:
        ele.append(x[0])
        
    # Save the file into csv files for using in the GO enrichment
    df_final = pd.DataFrame(ele)
    df_final.to_csv(filename, index=False, header=False)
    
    # Return the CAR score
    return car3_sorted