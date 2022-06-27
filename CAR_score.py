#!/usr/bin/env python
## -*- coding: utf-8 -*-

"""
Description:
author:    Luis AE Nagai
data:      2020-06-15
aim:       Provide a list of important genes from several centrality measures
input:     Centrality measure list of genes (columns: CM, rows: ranked genes)
output:    List of ranked genes by CAR score

last modified: 2022-06-23

TODO
1- 


"""


# ---------------------------------------------------------------
# Libraries
# ---------------------------------------------------------------
import argparse
import pandas as pd
from pathlib import Path


# ---------------------------------------------------------------
# Input
# ---------------------------------------------------------------
# For testing, lets just use a static file.  
#df1 = pd.read_csv("minus_allCMs_genes.csv")
parser = argparse.ArgumentParser(description="Calculate Centrality Measure Ranking scores")

# Instructions for the matrix input
parser.add_argument('--input', type=str,
                    help="A matrix of genes ranked by centralities. Column names should contain the centralities")

# Instructions for the name
parser.add_argument('--filename', type=str,
                    help="Insert the sample name for this experiment.")

# Instructions for the output
parser.add_argument('--outputdir', type=str,
                    help="Include the path for the output.")

args = parser.parse_args()



# ---------------------------------------------------------------
# Auxiliar functions
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
# Main CAR score function
# ---------------------------------------------------------------
def getCAR(df_input, file_name, output_dir):

    #output name
    suffix = "_CARscore.csv"
    suffix_gene = "_CARscore_genes.csv"
    outputfile_gene = Path(output_dir, file_name + suffix_gene)
    outputfile = Path(output_dir, file_name + suffix)
    print(outputfile)

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
    # on the centralities, include more weights. The more the gene appers
    # the better ranked it will be.
    car2 = {}
    centralitity_weight = 4
    car2 = {gene: ((df_input.values == gene).sum()) ** centralitity_weight for gene in gene_names}
    #car2
    

    ### Step3 
    # Divide each gene's sum of weights by the number of centralities identified with weights
    car3 = {}
    car3_sorted = {}
    car3 = {x:float(car[x])/car2[x] for x in gene_names}
    car3_sorted = sorted(car3.items(), key=lambda x: x[1])
    car3_sorted = pd.DataFrame(car3_sorted)
    car3_sorted.columns = ['genes', 'CAR.scores']
    

    ### Step4
    # Export a csv file with only names. It can be used in GO term enrichment
    car3_sorted.genes.to_csv(outputfile_gene, index=False, header=False)
    
    # Export the CAR score file 
    car3_sorted.to_csv(outputfile, index=False, header=True)
    
    # Return the CAR score
    return car3_sorted



if __name__ == '__main__':

    # Read input file
    df1 = pd.read_csv(args.input)

    # Calculate the CAR score
    result = getCAR(df1, args.filename, args.outputdir)
    
    # Print the whole list on the screen
    print(result.to_string())