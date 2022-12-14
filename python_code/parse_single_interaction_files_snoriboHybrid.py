import pandas as pd
import numpy as np
import os

# relative paths for data
OUTPUT_DIR = "C:/Users/vitto/Documents/python/data_mining_project/nessra_outputs/merged_interactions_files/single_expansions"
MERGED_INTERACTION_SNO = OUTPUT_DIR + "/merged_interactions_files_sno.csv"
MERGED_INTERACTION_RIBO = OUTPUT_DIR + "/merged_interactions_files_ribo.csv"

files = [MERGED_INTERACTION_RIBO, MERGED_INTERACTION_SNO]

# creating a dictionary for ensg-ensg storing
# key : [ensg1, ensg2], value : frel of the interaction
# in this way we do not add twice the same interaction
pairs_interactors = {}

for file in files:
    print("now considering file", file)
    # reading file content
    file_content = pd.read_csv(file, index_col=0)
    
    for index, row in file_content.iterrows():
        
        # get tuple of interactors
        interactors = (str(row.loc["x"]), str(row.loc["y"]))

        # if the interaction is already inside the keys of the dictionary
        # aka we already found the interaction
        # we compute the average between the current value in the dictionary
        # and the new frel found in the new file

        if interactors in pairs_interactors.keys():
            # this is a good metric since every pair should compare at most two times 
            # (in the expansion of RIBO and in the expansions of SNO)
            pairs_interactors[interactors] = (pairs_interactors[interactors] + row.loc["frel"])/2
        else:
            # otherwise we simply add the interaction and its relative frequency
            pairs_interactors[interactors] = row.loc["frel"]
    
# at this point we have a complete dictionary of the interactions starting from the two merged files
# creating basic structure for the final dataframe
column_names = ["x", "y", "frel"]
dataframe = pd.DataFrame(columns=column_names)

# converting each entry of the dictionary in a row to add into the final dataframe
for interactors, frel in pairs_interactors.items():
    new_row = {"x": interactors[0], "y": interactors[1], "frel" : frel}
    dataframe = dataframe.append(new_row, ignore_index = True)

# save dataframe into specific file
print("\tsaving merged file in", OUTPUT_DIR + "...")
dataframe.to_csv(OUTPUT_DIR + "/merged_interactions_files_sno_ribo.csv")
print("\t...done.")