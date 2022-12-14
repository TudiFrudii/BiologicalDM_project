import pandas as pd
import numpy as np
import os

def define_upper_lower_case (gene_name):
    if gene_name.startswith("ensg") or gene_name.startswith("nr"):
        return gene_name.upper()
    else:
        # cases for clusters or trna_genes
        return gene_name

# relative paths for data
INTERACTION_FILES_SNO = "C:/Users/vitto/Documents/python/data_mining_project/nessra_outputs/interactions_files_sno"
INTERACTION_FILES_RIBO = "C:/Users/vitto/Documents/python/data_mining_project/nessra_outputs/interactions_files_ribo"
OUTPUT_DIR = "C:/Users/vitto/Documents/python/data_mining_project/nessra_outputs/merged_interactions_files/single_expansions"

dir_paths = [INTERACTION_FILES_RIBO, INTERACTION_FILES_SNO]

# we are considering either SNO interactions file or RIBO interactions
for path in dir_paths:

    # creating a dictionary for ensg-ensg storing
    # key : [ensg1, ensg2], value : frel of the interaction
    # in this way we do not add twice the same interaction
    pairs_interactors = {}

    print("starting parsing of", path)

    for file in os.listdir(path):
        complete_file_path = os.path.join(path, file)

        # checking if it is a file
        if os.path.isfile(complete_file_path):

            # reading file content
            file_content = pd.read_csv(complete_file_path, skiprows=1, delimiter=",", index_col=0)
            
            # read each row of the file, update or insert values of the dictionary
            for index, row in file_content.iterrows():
                
                # converting to upper or keep it lowe case based on the label
                interactor_x = define_upper_lower_case(str(row.loc["x"]))
                interactor_y = define_upper_lower_case(str(row.loc["y"]))

                # pair of interactors
                interactors = (interactor_x, interactor_y)

                # if the interaction is already inside the keys of the dictionary
                # aka we already found the interaction
                # we compute the average between the current value in the dictionary
                # and the new frel found in the new file
                if interactors in pairs_interactors.keys():
                    # this is a good metric since every pair should compare at most two times 
                    # (in the expansion of x and in the expansions of y)
                    pairs_interactors[interactors] = (pairs_interactors[interactors] + row.loc["Frel"])/2
                else:
                    # otherwise we simply add the interaction and its relative frequency
                    pairs_interactors[interactors] = row.loc["Frel"]
                
    # at this point we have a complete dictionary of the interactions starting from the single
    # expansions (BUT using interaction files for simplicity)

    # creating basic structure for the final dataframe
    column_names = ["x", "y", "frel"]
    dataframe = pd.DataFrame(columns=column_names)

    # converting each entry of the dictionary in a row to add into the final dataframe
    for interactors, frel in pairs_interactors.items():
        new_row = {"x": interactors[0], "y": interactors[1], "frel" : frel}
        dataframe = dataframe.append(new_row, ignore_index = True)
    
    # save dataframe into specific file
    print("\tsaving merged file in", OUTPUT_DIR + "...")
    if path.endswith("ribo"):
        dataframe.to_csv(OUTPUT_DIR + "/merged_interactions_files_ribo.csv")
    else:
        dataframe.to_csv(OUTPUT_DIR + "/merged_interactions_files_sno.csv")
    print("\t...done.")