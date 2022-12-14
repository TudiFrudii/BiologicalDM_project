import pandas as pd
import numpy as np
import os

def lastWord(string):
    newstring = ""
    # calculating length of string
    length = len(string)

    # traversing from last
    for i in range(length-1, 0, -1):
        # if space is occurred then return
        if(string[i] == " "):
            # return reverse of newstring
            return newstring[::-1]
        else:
            newstring = newstring + string[i]

def define_upper_lower_case (gene_name):
    if gene_name.startswith("ensg") or gene_name.startswith("nr"):
        return gene_name.upper()
    else:
        # cases for clusters or trna_genes
        return gene_name

# relative paths for data
INTERACTION_FILES = "C:/Users/vitto/Documents/python/data_mining_project/nessra_outputs/interaction_files_netwoks"
OUTPUT_DIR = "C:/Users/vitto/Documents/python/data_mining_project/nessra_outputs/merged_interactions_files/networks"

# structure of the dataframe
column_names = ["x", "y", "frel"]

# define thresholds of interest
# thresholds = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
thresholds = [0]

for threshold in thresholds:

    for file in os.listdir(INTERACTION_FILES):
        complete_file_path = os.path.join(INTERACTION_FILES, file)

        # checking if it is a file
        if os.path.isfile(complete_file_path):
            
            graph = pd.DataFrame(columns=column_names)

            # get file proper name
            with open(complete_file_path) as f:
                first_line = f.readline()
                file_name = lastWord(first_line)[:-1]
                print(file_name)

            # reading file content
            file_content = pd.read_csv(complete_file_path, skiprows=1, delimiter=",", index_col=0)
            
            # individuate interacting gene with relative frequence, convert into lists
            gene_x = file_content["x"].tolist()
            gene_y = file_content["y"].tolist()
            frel_s = file_content["Frel"].tolist()
            intra = file_content["intra"].tolist()

            # same length of every list
            length = len(gene_x)

            # save actual values into the dataframe
            for i in range(length):
                if frel_s[i] > threshold:
                    x = define_upper_lower_case(str(gene_x[i]))
                    y = define_upper_lower_case(str(gene_y[i]))
                    new_entry = {"x": x, "y" : y, "frel": frel_s[i]}
                    graph = graph.append(new_entry, ignore_index = True)
        
        graph.to_csv(OUTPUT_DIR + "/" + file_name + "_" + str(threshold) + ".csv")