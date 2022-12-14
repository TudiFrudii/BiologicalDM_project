import pandas as pd
import numpy as np
import os
from itertools import chain
import itertools


def remove_duplicates(list):
    new_list = [] 
    [new_list.append(x) for x in list if x not in new_list]
    return new_list

def flatten_list(list):
     return [item for sublist in list for item in sublist]

def get_spliceosome_complete_map(spliceosome):
    big_map = {}
    # iterate through entire df
    for index, row in spliceosome.iterrows():

        alternative_names = []
        if isinstance(row.loc["Previous symbols"], str):
            alternative_names.append(row.loc["Approved symbol"].split(", "))
        if isinstance(row.loc["Previous symbols"], str):
            alternative_names.append(row.loc["Previous symbols"].split(", "))
        if isinstance(row.loc["Alias symbols"], str):
            alternative_names.append(row.loc["Alias symbols"].split(", "))
        
        alternative_names = flatten_list(alternative_names)

        big_map[row.loc["Ensembl gene ID"]] = alternative_names

    return big_map

def get_snoDB_map(dataframe):
    big_map = {}
    # iterate through entire df
    for index, row in dataframe.iterrows():

        alternative_names = []
        if isinstance(row.loc["gene_name"], str):
            alternative_names.append(row.loc["gene_name"].split(";"))
        if isinstance(row.loc["synonyms"], str):
            alternative_names.append(row.loc["synonyms"].split(";"))
        
        alternative_names = flatten_list(alternative_names)

        big_map[row.loc["ensembl_id"]] = alternative_names

    return big_map


# relative path for data
SNO_DB_PATH = "C:/Users/vitto/Documents/python/data_mining_project/data_mining_project_cutoff_definition/snoDB_All_V2.0.tsv"
SPLICEOSOME_PATH = "C:/Users/vitto/Documents/python/data_mining_project/data_mining_project_cutoff_definition/spliceosome.tsv"
SNO_LIST = "C:/Users/vitto/Documents/python/data_mining_project/exp_snorna.txt"
RIBO_LIST = "C:/Users/vitto/Documents/python/data_mining_project/ribosomal_list.txt"
ENS_TO_GENE_NAME = "C:/Users/vitto/Documents/python/data_mining_project/data_mining_project_cutoff_definition/ensembl_to_genenames.tsv"

# reading snoDB file
sno_db = pd.read_csv(SNO_DB_PATH, delimiter="\t", index_col=1)
sno_db_2 = pd.read_csv(SNO_DB_PATH, delimiter="\t")

# reading expanded_sno_list
temp = pd.read_csv(SNO_LIST, sep='\n', header=None)
list_snoRNAs = temp[0]

# reaging data from spliceosome table
spliceosome = pd.read_csv(SPLICEOSOME_PATH, delimiter="\t")

# reading file to associate ens to std gene names
ensg_to_gene_name = pd.read_csv(ENS_TO_GENE_NAME, sep="\t")

# create dataframe to contain canonical interaction of genes of interest
df_column_names = ["gene_of_interest", "len_target_list", "target_list", "len_ensg_target_list", "ensg_target_list", "target_std_name_to_ensg_map"]
genes_interactions = pd.DataFrame(columns=df_column_names)

# columns with target names
# target_types = ["rrna_targets", "snrna_targets", "lncrna_targets", "protein_coding_targets", "snorna_targets",
# "mirna_targets", "trna_targets", "ncrna_targets", "pseudogene_targets", "other_targets"]
target_types = ["snrna_targets", "lncrna_targets", "protein_coding_targets", "snorna_targets",
"mirna_targets", "trna_targets", "ncrna_targets", "pseudogene_targets", "other_targets"]

spliceosome_data_synonyms_map = get_spliceosome_complete_map(spliceosome=spliceosome)
snodb_data_synonyms_map = get_snoDB_map(dataframe=sno_db_2)

genes_without_ensmb = []


# for each expanded snoRNA
for interest_gene in list_snoRNAs:
    target_list = []

    # get interactors from the selected list of columns
    for target_type in target_types:
        # if sno_db.loc[interest_gene, target_type]:
        if isinstance(sno_db.loc[interest_gene, target_type], str):
            target_list.append(sno_db.loc[interest_gene, target_type].split(";"))
    # get a flattened list
    target_list = flatten_list(target_list)


    # fix target list
    for idx, target in enumerate(target_list):
        # avoid portion of th estring after "-" of 18s or 28s genes
        if target.startswith("18S-") or target.startswith("28S-"):
            target_list[idx] = target.split("-")[0]


    # remove duplicates from target_list
    target_list = remove_duplicates(target_list)

    # create map for association std_name with ensg name
    # this map is specific for every gene of interest (snoRNA expanded)
    # it aims to have as key the canonical interactor as std_name and as values 
    # a list of ensg_id associated to that interactor
    interest_gene_interaction_map_to_ensg = {}

    # convert target list into ensg_ids 
    total_ensg_target_list = []

    # for every gene that canonical
    for target in target_list:
        single_target_ensg_list = []
        found = False

        # converting as stefano said
        if target.startswith("U"):
            target = target.split("-")[0]
            if "." in target:
                target = target.split(".")[0] + "-" + target.split(".")[1]

        # list of ensamble ids for genename as main name
        row = (ensg_to_gene_name.loc[ensg_to_gene_name['Gene_name'] == target])
        total_ensg_target_list.append(row['Gene_stable_ID'].to_list())
        single_target_ensg_list.append(row['Gene_stable_ID'].to_list())
        if row['Gene_stable_ID'].to_list() != []:
            found = True
        
        # list of ensamble ids for genename as synonym
        row = (ensg_to_gene_name.loc[ensg_to_gene_name['Gene_Synonym'] == target])
        total_ensg_target_list.append(row['Gene_stable_ID'].to_list())
        single_target_ensg_list.append(row['Gene_stable_ID'].to_list())
        if row['Gene_stable_ID'].to_list() != []:
            found = True

        # looking for gene name in spliceosome dataset
        for entry in spliceosome_data_synonyms_map.values(): 
            # check if the target is present in values (synonyms) of this dictionary, if so, get key (ensg_id)
            if target in entry:
                value = {i for i in spliceosome_data_synonyms_map if spliceosome_data_synonyms_map[i]==entry}
                total_ensg_target_list.append(value)
                single_target_ensg_list.append(value)
                found = True

        # last resort: get ensg from snodb
        for entry in snodb_data_synonyms_map.values(): 
            # check if the target is present in values (synonyms) of this dictionary, if so, get key (ensg_id)
            if target in entry:
                # print("considering", target, "in list", entry)
                value = {i for i in snodb_data_synonyms_map if snodb_data_synonyms_map[i]==entry}
                total_ensg_target_list.append(value)
                single_target_ensg_list.append(value)
                found = True

        # if not found:
        #     genes_without_ensmb.append(target)

        # the snoRNA of inetest interacts with "target" that corresponds to a list of ensg values
        single_target_ensg_list = flatten_list(single_target_ensg_list)
        single_target_ensg_list = remove_duplicates(single_target_ensg_list)
        interest_gene_interaction_map_to_ensg[target] = single_target_ensg_list
 

    # get a flattened list without dubplicates
    total_ensg_target_list = flatten_list(total_ensg_target_list)
    total_ensg_target_list = remove_duplicates(total_ensg_target_list)
    
 
    # create new complete row
    new_row = { "gene_of_interest" : interest_gene, 
                "len_target_list" : len(target_list), 
                "target_list" : target_list, 
                "len_ensg_target_list" : len(total_ensg_target_list),
                "ensg_target_list" : total_ensg_target_list, 
                "target_std_name_to_ensg_map": interest_gene_interaction_map_to_ensg
                }
    # insert row into final dataframe
    genes_interactions = genes_interactions.append(new_row, ignore_index = True)


# genes_without_ensmb = remove_duplicates(genes_without_ensmb)
# print(len(genes_without_ensmb))
# print(genes_without_ensmb)

genes_interactions.to_csv("C:/Users/vitto/Documents/python/data_mining_project/data_mining_project_cutoff_definition/ensg_ensg_interactions_v2.csv")