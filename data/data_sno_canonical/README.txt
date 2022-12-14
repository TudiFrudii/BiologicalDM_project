This file refers to the canonical interactions that we know from literature between the snoRNAs that are subject of our study and any other genes.
The file used snoDB2 as reference to get the canonical interaction. In particular, we selected columns named 

snrna_targets, lncrna_targets, protein_coding_targets, snorna_targets, mirna_targets, trna_targets, ncrna_targets, pseudogene_targets, other_targets

in order to get a list of interactors of interest. Note that the column "rrna_targets" has not been taken into account because would have no use in our study since the dataset of the expansion is ribodepleted.

From this list of canonical interactors we looked for the ENSAMBLE id on several souces that include
 - snoDB_All_V2.0.tsv
 - spliceosome.tsv
 - ensmbl_to_genomes.tsv
Still, some canonical genes miss a specific ENSAMBLE.

The structure of the file consists in 397 rows x 7 columns, the column describe the follogwing fiels:
1) index
2) ENSABMLE id of the snoRNA of interest
3) number of canonical interactors with standard name
4) list of canonical interactors with standard name
5) number of ENSAMBLE id of canonical interactors
6) list of ENSAMBLE id of canonical interactors
7) map of canonical interactors name to all the ENSAMBLE ids associated to it
