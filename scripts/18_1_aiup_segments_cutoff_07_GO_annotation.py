#Functional annotation of the indentified proteins with predicted disordered regions 
#Predicted by AIUPred with a modfified cutoff from 0.3 to 0.7 to predict a region as redox-sensitive disordered

import pandas as pd
from goscripts import enrichment_stats, obo_tools, gaf_parser
from contextlib import redirect_stdout


#Create a set of accession numbers from the tsv file created with predicted disordered non-domain regions
#using the AIUPred cutoff 0.7
proteins_of_interest= set(pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/17_aiupred_modified_cutoff_07_segments_predicted_scores.tsv", sep="\t").iloc[:, 0])

#Perform functional enrichment on the proteins of interest:
with open("/home/guest/Internship/results/Interproscan_Pfam/18_annotation_script_output.txt", "w") as file:
    with redirect_stdout(file):
        GO_dict=obo_tools.importOBO("/home/guest/Internship/data/go-basic.obo", ignore_part_of=True)
        gaf_Dict= gaf_parser.importGAF("/home/guest/Internship/data/HUMAN-uniprot.gaf",geneSet=None)
        gaf_Subset= gaf_parser.importGAF("/home/guest/Internship/data/HUMAN-uniprot.gaf", proteins_of_interest)

        functional_annotation_dict=enrichment_stats.enrichmentAnalysis(GOdict=GO_dict,gafDict=gaf_Dict,gafSubset=gaf_Subset)

#Create a file for the result:
functional_annotation_df = pd.DataFrame(functional_annotation_dict)
functional_annotation_df.index.name = "GO_ID"
functional_annotation_df=functional_annotation_df.reset_index()
functional_annotation_df.to_csv("/home/guest/Internship/results/Interproscan_Pfam/18_aiupred_cutoff_07_segments_functional_enrichment_results.tsv", sep="\t", index=False)

#Parsing the accession numbers and protein names from the gaf file:creating a dataframe
gaf_file = pd.read_csv("/home/guest/Internship/data/HUMAN-uniprot.gaf",sep="\t",comment="!",header=None, dtype=str)

proteins_acc_name = gaf_file.iloc[:, [1, 9]].copy()

proteins_acc_name.columns = ["accession", "protein_name"]

#Creating a dictionary which contains the interested proteins' accesion number,name,GO ID-s
proteins_acc_name_GO_dict = {}

for _, row in proteins_acc_name.iterrows():
    acc = row["accession"]
    name = row["protein_name"]
    if acc in gaf_Subset:
        proteins_acc_name_GO_dict[acc] = {
            "name": name,
            "GOs": gaf_Subset[acc]
        }

#Creating a file with the interested acc numbers,protein names and GO numbers:
result_file=open("/home/guest/Internship/results/Interproscan_Pfam/18_aiupred_cutoff_07_segments_GO_annotations.tsv","w")
result_file.write("Accession_number\tProtein_name\tGO_ID\n")

for acc,values in proteins_acc_name_GO_dict.items():
    name=proteins_acc_name_GO_dict[acc]["name"]
    GOs=proteins_acc_name_GO_dict[acc]["GOs"]
    result_file.write(f"{acc}\t{name}\t{GOs}\n")

result_file.close()

