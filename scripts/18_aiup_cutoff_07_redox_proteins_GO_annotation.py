#Looking for proteins predicted by AIUPred that have role in biological processes related to redox regulation
#Inspect predicted redox-sensitive proteins by AIUPred with the modified 0.7 cutoff value

import pandas as pd
import ast
from goatools import obo_parser

#Read in the GO annotation file of the predicted redox-sensitive disordered proteins:
GO_file=pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/18_aiupred_cutoff_07_segments_GO_annotations.tsv",sep="\t")

#Create a list from the GO IDs of interest
GO_list=[]

for _,row in GO_file.iterrows():
    GO_IDs=ast.literal_eval(row["GO_ID"])
    for GOs in GO_IDs:
        if GOs not in GO_list:
            GO_list.append(GOs)

#Read in the OBO with GODag:
obo_file=obo_parser.GODag("/home/guest/Internship/data/go-basic.obo")
#print(obo_file)

#Filter the file (object) so it only contains the GO_IDs and their information from the GO_IDs of interest
filtered_obo_file={
    go_id: obo_file[go_id]
    for go_id in GO_list
    if go_id in obo_file
}
#print(filtered_obo_file)

#Create a dictionary from the original obo file containing the GO ID and the definition text
obo_file_definitions = {}

with open("/home/guest/Internship/data/go-basic.obo") as file:
    GO_id = None

    for line in file:
        line = line.strip()

        if line.startswith("id: GO:"):
            GO_id = line[4:]

        elif line.startswith("def: ") and GO_id:
            obo_file_definitions[GO_id] = line[5:]

#Combine the two dictionaries together so it contains all the information about the GOs of interest:

combined_GO_dict={}

for go_id,definition in obo_file_definitions.items():
    if go_id in filtered_obo_file:
        term=filtered_obo_file[go_id]

        combined_GO_dict[go_id]={
            "name": term.name,
            "definition": definition,
            "namespace": term.namespace}

#Create a dictionary containing only those GO IDs and related information that contains the word "redox" in their definition
redox_GO_dict={}

for go_id,values in combined_GO_dict.items():
    definition=combined_GO_dict[go_id]["definition"]
    if "redox" in definition.lower()or "oxidat" in definition.lower():
        redox_GO_dict[go_id]=values

#Go over the GO annotation file, and collect those proteins that has the GO ID-s with the redox functions:
redox_proteins_dict={}

for _,row in GO_file.iterrows():
    acc=row["Accession_number"]
    name=row["Protein_name"]
    GO_IDs=ast.literal_eval(row["GO_ID"])
    for GOs in GO_IDs:
        if GOs in redox_GO_dict:
            definition=redox_GO_dict[GOs]["definition"]
            namespace=redox_GO_dict[GOs]["namespace"]
            redox_proteins_dict[acc]={"name":name,
                                      "relevant GO_ID": GOs,
                                      "namespace":namespace,
                                      "definition": definition,
                                      "GO_IDs":GO_IDs}

#Create a file for the results:
result_file=open("/home/guest/Internship/results/Interproscan_Pfam/18_aiupred_cutoff_07_redox_proteins.tsv","w")
result_file.write("Accession_number\tProtein_name\tRelevant_GO_ID\tGO_category\tDefinition\tOther_GO_IDs\n")


for acc, info in redox_proteins_dict.items():
    name=redox_proteins_dict[acc]["name"]
    relevant=redox_proteins_dict[acc]["relevant GO_ID"]
    namespace=redox_proteins_dict[acc]["namespace"]
    definition=redox_proteins_dict[acc]["definition"]
    GOs=redox_proteins_dict[acc]["GO_IDs"]
    result_file.write(f"{acc}\t{name}\t{relevant}\t{namespace}\t{definition}\t{GOs}\n")

result_file.close()  


