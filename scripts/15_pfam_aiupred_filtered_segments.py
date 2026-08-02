#Collecting the found redox-dependent disordered non-domain regions(segments) predicted by Aiupred
#which have a predicted score above 0.5

import pandas as pd
import ast

#Read in the tsv file:
tsv_file=pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_aiupred_segments_predicted_scores.tsv",sep='\t')

#Create a new files for the results:
##For the regions:
region_file=open("/home/guest/Internship/results/Interproscan_Pfam/15_pfam_aiupred_filtered_segments_predicted_scores.tsv", "w")
region_file.write("Accession_number\tStart\tEnd\tPredicted_score\n")

protein_file=open("/home/guest/Internship/results/Interproscan_Pfam/15_pfam_aiupred_filtered_proteins_predicted_scores.tsv", "w")
protein_file.write("Accession_number\tPredicted_score\n")

#Adding data to the result files:
acc_list=[]

n=0
m=0
for _,row in tsv_file.iterrows():
    acc=row["Accession_number"]
    start=row["Start"]
    end=row["End"]
    pred_score=row["Predicted_score"]
    if pred_score >0.5:
        region_file.write(f"{acc}\t{start}\t{end}\t{pred_score}\n")
        n +=1
        if acc not in acc_list:
            acc_list.append(acc)
            protein_file.write(f"{acc}\t{pred_score}\n")
        else:
            continue
        m+=1
    
region_file.close()
protein_file.close()

print(n)
print(m)
