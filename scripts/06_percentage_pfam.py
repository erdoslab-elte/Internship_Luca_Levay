#Counting how many percentage of the found disordered regions overlaps more than 50% with the annoteted domains
#Pfam

import pandas as pd

#Read the Iupred tsv file
Iupred_tsv=pd.read_csv("/home/guest/Internship/results/01_iupred_disordered_regions.tsv",sep="\t")

#Read the Iupred tsv file
AIUpred_tsv=pd.read_csv("/home/guest/Internship/results/01_aiupred_disordered_regions.tsv",sep="\t")

#Read the created overlap tsv files:
iup_overlap_above_tsv=pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_iupred_overlap_above_50.tsv",sep="\t")
iup_overlap_below_tsv=pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_iupred_overlap_below_50.tsv",sep="\t")
aiup_overlap_above_tsv=pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_aiupred_overlap_above_50.tsv",sep="\t")
aiup_overlap_below_tsv=pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_aiupred_overlap_below_50.tsv",sep="\t")

#Counting how many percentage of the found disordered regions overlaps more than 50% with the annoteted domains
#For Iupred
#count the number of regions:
iup_num_disord_regions=len(Iupred_tsv)

iup_num_overlap_above=len(iup_overlap_above_tsv)
iup_num_overlap_below=len(iup_overlap_below_tsv)
#count the percentage:
iup_percentage_above= iup_num_overlap_above / iup_num_disord_regions *100
iup_percentage_below= iup_num_overlap_below / iup_num_disord_regions *100
#print(iup_percentage_above)
#print(iup_percentage_below)

#For AIUpred
#count the number of regions:
aiup_num_disord_regions=len(AIUpred_tsv)

aiup_num_overlap_above=len(aiup_overlap_above_tsv)
aiup_num_overlap_below=len(aiup_overlap_below_tsv)
#count the percentage:
aiup_percentage_above= aiup_num_overlap_above / aiup_num_disord_regions *100
aiup_percentage_below= aiup_num_overlap_below / aiup_num_disord_regions *100
#print(aiup_percentage_above)
#print(aiup_percentage_below)

#Summarize the results in a txt file
with open("/home/guest/Internship/results/Interproscan_Pfam/06_pfam_percentage_overlapping_regions.txt","w") as file:
    file.write(f"Iupred2a: {iup_percentage_above}% of the found disordered regions overlaps (equal to/)more than 50% with the annotated domains\n")
    file.write(f"Iupred2a: {iup_percentage_below}% of the found disordered regions overlaps less than 50% with the annotated domains\n")
    file.write(f"Aiupred: {aiup_percentage_above}% of the found disordered regions overlaps (equal to/)more than 50% with the annotated domains\n")
    file.write(f"Aiupred: {aiup_percentage_below}% of the found disordered regions overlaps less than 50% with the annotated domains\n")

