#Collect those proteins that are found both by Iupred2a and Aiupred
import pandas as pd

#Reading in the two filtered tsv files that contains the found proteins:
pfam_proteins_iup_file= pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/09_pfam_iupred_found_proteins.tsv", sep="\t")
pfam_proteins_aiup_file= pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/09_pfam_aiupred_found_proteins.tsv", sep="\t")

#Creating a list that contains the accession numbers of proteins found by Iupred2a:
iup_proteins=[]

for iup_index,iup_row in pfam_proteins_iup_file.iterrows():
    acc_iup=iup_row["Accession_number"]
    if acc_iup not in iup_proteins:
        iup_proteins.append(acc_iup)

#Creating a list that contains the accession numbers of proteins found by Aiupred:
aiup_proteins=[]

for aiup_index,aiup_row in pfam_proteins_aiup_file.iterrows():
    acc_aiup=aiup_row["Accession_number"]
    if acc_aiup not in aiup_proteins:
        aiup_proteins.append(acc_aiup)

#Creating a file for the found common proteins:
num_found_proteins= 0
common_found_proteins=open("/home/guest/Internship/results/Interproscan_Pfam/10_pfam_common_found_proteins.tsv","w")
common_found_proteins.write(f"Accession_number\n")

#Iterate over the two files with found proteins to find common ones:
for acc_iup in iup_proteins:
    for acc_aiup in aiup_proteins:
        if acc_iup == acc_aiup:
            common_found_proteins.write(f"{acc_iup}\n")

common_found_proteins.close()


common_found_proteins_file= pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/10_pfam_common_found_proteins.tsv", sep="\t")
num_found_proteins= len(common_found_proteins_file)

#Create a txt file for the number of found common proteins:
with open("/home/guest/Internship/results/Interproscan_Pfam/10_pfam_number_of_common_proteins.txt","w") as file:
    file.write(f'{num_found_proteins} proteins can be found both by Iupred2a and Aiupred')