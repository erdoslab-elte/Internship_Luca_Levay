# Modify the files that contains the regions that overlaps more than 50% with the annotated domains:Pfam and Gene3D, Iupred and Aiupred
#Filter on the found domains
import pandas as pd

#Read in the tsv files:
#Pfam
pfam_iup_overlap_file = pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_iupred_overlap_above_50.tsv",sep="\t")
pfam_aiup_overlap_file = pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_aiupred_overlap_above_50.tsv",sep="\t")

#Gene3D
gene3d_iup_overlap_file = pd.read_csv("/home/guest/Internship/results/Interproscan_Gene3D/05_gene3d_iupred_overlap_above_50.tsv",sep="\t")
gene3d_aiup_overlap_file = pd.read_csv("/home/guest/Internship/results/Interproscan_Gene3D/05_gene3d_aiupred_overlap_above_50.tsv",sep="\t")

#Drop duplicates: filter on the acession numbers and domain names
#Pfam
pfam_filtered_iup_overlap_file=pfam_iup_overlap_file.drop_duplicates(subset=["Accession_number", "Domain"])
pfam_filtered_aiup_overlap_file=pfam_aiup_overlap_file.drop_duplicates(subset=["Accession_number" , "Domain"])

#Gene3D
gene3d_filtered_iup_overlap_file=gene3d_iup_overlap_file.drop_duplicates(subset=["Accession_number", "Domain"])
gene3d_filtered_aiup_overlap_file=gene3d_aiup_overlap_file.drop_duplicates(subset=["Accession_number" , "Domain"])

#Creating new files from the modified versions:
#Pfam
pfam_filtered_iup_overlap_file.to_csv("~/Internship/results/Interproscan_Pfam/08_pfam_iupred_found_domains.tsv", sep="\t", index=False)
pfam_filtered_aiup_overlap_file.to_csv("~/Internship/results/Interproscan_Pfam/08_pfam_aiupred_found_domains.tsv",sep="\t",index=False)
#Gene3D
gene3d_filtered_iup_overlap_file.to_csv("~/Internship/results/Interproscan_Gene3D/08_gene3d_iupred_found_domains.tsv", sep="\t", index=False)
gene3d_filtered_iup_overlap_file.to_csv("~/Internship/results/Interproscan_Gene3D/08_gene3d_aiupred_found_domains.tsv", sep="\t", index=False)

#Counting how many domains are found
#info:domains which overlaps more than 50% with the disordered regions predicted by Iupred and AIUpred

#Counting
pfam_iup_num_domains=len(pfam_filtered_iup_overlap_file)
pfam_aiup_num_domains=len(pfam_filtered_aiup_overlap_file)

gene3d_iup_num_domains=len(gene3d_filtered_iup_overlap_file)
gene3d_aiup_num_domains=len(gene3d_filtered_aiup_overlap_file)

#Create txt files for the result
#Pfam
with open("/home/guest/Internship/results/Interproscan_Pfam/08_pfam_number_of_domains.txt","w") as file_pfam:
    file_pfam.write("Number of found domains annotated by Pfam that overlaps more than 50% with disordered regions predicted by Iupred2a: {0}\n".format(pfam_iup_num_domains))
    file_pfam.write("Number of found domains annotated by Pfam that overlaps more than 50% with disordered regions predicted by Aiupred: {0}\n".format(pfam_aiup_num_domains))

#Gene3D
with open("/home/guest/Internship/results/Interproscan_Gene3D/08_gene3d_number_of_domains.txt","w") as file_gene3d:
    file_gene3d.write("Number of found domains annotated by Gene3D that overlaps more than 50% with disordered regions predicted by Iupred2a: {0}\n".format(gene3d_iup_num_domains))
    file_gene3d.write("Number of found domains annotated by Gene3D that overlaps more than 50% with disordered regions predicted by Aiupred: {0}\n".format(gene3d_aiup_num_domains))