#Creating histograms

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

##Reading in the files, and count the data:
scr01_disord_regions_iupred=pd.read_csv("/home/guest/Internship/results/01_iupred_disordered_regions.tsv",sep="\t")
scr01_disord_regions_aiupred=pd.read_csv("/home/guest/Internship/results/01_aiupred_disordered_regions.tsv",sep="\t")

pfam_annotated_domains=len(pd.read_csv("/home/guest/Internship/data/Interproscan_Pfam/UP000005640_9606.fasta.tsv",sep="\t"))
gene3d_annotated_domains=len(pd.read_csv("/home/guest/Internship/data/Interproscan_Gene3D/UP000005640_9606.fasta.tsv",sep="\t"))

scr04_pfam_iupred_above_50=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_iupred_overlap_above_50.tsv",sep="\t"))
scr04_pfam_iupred_below_50=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_iupred_overlap_below_50.tsv",sep="\t"))
scr04_pfam_aiupred_above_50=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_aiupred_overlap_above_50.tsv",sep="\t"))
scr04_pfam_aiupred_below_50=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_aiupred_overlap_below_50.tsv",sep="\t"))

scr05_gene3d_iupred_above_50=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Gene3D/05_gene3d_iupred_overlap_above_50.tsv",sep="\t"))
scr05_gene3d_iupred_below_50=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Gene3D/05_gene3d_iupred_overlap_below_50.tsv",sep="\t"))
scr05_gene3d_aiupred_above_50=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Gene3D/05_gene3d_aiupred_overlap_above_50.tsv",sep="\t"))
scr05_gene3d_aiupred_below_50=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Gene3D/05_gene3d_aiupred_overlap_below_50.tsv",sep="\t"))

scr08_pfam_domains_iupred=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/08_pfam_iupred_found_domains.tsv",sep="\t"))
scr08_pfam_domains_aiupred=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/08_pfam_aiupred_found_domains.tsv",sep="\t"))
scr08_gene3d_domains_iupred=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Gene3D/08_gene3d_iupred_found_domains.tsv",sep="\t"))
scr08_gene3d_domains_aiupred=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Gene3D/08_gene3d_aiupred_found_domains.tsv",sep="\t"))

scr09_pfam_proteins_iupred=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/09_pfam_iupred_found_proteins.tsv",sep="\t"))
scr09_pfam_proteins_aiupred=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/09_pfam_aiupred_found_proteins.tsv",sep="\t"))
scr09_gene3d_proteins_iupred=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Gene3D/09_gene3d_iupred_found_proteins.tsv",sep="\t"))
scr09_gene3d_proteins_aiupred=len(pd.read_csv("/home/guest/Internship/results/Interproscan_Gene3D/09_gene3d_aiupred_found_proteins.tsv",sep="\t"))

#Processing the disordered regions tsv files:
#counting the lines=found regions:
numb_disord_regions_iupred=len(scr01_disord_regions_iupred)
numb_disord_regions_aiupred=len(scr01_disord_regions_aiupred)

#drop duplicates: filter on the acession numbers:
filtered_proteins_iup=scr01_disord_regions_iupred.drop_duplicates(subset=["accession"])
filtered_proteins_aiup=scr01_disord_regions_aiupred.drop_duplicates(subset=["accession"])

#counting the lines=found proteins:
numb_proteins_iup=len(filtered_proteins_iup)
numb_proteins_aiup=len(filtered_proteins_aiup)

##Creating plots:
#How many regions are found by Iupred2a and Aiupred
regions= [numb_disord_regions_iupred,numb_disord_regions_aiupred]
tools=['Iupred2a','Aiupred']

plt.bar(tools,regions)
plt.title('Number of predicted disordered regions Iupred2a vs Aiupred')
plt.xlabel('Tools')
plt.ylabel('Number of regions')
plt.savefig("/home/guest/Internship/results/11_barplot_number_disordered_regions.png", dpi=300, bbox_inches="tight")
plt.show()

#How many proteins are found by Iupred2a vs Aiupred
proteins= [numb_proteins_iup,numb_proteins_aiup]
tools=['Iupred2a','Aiupred']

plt.bar(tools,proteins)
plt.title('Number of proteins with disordered regions\n predicted by Iupred2a vs Aiupred')
plt.xlabel('Tools')
plt.ylabel('Number of proteins')
plt.savefig("/home/guest/Internship/results/11_barplot_number_proteins_with_disordered_regions.png", dpi=300, bbox_inches="tight")
plt.show()

#Region length distribution of them: already created within the 01_script

# #Pfam vs Gene3D annotated domains/regions?
# annotations= [pfam_annotated_domains,gene3d_annotated_domains]
# databases=['Pfam','Gene3D']

# plt.bar(databases,annotations)
# plt.title('Number of annotated domains/regions by Pfam vs Gene3D')
# plt.xlabel('Databases')
# plt.ylabel('Number of annotated domains/regions')
# plt.show()

#Pfam: How many overlapping regions: more than 50% vs less than 50% Iupred vs Aiupred
df= pd.DataFrame({
    'Overlap above 50%': [scr04_pfam_iupred_above_50,scr04_pfam_aiupred_above_50],
    'Overlap below 50%': [scr04_pfam_iupred_below_50,scr04_pfam_aiupred_below_50]
}, index=['Iupred2a','Aiupred'])

df.plot(kind='bar')
plt.title('Number of overlapping regions by Pfam')
plt.xticks(rotation=0)
plt.ylabel('Number of overlapping regions')
plt.savefig("/home/guest/Internship/results/11_barplot_pfam_number_overlapping_regions.png", dpi=500)
plt.show()

#Gene3D:How many overlapping regions: more than 50% vs less than 50% Iupred vs Aiupred
df= pd.DataFrame({
    'Overlap above 50%': [scr05_gene3d_iupred_above_50,scr05_gene3d_aiupred_above_50],
    'Overlap below 50%': [scr05_gene3d_iupred_below_50,scr05_gene3d_aiupred_below_50]
}, index=['Iupred2a','Aiupred'])

df.plot(kind='bar')
plt.title('Number of overlapping regions by Gene3D')
plt.xticks(rotation=0)
plt.ylabel('Number of overlapping regions')
plt.savefig("/home/guest/Internship/results/11_barplot_gene3d_number_overlapping_regions.png", dpi=500)
plt.show()

#Pfam vs Gene3D: how many overlapping regions more than 50% Aiupred vs Iupred

# data= pd.DataFrame({
#     'Pfam': [scr04_pfam_iupred_above_50,scr04_pfam_aiupred_above_50],
#     'Gene3D': [scr05_gene3d_iupred_above_50,scr05_gene3d_aiupred_above_50]
# }, index=['Iupred2a','Aiupred'])

# data.plot(kind='bar')
# plt.title('Regions that overlaps more than 50%')
# plt.ylabel('Number of overlapping regions above 50%')
# plt.show()

#Pfam vs Gene3D:How many found domains by Aiupred vs Iupred
domains= pd.DataFrame({
    'Iupred2a': [scr08_pfam_domains_iupred,scr08_gene3d_domains_iupred],
    'Aiupred': [scr08_pfam_domains_aiupred,scr08_gene3d_domains_aiupred]
}, index=['Pfam','Gene3D'])

domains.plot(kind='bar')
plt.title('Annotated domains overlapping with disordered regions\n predicted by Iupred2a and Aiupred')
plt.xticks(rotation=0)
plt.ylabel('Number of found domains')
plt.savefig("/home/guest/Internship/results/11_barplot_number_annotated_domains.png", dpi=500)
plt.show()

#Pfam vs Gene3D: How many found proteins by Aiupred vs Iupred
proteins= pd.DataFrame({
    'Iupred2a': [scr09_pfam_proteins_iupred,scr09_gene3d_proteins_iupred],
    'Aiupred': [scr09_pfam_proteins_aiupred,scr09_gene3d_proteins_aiupred]
}, index=['Pfam','Gene3D'])

proteins.plot(kind='bar')
plt.title('Annotated proteins with regions overlapping with disordered regions\n predicted by Iupred2a and Aiupred')
plt.xticks(rotation=0)
plt.ylabel('Number of found proteins')
plt.savefig("/home/guest/Internship/results/11_barplot_number_annotated_proteins.png", dpi=500)
plt.show()