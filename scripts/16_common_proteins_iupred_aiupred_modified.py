#Collect those proteins that are found both by Iupred2a and Aiupred with the modified cut-off value in Aiupred
#Collect those proteins that are found both by Iupred2a and Aiupred
import pandas as pd
import matplotlib.pyplot as plt

#Read in the files
iupred_file=pd.read_csv("/home/guest/Internship/results/01_iupred_disordered_regions.tsv",sep="\t")
modified_aiupred_file=pd.read_csv("/home/guest/Internship/results/16_aiupred_modified_disordered_regions.tsv",sep="\t")
original_common_protein_file=pd.read_csv("/home/guest/Internship/results/10_common_found_proteins.tsv",sep="\t")
                                         
#Creating a list from disordered region predicted by Iupred2a
iupred_acc_list=[]

for _,iup_row in iupred_file.iterrows():
    iup_acc=iup_row["accession"]
    iupred_acc_list.append(iup_acc)

#Creating a list from disordered region predicted by Aiupred with the modified cut-off value
modified_aiupred_acc_list=[]

for _,aiup_row in modified_aiupred_file.iterrows():
    aiup_acc=aiup_row["Accession_number"]
    modified_aiupred_acc_list.append(aiup_acc)


#Creating three lists:accession numbers of proteins with disordered regions: one with common proteins found by both of the tools, one with found only by Iupred2a, one found only by Aiupred
common_acc_list = []
only_iupred_list = []
only_aiupred_list = []

for iup_acc in iupred_acc_list:
    if iup_acc in modified_aiupred_acc_list:
        if iup_acc not in common_acc_list:
            common_acc_list.append(iup_acc)
    else:
        if iup_acc not in only_iupred_list:
            only_iupred_list.append(iup_acc)

for aiup_acc in modified_aiupred_acc_list:
    if aiup_acc not in iupred_acc_list:
        if aiup_acc not in only_aiupred_list:
            only_aiupred_list.append(aiup_acc)

#Counting the number of found proteins with the modified Aiupred setting:
c = len(common_acc_list)
oi = len(only_iupred_list)
oai = len(only_aiupred_list)

#Counting the number of found proteins with the original Aiupred setting:
found_proteins=original_common_protein_file.count()

original_c=(found_proteins["Common proteins"])
original_oi=(found_proteins["Only by Iupred2a"])
original_oai=(found_proteins["Only by Aiupred"])

#Creating a txt file with the number of found proteins in the 3 categories
with open("/home/guest/Internship/results/16_modified_number_common_proteins.txt","w") as file:
    file.write(f"Three categories of found proteins by Iupred2a and Aiupred with the modified cut-off value\n")
    file.write(f"Number of proteins found only by Iupred2a: {oi}\n")
    file.write(f"Number of proteins found only by modified Aiupred: {oai}\n")
    file.write(f"Number of common proteins found both with Iupred2a and modified Aiupred: {c}\n")


#Creating a tsv file for the 3 categories of the found proteins:
df = pd.DataFrame({
    "Only by Iupred2a": pd.Series(only_iupred_list),
    "Common proteins": pd.Series(common_acc_list),
    "Only by Aiupred": pd.Series(only_aiupred_list),
})

df.to_csv("/home/guest/Internship/results/16_modified_common_found_proteins.tsv", sep="\t", index=False)

#Creating a summary barplot:
proteins= [oi,c,oai]
categories=['Only by Iupred2a','Common proteins','Only by Aiupred']

fig,ax=plt.subplots()

bars = ax.bar(categories, proteins,color=['blue', 'orange', 'green'])
ax.bar_label(bars, padding=3)
ax.margins(y=0.2)
plt.title('Three categories of proteins with disordered regions\n found by Iupred2a vs modified Aiupred')
plt.xticks(rotation=0)
plt.ylabel('Number of proteins')
plt.savefig("/home/guest/Internship/results/16_barplot_modified_common_found_proteins.png", dpi=300, bbox_inches="tight")
plt.show()

#Creating a comparison barplot of the 3 categories of found proteins with the original Aiupred and the modified Aiupred tool
df= pd.DataFrame({
    'Only by Iupred2a': [original_oi,oi],
    'Common proteins': [original_c,c],
    'Only by Aiupred': [original_oai,oai]
}, index=['Original Aiupred','Modified Aiupred'])

fig,ax=plt.subplots()
ax.margins(y=0.2)
df.plot(kind='bar',ax=ax, color=['blue','orange','green'])

for container in ax.containers:
    ax.bar_label(container)
    
plt.title('Comparison of found proteins\nby Iupred2a and original Aiupred vs modified Aiupred settings')
plt.xticks(rotation=0)
plt.ylabel('Number of proteins')
plt.savefig("/home/guest/Internship/results/16_barplot_comparison_common_found_proteins.png", dpi=500)
plt.show()

