#Collect those proteins that are found both by Iupred2a and Aiupred
import pandas as pd
import matplotlib.pyplot as plt

#Read in the files
iupred_file=pd.read_csv("/home/guest/Internship/results/01_iupred_disordered_regions.tsv",sep="\t")
aiupred_file=pd.read_csv("/home/guest/Internship/results/01_aiupred_disordered_regions.tsv",sep="\t")

#Creating a list from disordered region predicted by Iupred2a
iupred_acc_list=[]

for _,iup_row in iupred_file.iterrows():
    iup_acc=iup_row["accession"]
    iupred_acc_list.append(iup_acc)

#Creating a list from disordered region predicted by Aiupred
aiupred_acc_list=[]

for _,aiup_row in aiupred_file.iterrows():
    aiup_acc=aiup_row["accession"]
    aiupred_acc_list.append(aiup_acc)

#Creating three lists:accession numbers of proteins with disordered regions: one with common proteins found by both of the tools, one with found only by Iupred2a, one found only by Aiupred
common_acc_list = []
only_iupred_list = []
only_aiupred_list = []

for iup_acc in iupred_acc_list:
    if iup_acc in aiupred_acc_list:
        if iup_acc not in common_acc_list:
            common_acc_list.append(iup_acc)
    else:
        if iup_acc not in only_iupred_list:
            only_iupred_list.append(iup_acc)

for aiup_acc in aiupred_acc_list:
    if aiup_acc not in iupred_acc_list:
        if aiup_acc not in only_aiupred_list:
            only_aiupred_list.append(aiup_acc)

#Counting the number of found proteins
c = len(common_acc_list)
oi = len(only_iupred_list)
oai = len(only_aiupred_list)

#Creating a txt file with the number of found proteins in the 3 categories
with open("/home/guest/Internship/results/10_number_common_proteins.txt","w") as file:
    file.write(f"Number of proteins found only by Iupred2a: {oi}\n")
    file.write(f"Number of proteins found only by Aiupred: {oai}\n")
    file.write(f"Number of common proteins found both with Iupred2a and Aiupred: {c}\n")


#Creating a tsv file for the 3 categories of the found proteins:
df = pd.DataFrame({
    "Common proteins": pd.Series(common_acc_list),
    "Only by Iupred2a": pd.Series(only_iupred_list),
    "Only by Aiupred": pd.Series(only_aiupred_list),
})

df.to_csv("/home/guest/Internship/results/10_common_found_proteins.tsv", sep="\t", index=False)

proteins= [c,oi,oai]
categories=['Common proteins','Only by Iupred2a', 'Only by Aiupred']

#Creating a summary barplot:
plt.bar(categories,proteins)
plt.title('Three categories of proteins with disordered regions\n found by Iupred2a vs Aiupred')
plt.xticks(rotation=0)
plt.ylabel('Number of proteins')
plt.savefig("/home/guest/Internship/results/10_barplot_number_common_found_proteins.png", dpi=300, bbox_inches="tight")
plt.show()
