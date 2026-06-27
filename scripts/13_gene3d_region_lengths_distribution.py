#Analyzing the distirbuiton of the lengths of the predicted redox sensitive regions
import pandas as pd
import ast
import numpy as np
import matplotlib.pyplot as plt

#Kérdés: rá kéne e szűrni az ugyanazokhoz a doménekhez tartozó régiókra?
# Sztem nem, mert most az összes régiót nézzük ami átfed valamivel

#Gene3D
#Reading in the files: predicted by Iupred2a and Aiupred overlapping more than 50% and less than 50%
#Iupred2a
gene3d_iup_above = pd.read_csv("/home/guest/Internship/results/Interproscan_Gene3D/05_gene3d_iupred_overlap_above_50.tsv",sep="\t")
gene3d_iup_below = pd.read_csv("/home/guest/Internship/results/Interproscan_Gene3D/05_gene3d_iupred_overlap_below_50.tsv",sep="\t")

#Aiupred
gene3d_aiup_above = pd.read_csv("/home/guest/Internship/results/Interproscan_Gene3D/05_gene3d_aiupred_overlap_above_50.tsv",sep="\t")
gene3d_aiup_below = pd.read_csv("/home/guest/Internship/results/Interproscan_Gene3D/05_gene3d_aiupred_overlap_below_50.tsv",sep="\t")

#Create a file for the results:

##Calculating region lengths:
#Iupred2a

iup_above_region_lengths=[]
for index, iup_above_row in gene3d_iup_above.iterrows():
    iup_above_acc=iup_above_row['Accession_number']
    iup_above_regions= ast.literal_eval(iup_above_row["Region"])
    iup_above_start,iup_above_end=iup_above_regions
    iup_above_region_length=iup_above_end - iup_above_start + 1
    iup_above_region_lengths.append(iup_above_region_length)
    
iup_below_region_lengths=[]
for index, iup_below_row in gene3d_iup_below.iterrows():
    iup_below_acc=iup_below_row['Accession_number']
    iup_below_regions= ast.literal_eval(iup_below_row["Region"])
    iup_below_start,iup_below_end=iup_below_regions
    iup_below_region_length=iup_below_end - iup_below_start + 1
    iup_below_region_lengths.append(iup_below_region_length)
    
#Aiupred
aiup_above_region_lengths=[]
for index, aiup_above_row in gene3d_aiup_above.iterrows():
    aiup_above_acc=aiup_above_row['Accession_number']
    aiup_above_regions= ast.literal_eval(aiup_above_row["Region"])
    aiup_above_start,aiup_above_end=aiup_above_regions
    aiup_above_region_length=aiup_above_end - aiup_above_start + 1
    aiup_above_region_lengths.append(aiup_above_region_length)
    
aiup_below_region_lengths=[]
for index, aiup_below_row in gene3d_aiup_below.iterrows():
    aiup_below_acc=aiup_below_row['Accession_number']
    aiup_below_regions= ast.literal_eval(aiup_below_row["Region"])
    aiup_below_start,aiup_below_end=aiup_below_regions
    aiup_below_region_length=aiup_below_end - aiup_below_start + 1
    aiup_below_region_lengths.append(aiup_below_region_length)  

#Filtering the outliers:IQR method (boxplot rule)
data = [iup_above_region_lengths, iup_below_region_lengths, aiup_above_region_lengths, aiup_below_region_lengths]

filtered_data=[]

for list in data: 
    q1 = np.percentile(list, 10)
    q3 = np.percentile(list, 90)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    filtered_list = [x for x in list if lower_bound <= x <= upper_bound]
    filtered_data.append(filtered_list)

#Creating violin plts for that
plt.figure(figsize=(8, 6))
vp=plt.violinplot(filtered_data,showmeans=True)
plt.xticks([1, 2, 3, 4],
           ["Iupred2a domains", "Iupred2a non-domains", "Aiupred domains 3", "Aiupred non-domains"])
plt.ylabel("Value")
plt.title("Comparison of distributions of region lengths")

colors=['blue','red','blue','red']

for i, body in enumerate(vp["bodies"]):
    body.set_facecolor(colors[i])
    body.set_edgecolor('black')
    body.set_alpha(0.5)
    body.set_linewidth(1.5)

quantilies = [np.percentile(d, [25, 50, 75]) for d in filtered_data]

for i, q in enumerate(quantilies):
    plt.scatter([i+1]*3, q, color="white", s=20,marker="s")

vp['cmins'].set_color('black')
vp['cmaxes'].set_color('black')
vp['cbars'].set_color('black')
vp['cmeans'].set_color('white')
plt.show()

print(vp.keys())
# #Creating histograms for that:

# fig, axes = plt.subplots(2, 2, figsize=(10, 8))

# axes[0, 0].hist(filtered_data[0], bins=10)
# axes[0, 0].set_title("Iupred2a domains")

# axes[0, 1].hist(filtered_data[1], bins=10)
# axes[0, 1].set_title("Iupred2 a non-domains")

# axes[1, 0].hist(filtered_data[2], bins=10)
# axes[1, 0].set_title("Aiupred domains")

# axes[1, 1].hist(filtered_data[3], bins=10)
# axes[1, 1].set_title("Aiupred non-domains")

# plt.tight_layout()
# plt.show()