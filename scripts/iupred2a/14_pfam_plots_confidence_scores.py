#Creating violinplots to illustrate the distribution of redox potential confidence scores of the domain and non-domain regions
#predicted by Iupred2a and Aiupred

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#Load in the txt files:
#Iupred2a domains:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_iupred_domains_confidence_scores.txt", "r") as iup_domain_file:
    iup_domain_data = [float(line.rstrip("\n")) for line in iup_domain_file]
#print(iup_domain_data)

#Iupred2a non-domains/regions:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_iupred_regions_confidence_scores.txt", "r") as iup_region_file:
    iup_region_data = [float(line.rstrip("\n")) for line in iup_region_file]
#print(iup_region_data)

#Aiupred domains:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_aiupred_domains_confidence_scores.txt", "r") as aiup_domain_file:
    aiup_domain_data = [float(line.rstrip("\n")) for line in aiup_domain_file]
# #print(aiup_domain_data)

#Aiupred non-domains/regions:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_aiupred_regions_confidence_scores.txt", "r") as aiup_region_file:
    aiup_region_data = [float(line.rstrip("\n")) for line in aiup_region_file]
# print(aiup_region_data)

#CFiltering the data:
data=[iup_domain_data,iup_region_data,aiup_domain_data,aiup_region_data]

filtered_data=[]

for list in data: 
    q1 = np.percentile(list, 30)
    q3 = np.percentile(list, 70)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    filtered_list = [x for x in list if lower_bound <= x <= upper_bound]
    filtered_data.append(filtered_list)

#Creating violinplot
# plt.figure(figsize=(8, 6))
# vp=plt.violinplot(filtered_data,showmeans=True)
# plt.xticks([1, 2, 3, 4],
#            ["Iupred2a domains", "Iupred2a non-domains", "Aiupred domains", "Aiupred non-domains"])
# plt.ylabel("Confidence scores")
# plt.title("Distribution of predicted scores for disordered regions")

# colors=['blue','red','blue','red']

# for i, body in enumerate(vp["bodies"]):
#     body.set_facecolor(colors[i])
#     body.set_edgecolor('black')
#     body.set_alpha(0.5)
#     body.set_linewidth(1.5)

# vp['cmins'].set_color('black')
# vp['cmaxes'].set_color('black')
# vp['cbars'].set_color('black')
# vp['cmeans'].set_color('white')
# plt.savefig("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_violinplot_confidence_scores.png", dpi=300, bbox_inches="tight")
# plt.show()

#Creating boxplot
plt.figure(figsize=(8, 6))
bp=plt.boxplot(filtered_data,showmeans=True,patch_artist= True,notch='True')
plt.xticks([1, 2, 3, 4],
           ["Iupred2a domains", "Iupred2a non-domains", "Aiupred domains", "Aiupred non-domains"])
plt.ylabel("Confidence scores")
plt.title("Distribution of predicted scores for disordered regions")

colors=['blue','red','blue','red']

for i, box in enumerate(bp["boxes"]):
    box.set_facecolor(colors[i])
    box.set_alpha(0.5)
    box.set_linewidth(1.5)

for median in bp['medians']:
    median.set(color ='yellow',
               linewidth = 2)
    
# bp['cmins'].set_color('black')
# bp['cmaxes'].set_color('black')
# bp['cbars'].set_color('black')
#bp['means'].set_color('white')
plt.savefig("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_boxplot_confidence_scores.png", dpi=300, bbox_inches="tight")
plt.show()

# dict_keys([
#     'whiskers',
#     'caps',
#     'boxes',
#     'medians',
#     'fliers',
#     'means'
# ])