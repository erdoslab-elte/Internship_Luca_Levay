#Creating plots to illustrate the distribution of (normalized) redox-dependent predicted scores of the domain and non-domain regions
#predicted by Iupred2a and Aiupred

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#Load in the txt files:
#Iupred2a domains:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_iupred_domains_predicted_scores.txt", "r") as iup_domain_file:
    iup_domain_data = [float(line.rstrip("\n")) for line in iup_domain_file]
#print(iup_domain_data)

#Iupred2a non-domains/segments:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_iupred_segments_predicted_scores.txt", "r") as iup_segment_file:
    iup_segment_data = [float(line.rstrip("\n")) for line in iup_segment_file]
#print(iup_segment_data)

#Aiupred domains:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_aiupred_domains_predicted_scores.txt", "r") as aiup_domain_file:
    aiup_domain_data = [float(line.rstrip("\n")) for line in aiup_domain_file]
# #print(aiup_domain_data)

#Aiupred non-domains/segments:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_aiupred_segments_predicted_scores.txt", "r") as aiup_segment_file:
    aiup_segment_data = [float(line.rstrip("\n")) for line in aiup_segment_file]
# print(aiup_segment_data)

#Filtering the data:
data=[iup_domain_data,iup_segment_data,aiup_domain_data,aiup_segment_data]

# filtered_data=[]

# for list in data: 
#     q1 = np.percentile(list, 30)
#     q3 = np.percentile(list, 70)
#     iqr = q3 - q1
#     lower_bound = q1 - 1.5 * iqr
#     upper_bound = q3 + 1.5 * iqr
#     filtered_list = [x for x in list if lower_bound <= x <= upper_bound]
#     filtered_data.append(filtered_list)

#Creating violinplot
plt.figure(figsize=(8, 6))
vp=plt.violinplot(data,showmeans=True)
plt.xticks([1, 2, 3, 4],
           ["Iupred2a domains", "Iupred2a non-domains", "Aiupred domains", "Aiupred non-domains"])
plt.ylabel("Normalized predicted scores")
plt.title("Distribution of predicted scores for disordered regions")

colors=['blue','red','blue','red']

for i, body in enumerate(vp["bodies"]):
    body.set_facecolor(colors[i])
    body.set_edgecolor('black')
    body.set_alpha(0.5)
    body.set_linewidth(1.5)

vp['cmins'].set_color('black')
vp['cmaxes'].set_color('black')
vp['cbars'].set_color('black')
vp['cmeans'].set_color('white')
plt.savefig("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_violinplot_predicted_scores.png", dpi=300, bbox_inches="tight")
plt.show()

#Creating boxplot
plt.figure(figsize=(8, 6))
bp=plt.boxplot(data,showmeans=True,patch_artist= True,notch='True')
plt.xticks([1, 2, 3, 4],
           ["Iupred2a domains", "Iupred2a non-domains", "Aiupred domains", "Aiupred non-domains"])
plt.ylabel("Normalized predicted scores")
plt.title("Distribution of predicted scores for disordered regions")

colors=['blue','red','blue','red']

for i, box in enumerate(bp["boxes"]):
    box.set_facecolor(colors[i])
    box.set_alpha(0.5)
    box.set_linewidth(1.5)

for median in bp['medians']:
    median.set(color ='yellow',
               linewidth = 2)
    
plt.savefig("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_boxplot_predicted_scores.png", dpi=300, bbox_inches="tight")
plt.show()

# dict_keys([
#     'whiskers',
#     'caps',
#     'boxes',
#     'medians',
#     'fliers',
#     'means'
# ])