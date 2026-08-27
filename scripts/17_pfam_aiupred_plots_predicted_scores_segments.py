#Creating plots to illustrate the distribution of (normalized) redox-dependent predicted scores of the non-domain regions
#predicted by Aiupred with the original and modified cutoff value of the get_redox_region function

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def txt_to_list(number):
    with open(f"/home/guest/Internship/results/Interproscan_Pfam/17_aiupred_modified_cutoff_{number}_segments_predicted_scores.txt","r") as txt_file:
        data=[float(line.rstrip("\n")) for line in txt_file]
    return data

#Creating list from the predicted scores from the original set up
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_aiupred_segments_predicted_scores.txt", "r") as aiup_segment_file:
    data_cutoff_original = [float(line.rstrip("\n")) for line in aiup_segment_file]

#Creating lists from the predicted scores with the modified cutoff values:
numbers=["04","05","06","07"]
lists=[]

for number in numbers:
    lists.append(txt_to_list(number))

data_cutoff_04=lists[0]
data_cutoff_05=lists[1]
data_cutoff_06=lists[2]
data_cutoff_07=lists[3]  

#Creating violinplot
data=[data_cutoff_original,data_cutoff_04,data_cutoff_05,data_cutoff_06,data_cutoff_07]

plt.figure(figsize=(8, 6))
vp=plt.violinplot(data,showmeans=True)
plt.xticks([1, 2, 3, 4,5],
           ["Original 0.3", "Cutoff 0.4", "Cutoff 0.5", "Cutoff 0.6", "Cutoff 0.7"],rotation=0)
plt.ylabel("Normalized predicted scores")
plt.ylim(0, 1)
plt.yticks(np.arange(0, 1.1, 0.1))
plt.title("Distribution of predicted scores of non-domain regions by Aiupred\nwith different cut-off values")

colors=['blue','red','green','yellow','orange']

for i, body in enumerate(vp["bodies"]):
    body.set_facecolor(colors[i])
    body.set_edgecolor('black')
    body.set_alpha(0.5)
    body.set_linewidth(1.5)

vp['cmins'].set_color('black')
vp['cmaxes'].set_color('black')
vp['cbars'].set_color('black')
vp['cmeans'].set_color('white')
plt.savefig("/home/guest/Internship/results/Interproscan_Pfam/17_violinplot_modified_cutoffs_predicted_scores.png", dpi=300, bbox_inches="tight")
plt.show()
