#Creating a byrplot to visualize the number of found disordered redox-dependent regions 
#predicted by Aiupred with the modified cut off value: 0.4
#Compare it with the number of regions predicted by Iupred2a

import pandas as pd
import matplotlib.pyplot as plt

#Read in the tsv files and count the number of regions in it
aiupred_number_regions = len(pd.read_csv("/home/guest/Internship/results/16_aiupred_modified_cutoff_04_disordered_regions.tsv", sep="\t"))
iupred_number_regions=len(pd.read_csv("/home/guest/Internship/results/01_iupred_disordered_regions.tsv", sep="\t"))

#Creating a barplot:
#How many regions are found by Iupred2a and Aiupred with the modifed cut-off value
regions= [iupred_number_regions,aiupred_number_regions]
tools=['Iupred2a\ncut-off value:0.3','Aiupred_modified\ncut-off value:0.4']

fig,ax=plt.subplots()

bars = ax.bar(tools, regions)
#plt.bar(tools,regions)
ax.bar_label(bars, padding=3)
ax.margins(y=0.2)
plt.title('Number of predicted disordered regions\nIupred2a vs Aiupred with modified cut-off value')
plt.xticks(rotation=0)
plt.ylabel('Number of regions')
plt.savefig("/home/guest/Internship/results/16_barplot_modified_cutoff_04_disordered_regions.png", dpi=300, bbox_inches="tight")
plt.show()