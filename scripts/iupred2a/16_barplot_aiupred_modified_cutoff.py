#Creating a byrplot to visualize the number of found disordered redox-dependent regions 
#predicted by Aiupred with the modified cut off value
#Compare it with the number of regions predicted by Iupred2a

import pandas as pd
import matplotlib.pyplot as plt

#Read in the tsv files and count the number of regions in it
aiupred_number_regions = len(pd.read_csv("/home/guest/Internship/results/16_aiupred_modified_disordered_regions.tsv", sep="\t"))
iupred_number_regions=len(pd.read_csv("/home/guest/Internship/results/01_iupred_disordered_regions.tsv", sep="\t"))

#Creating a barplot:
#How many regions are found by Iupred2a and Aiupred with the modifed cut-off value
regions= [iupred_number_regions,aiupred_number_regions]
tools=['Iupred2a','Aiupred_modified']

plt.bar(tools,regions)
plt.title('Number of predicted disordered regions Iupred2a vs Aiupred with modified cut-off value')
plt.xlabel('Tools')
plt.ylabel('Number of regions')
plt.savefig("/home/guest/Internship/results/16_barplot_modified_number_disordered_regions.png", dpi=300, bbox_inches="tight")
plt.show()