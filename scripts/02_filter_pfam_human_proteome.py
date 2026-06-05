# Filter the annotated human proteome, so it only contains the acession number and the
#start and end position of the domains

#Pfam file

import pandas as pd

file_location = "/home/guest/Internship/data/Interproscan_Pfam/UP000005640_9606.fasta.tsv"

df=pd.read_csv(file_location,sep="\t",header=None)
filtered_df = df[[1,8,9,7]]
filtered_df.columns=["accession","start","end","domain"]
filtered_df.to_csv("/home/guest/Internship/results/Interproscan_Pfam/02_pfam_filtered_UP000005640_9606.fasta.tsv",sep="\t",index=False)
