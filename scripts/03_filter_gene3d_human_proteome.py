# Filter the annotated human proteome, so it only contains the acession number and the
#Gene3D file

import pandas as pd

file_location = "/home/guest/Internship/data/Interproscan_Gene3D/UP000005640_9606.fasta.tsv"

df=pd.read_csv(file_location,sep="\t",header=None)
filtered_df = df[[1,8,9,7]]
filtered_df.columns=["accession","start","end","domain"]
filtered_df.to_csv("/home/guest/Internship/results/Interproscan_Gene3D/03_gene3d_filtered_UP000005640_9606.fasta.tsv",sep="\t",index=False)