#Run the Iupred2a and the AIUpred tools on the Human proteome
#ACompare how many disorganised regions are found by the two tools

from iupred2a_lib import iupred, get_redox_regions
import gzip
from aiupred import AIUPred
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

predictor = AIUPred()

def iupred_recdox(seq):
    pred, _ = iupred(seq)

    mutated_seq = seq.replace('C', 'S')

    mutated_prediction, _ = iupred(mutated_seq)
    regions = get_redox_regions(mutated_prediction, pred)
    return regions

def aiupred_recdox(seq):
    pred = predictor.predict_disorder(seq)

    mutated_seq = seq.replace('C', 'S')

    mutated_prediction = predictor.predict_disorder(mutated_seq)
    regions = get_redox_regions(mutated_prediction, pred)
    return regions

def multi_fasta_reader(file_location):
    fasta_dat = {}
    header = None
    with gzip.open(file_location, 'rt') as file_handler:
        for line in file_handler:
            if line.startswith('>'):
                header = line.split('|')[1]
                fasta_dat[header] =''
            elif line.strip() and header:
                fasta_dat[header] += line.strip()
    return fasta_dat


fasta_data = multi_fasta_reader('/home/guest/Internship/data/UP000005640_9606.fasta.gz')

#n = 0

iupred_num_of_reg = 0
aiupred_num_of_reg = 0
iupred_region_lengths = []
aiupred_region_lengths = []
iupred_num_of_prot = 0
aiupred_num_of_prot = 0

results_regions_iupred = []
results_regions_aiupred = []

for acc, seq in fasta_data.items():
    try:
        #print(acc, iupred_recdox(seq), aiupred_recdox(seq))
        iupred_regions = iupred_recdox(seq)
        if len(iupred_regions)!=0:
        # print(len(iupred_recdox(seq)))
            iupred_num_of_reg = iupred_num_of_reg + len(iupred_regions)
            iupred_num_of_prot += 1
            for start,end in iupred_regions.items():
                region_length= end - start + 1
                iupred_region_lengths.append(region_length)
                results_regions_iupred.append({"accession": acc,
                                    "start": start,
                                    "end":end})
        aiupred_regions = aiupred_recdox(seq)
        if len(aiupred_regions)!=0:
     #print(len(aiupred_recdox(seq)))
            aiupred_num_of_reg = aiupred_num_of_reg + len(aiupred_regions)
            aiupred_num_of_prot += 1
            for start,end in aiupred_regions.items():
                region_length= end - start
                aiupred_region_lengths.append(region_length)
                results_regions_aiupred.append({"accession": acc,
                                    "start": start,
                                    "end":end})
    except IndexError:
        continue
    #n += 1
    #if n>20:
        #break

#Creating tsv files that contains the found dirosdered regions: accession number, start and end position of the regions       
df_iupred=pd.DataFrame(results_regions_iupred)
df_iupred.to_csv("/home/guest/Internship/results/01_iupred_disordered_regions.tsv",sep="\t",index=False)

df_aiupred=pd.DataFrame(results_regions_aiupred)
df_aiupred.to_csv("/home/guest/Internship/results/01_aiupred_disordered_regions.tsv",sep="\t",index=False)

#Summarize the results in a txt file
with open("/home/guest/Internship/results/01_number_of_predicted_regions_proteins.txt","w") as file:
    file.write("Number of disordered regions predicted by Iupred2a: {0}\n".format(iupred_num_of_reg))
    file.write("Number of disordered regions predicted by Aiupred: {0}\n".format(aiupred_num_of_reg))
    file.write("Number of proteins containing disordered regions predicted by Iupred2a: {0}\n".format(iupred_num_of_prot))
    file.write("Number of proteins containing disordered regions predicted by Aiupred: {0}\n".format(aiupred_num_of_prot))


#print("Lengths of disordered regions by Iupred2a: {0}".format(iupred_region_lengths))
#print("Lengths of disordered regions by Aiupred: {0}".format(aiupred_region_lengths))

iupred_region_lengths = np.array(iupred_region_lengths)
aiupred_region_lengths = np.array(aiupred_region_lengths)

iupred_bins= np.logspace(np.log10(iupred_region_lengths.min()),np.log10(iupred_region_lengths.max()),60)
aiupred_bins= np.logspace(np.log10(aiupred_region_lengths.min()),np.log10(aiupred_region_lengths.max()),60)

#Creating histograms to present the distribution of the disogranized regions'lengths

plt.figure()
plt.hist(iupred_region_lengths, bins=iupred_bins)
plt.xscale("log")
plt.xlabel("Region lengths")
plt.ylabel("Frequency")
plt.title("Region lenghts by Iupred2a")
plt.savefig("home/guest/Internship/results/01_histogram_regions_length_iupred.png",dpi=300)

plt.figure()
plt.hist(aiupred_region_lengths, bins=aiupred_bins)
plt.xscale("log")
plt.xlabel("Region lengths")
plt.ylabel("Frequency")
plt.title("Region lenghts by Aiupred")
plt.savefig("/home/guest/Internship/results/01_histogram_regions_length_aiupred.png",dpi=300)

plt.show()