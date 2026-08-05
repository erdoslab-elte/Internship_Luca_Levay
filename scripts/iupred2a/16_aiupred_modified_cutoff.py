#Increasing the cut off that Aiupred uses for predicting disordered regions
#What I modified: only the 0.5 value: increased it (by 20%) to 0.6
#I didn't modify the difference between them also not the values at the expanding part

import pandas as pd
from iupred2a_lib import get_redox_regions
from aiupred import AIUPred
import time

start_time=time.time()

#Creating a function to read a multifasta file:
def multi_fasta_reader(file_location):
    fasta_dat = {}
    header = None
    with open(file_location) as file_handler:
        for line in file_handler:
            if line.startswith('>'):
                header = line.split('|')[1]
                fasta_dat[header] =''
            elif line.strip() and header:
                fasta_dat[header] += line.strip()
    return fasta_dat

#Defining Aiupred functions for prediction of disordered regions
predictor = AIUPred()

def aiupred_score(seq):
    prediction = predictor.predict_disorder(seq)
    return prediction

def aiupred_redox_score(seq):
    return aiupred_score(seq.replace("C", "S"))

#Modification of get_redox_regions function: increase the cut off vlaue from 0.5 to 0.6
def get_redox_regions(redox_values, iupred_values):
    """
    Calculate the redox sensitive regions
    :param redox_values: Redox Y coordinates
    :param iupred_values: IUPred Y coordiantes
    :return:
    """
    patch_loc = {}
    trigger = False
    opening_pos = []
    start, end = 0, 0
    counter = 0
    # Calculate possible position
    for idx, redox_val in enumerate(redox_values):
        if redox_val > 0.5 > iupred_values[idx] and redox_val - iupred_values[idx] > 0.4:
            opening_pos.append(idx)
    # Filter out where not enough possible position is found
    # Enlarge region where enough position is found
    for idx, redox_val in enumerate(redox_values):
        if redox_val - iupred_values[idx] > 0.15 and redox_val >= 0.35:
            if not trigger:
                start = idx
                trigger = True
            if idx in opening_pos:
                counter += 1
            end = idx
        else:
            trigger = False
            if end - start > 14 and counter > 2:
                patch_loc[start] = end
            counter = 0
    if end - start > 14 and counter > 2:
        patch_loc[start] = end
    # Combine close regions
    deletable = []
    for start, end in patch_loc.items():
        for start2, end2 in patch_loc.items():
            if start != start2 and start2 - end < 10 and start2 > start:
                patch_loc[start] = end2
                deletable.append(start2)
    for start in deletable:
        del patch_loc[start]
    return patch_loc

#Reading in the human proteome file:
human_file = multi_fasta_reader('/home/guest/Internship/data/UP000005640_9606.fasta')

#Reading in the original Aiupred tsv file with the predicted disordered regions and creating a list from the accession numbers:
aiup_acc_list=[]

with open("/home/guest/Internship/results/01_aiupred_disordered_regions.tsv") as aiup_file:
    next(aiup_file)
    for line in aiup_file:
        acc=line.rstrip().split("\t")[0]
        if acc not in aiup_acc_list:
            aiup_acc_list.append(acc)

#Creating a new file for the results:
result_file=open("/home/guest/Internship/results/16_aiupred_modified_disordered_regions.tsv","w")
result_file.write("Accession_number\tStart\tEnd\n")

#Adding data to the result file:
#m=0
number_of_regions=0  

for human_acc,seq in human_file.items():
    for aiup_acc in aiup_acc_list:
        if aiup_acc==human_acc:
            aiup_score=aiupred_score(seq)
            redox_score=aiupred_redox_score(seq)
            aiup_redox_region=get_redox_regions(redox_score,aiup_score)
            for start,end in aiup_redox_region.items():
                result_file.write(f"{aiup_acc}\t{start}\t{end}\n")
                #print(aiup_acc,start,end)
                number_of_regions +=1
                #m +=1
        # if m >= 1:
        #     break
result_file.close()

#Creating a txt file for the number of predicted disordered regions by Aiupred:
with open("/home/guest/Internship/results/16_aiupred_modified_number_disordered_regions.txt","w") as file:
    file.write(f"Number of disordered regions predicted by Aiupred with the modified cut-off value: {number_of_regions}")

print(number_of_regions)

end_time=time.time()
elapsed_seconds = end_time - start_time

elapsed_hours = elapsed_seconds / 3600

print(f"Script completed in {elapsed_hours:.2f} hours.")