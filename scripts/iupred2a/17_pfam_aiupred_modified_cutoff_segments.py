#Calculation of normalized redox-state dependent disorder prediction scores
# for the non-domain regions identified by AIUPred using modified cutoff values

from aiupred import AIUPred
import pandas as pd
from ast import literal_eval
import time

start_time=time.time()

#Creating functions to calculate the redox-dependent disorder prediction scores:
predictor = AIUPred()

def aiupred_score(seq):
    prediction = predictor.predict_disorder(seq)
    return prediction

def aiupred_redox_score(seq):
    return aiupred_score(seq.replace("C", "S"))

#Modifying the get_redox_regions function: adding the cutoff parameter
#so the cutoff can be changed as part of the parameters of the function
#cutoff the difference between the calcualted aupred scores and redox scores
#default cutoff value is 0.3

def modified_get_redox_regions(redox_values, iupred_values, cutoff):
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
        if redox_val > 0.5 > iupred_values[idx] and redox_val - iupred_values[idx] > cutoff:
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

#Function for creating a list of accession numbers:
def create_acc_list_dict(file):
    list=[]
    dict={}
    for _,row in file.iterrows():
        acc=row["Accession_number"]
        region=literal_eval(row["Region"])
        if acc not in list:
            list.append(acc)
            dict[acc]=[]
        dict[acc].append(region)
    return list, dict

#Function for creating a dictionary based on same accession numbers:
#filter it so it only contains sequences that are at maximum 1000 aminoacids long

def create_id_dict(fasta_file,acc_list):
    found_ids_dict={}
    for id,seq, in fasta_file.items():
        for acc in acc_list:
            if id == acc and acc not in found_ids_dict:
                if len(seq) <=1000:
                    found_ids_dict[acc]=""
                    found_ids_dict[acc]=seq
                else:continue
    return found_ids_dict

#Read in the multifasta file (with the defined multifasta function)
fasta_data = multi_fasta_reader('/home/guest/Internship/data/UP000005640_9606.fasta')

#Read in the Aiupred predicted region file:
##Non-domains/Segments:called here segments: regions that overlap less than 50% with annotated domains:
aiup_segment_file=pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_aiupred_overlap_below_50.tsv",sep="\t")

#Create a list of the accesion numbers of the predicted redox-state dependent disordered regions by Aiupred:
#For the non-domain regions:
aiup_segment_list=create_acc_list_dict(aiup_segment_file)[0]

#Create a dictionary of the proteins' accesion numbers and the predicted redox-dependent disordered regions by Aiupred:
#For the non-domain regions:
aiup_below_dict=create_acc_list_dict(aiup_segment_file)[1]

#Save the proteins and their sequences in a dictionary
aiup_segment_found_ids=create_id_dict(fasta_data,aiup_segment_list)

#Creating dictionaries for the calculated aiupred scores and redox aiupred socres for the selected proteins
aiup_scores = {}
aiup_redox_scores = {}

for id, seq in aiup_segment_found_ids.items():
    aiup_scores[id] = aiupred_score(seq)
    aiup_redox_scores[id] = aiupred_redox_score(seq)


#Calculating the aiupred scores for the sequences: calculating the difference of the two scores for a residual,then adding them together
##Calculating for SEGMENTS predicted by AIUPRED:

def predicted_score_calculation(aiup_scores,aiup_redox_scores,aiup_below_dict,cutoff):

    list_score_diff=[]
    dict_score_diff={}

    for id,scores in aiup_scores.items():
        redox_scores= aiup_redox_scores[id]
        redox_regions=modified_get_redox_regions(redox_scores,scores,cutoff)
        redox_regions_list=[[red_start,red_end] for red_start,red_end in redox_regions.items()]
        for redox_region in redox_regions_list:
            for acc,below_regions in aiup_below_dict.items():
                for region in below_regions:
                    if acc== id and region == redox_region:
                        for start,end in below_regions:
                            aiup_score_region=scores[start:end]
                            aiup_redox_score_region=redox_scores[start:end]
                            norm_score_diff= sum(abs(aiup_redox_score_region[i]-aiup_score_region[i]) for i, _ in enumerate(aiup_redox_score_region))/(end-start)
                            list_score_diff.append(norm_score_diff)
                            if id not in dict_score_diff:
                                dict_score_diff[id]={}
                            dict_score_diff[id][(start,end)]= norm_score_diff
    return list_score_diff,dict_score_diff

modified_cutoff_04=predicted_score_calculation(aiup_scores,aiup_redox_scores,aiup_below_dict,0.4)
modified_cutoff_05=predicted_score_calculation(aiup_scores,aiup_redox_scores,aiup_below_dict,0.5)
modified_cutoff_06=predicted_score_calculation(aiup_scores,aiup_redox_scores,aiup_below_dict,0.6)
modified_cutoff_07=predicted_score_calculation(aiup_scores,aiup_redox_scores,aiup_below_dict,0.7)

#Create text files to store these lists:
def create_txt(modified_cutoff,cutoff_number):
    with open(f"/home/guest/Internship/results/Interproscan_Pfam/17_aiupred_modified_cutoff_{cutoff_number}_segments_predicted_scores.txt","w") as file:
        for numbers in modified_cutoff[0]:
            file.write(f"{numbers}\n")

#Create tsv files to store these results:
def create_tsv(modified_cutoff,cutoff_number):
    with open(f"/home/guest/Internship/results/Interproscan_Pfam/17_aiupred_modified_cutoff_{cutoff_number}_segments_predicted_scores.tsv","w") as tsv:
        tsv.write("Accession_number\tStart\tEnd\tPredicted_score\n")
        for acc,regions in modified_cutoff[1].items():
            for (start,end), score in regions.items():
                tsv.write(f"{acc}\t{start}\t{end}\t{score}\n")

create_txt(modified_cutoff_04,"04")
create_txt(modified_cutoff_05,"05")
create_txt(modified_cutoff_06,"06")
create_txt(modified_cutoff_07,"07")

create_tsv(modified_cutoff_04,"04")
create_tsv(modified_cutoff_05,"05")
create_tsv(modified_cutoff_06,"06")
create_tsv(modified_cutoff_07,"07")


end_time=time.time()
elapsed_seconds = end_time - start_time

elapsed_hours = elapsed_seconds / 3600

print(f"Script completed in {elapsed_hours:.2f} hours.")