#Calculation of redox-state dependent disorder prediction scores for the domain and non-domain segments identified by Iupred2a
#Comparing the distribution of the prediction scores between domain and non-domain segments

import pandas as pd
from iupred2a_lib import iupred,iupred_redox,get_redox_regions
import time
from ast import literal_eval

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

#Read in the Iupred2a predicted segment files:
##Domains: segments that overlap more than 50% with annotated domains:
iup_domain_file=pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_iupred_overlap_above_50.tsv",sep="\t")

##Non-domains/Segments: segments that overlap less than 50% with annotated domains:
iup_segment_file=pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_iupred_overlap_below_50.tsv",sep="\t")

#Create a list of the accesion numbers of the predicted redox-dependent disordered regions by Iupred2a:
iup_domain_list = create_acc_list_dict(iup_domain_file)[0]
iup_segment_list=create_acc_list_dict(iup_segment_file)[0]

#Create a dictionary of the accesion numbers and the predicted redox-dependent disordered regions by Iupred2a:
iup_above_dict = create_acc_list_dict(iup_domain_file)[1]
iup_below_dict=create_acc_list_dict(iup_segment_file)[1]

# #Save the predicted redox sensitive segments' sequences in a dictionary
iup_domain_found_ids=create_id_dict(fasta_data,iup_domain_list)
iup_segment_found_ids=create_id_dict(fasta_data,iup_segment_list)

# #Calculating the iupred scores for the sequences: calculating the difference of the two scores for a residual,then adding them together
##Calculating for DOMAINS predicted by IUPRED2A:

list_domain_score_diff=[]
dict_domain_score_diff={}

#n=0
for id,iup_seq in iup_domain_found_ids.items():
        for acc,regions in iup_above_dict.items():
            if acc==id:
                iup_score=iupred(iup_seq)[0]
                iup_redox_score=iupred_redox(iup_seq)[0]
                redox_regions=get_redox_regions(iup_redox_score, iup_score)
                redox_regions_list=[[key,value] for key,value in redox_regions.items()]
                #print(id,redox_regions_list)
                #m=0
                for redox_region in redox_regions_list:
                    for region in regions:
                        if region == redox_region:
                            for start,end in regions:
                                iup_score_region=iup_score[start:end]
                                iup_redox_score_region=iup_redox_score[start:end]
                                norm_score_diff= sum(abs(iup_redox_score_region[i]-iup_score_region[i]) for i, _ in enumerate(iup_redox_score_region))/(end-start)
                                list_domain_score_diff.append(norm_score_diff)
                                if id not in dict_domain_score_diff:
                                    dict_domain_score_diff[id]={}
                                #if (start,end) not in dict_domain_score_diff[id]:
                                dict_domain_score_diff[id][(start, end)] = norm_score_diff
                                #print("result:" + id, start, end, norm_score_diff)
                                #m +=1
        #                         if m> 4:
        #                             break
        # n += 1
        # if n >3:
        #     break

# print(list_domain_score_diff)
# print(dict_domain_score_diff)

##Calculating for SEGMENTS/NON-DOMAINS predicted by IUPRED2A:

list_segment_score_diff=[]
dict_segment_score_diff={}

#n=0
for id,iup_seq in iup_segment_found_ids.items():
        for acc,regions in iup_below_dict.items():
            if acc==id:
                iup_score=iupred(iup_seq)[0]
                iup_redox_score=iupred_redox(iup_seq)[0]
                redox_regions=get_redox_regions(iup_redox_score, iup_score)
                redox_regions_list=[[key,value] for key,value in redox_regions.items()]
                #print(id,redox_regions_list)
                #m=0
                for redox_region in redox_regions_list:
                    for region in regions:
                        if region == redox_region:
                            for start,end in regions:
                                iup_score_region=iup_score[start:end]
                                iup_redox_score_region=iup_redox_score[start:end]
                                norm_score_diff= sum(abs(iup_redox_score_region[i]-iup_score_region[i]) for i, _ in enumerate(iup_redox_score_region))/(end-start)
                                list_segment_score_diff.append(norm_score_diff)
                                if id not in dict_segment_score_diff:
                                    dict_segment_score_diff[id]={}
                                dict_segment_score_diff[id][(start, end)] = norm_score_diff
                                #print("result:" + id, start, end, norm_score_diff)
        #                         m +=1
        #                         if m> 4:
        #                             break
        # n += 1
        # if n >2:
        #     break

# print(list_segment_score_diff)
# print(dict_segment_score_diff)

# #Create text files to store these lists:
#For the domains:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_iupred_domains_predicted_scores.txt", "w") as file:
    for numbers in list_domain_score_diff:
        file.write(f"{numbers}\n")

#For the segments:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_iupred_segments_predicted_scores.txt", "w") as file:
    for numbers in list_segment_score_diff:
        file.write(f"{numbers}\n")

#Create tsv files for the accession numbers and the sum of the confidence scores belonging to them
#For the domains:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_iupred_domains_predicted_scores.tsv","w") as tsv:
    tsv.write("Accession_number\tStart\tEnd\tPredicted_score\n")
    for acc, regions in dict_domain_score_diff.items():
        for (start, end), score in regions.items():
            tsv.write(f"{acc}\t{start}\t{end}\t{score}\n")

#For the non-domains:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_iupred_segments_predicted_scores.tsv","w") as tsv:
    tsv.write("Accession_number\tStart\tEnd\tPredicted_score\n")
    for acc, regions in dict_segment_score_diff.items():
        for (start, end), score in regions.items():
            tsv.write(f"{acc}\t{start}\t{end}\t{score}\n")

end_time=time.time()
elapsed_seconds = end_time - start_time
elapsed_hours = elapsed_seconds / 3600

print(f"Script completed in {elapsed_hours:.2f} hours.")