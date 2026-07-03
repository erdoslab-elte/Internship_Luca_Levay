#Calculation of redox-state dependent disorder prediction scores for the domain and non-domain regions identified by Iupred2a
#Comparing the distribution of the prediction scores between domain and non-domain regions

import pandas as pd
from iupred2a_lib import iupred,iupred_redox
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

#Read in the multifasta file (with the defined multifasta function)
fasta_data = multi_fasta_reader('/home/guest/Internship/data/UP000005640_9606.fasta')

#Read in the Iupred2a predicted region files:
##Domains: regions that overlap more than 50% with annotated domains:
iup_domain_file=pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_iupred_overlap_above_50.tsv",sep="\t")

##Non-domains/Regions: regions that overlap less than 50% with annotated domains:
iup_region_file=pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_iupred_overlap_below_50.tsv",sep="\t")

#Create a list of the accesion numbers of the predicted redox-state dependent disordered regions by Iupred2a:
iup_domain_data = []
for _,iup_row in iup_domain_file.iterrows():
    iup_acc=iup_row["Accession_number"]
    if iup_acc not in iup_domain_data:
        iup_domain_data.append(iup_acc)
# print(f"iup domain data: {iup_domain_data[:3]}")

iup_region_data = []
for _,iup_row in iup_region_file.iterrows():
    iup_acc=iup_row["Accession_number"]
    if iup_acc not in iup_region_data:
        iup_region_data.append(iup_acc)
# print(f"iup region data: {iup_region_data[:3]}")

#Save the predicted redox sensitive regions' sequences in a dictionary
iup_domain_found_ids={}
iup_region_found_ids={}

#n=0

for id,seq in fasta_data.items():
    for iup_dom_acc in iup_domain_data:
        if id == iup_dom_acc and iup_dom_acc not in iup_domain_found_ids:
            if len(seq) <= 1000:
                iup_domain_found_ids[iup_dom_acc]=""
                iup_domain_found_ids[iup_dom_acc]=seq
            else: continue
    for iup_reg_acc in iup_region_data:
        if id == iup_reg_acc and iup_reg_acc not in iup_region_found_ids:
            if len(seq) <= 1000:
                iup_region_found_ids[iup_reg_acc]=""
                iup_region_found_ids[iup_reg_acc]=seq
            else: continue     
    # n += 1
    # if n >20:
    #     break
#print(iup_domain_found_ids)
#print(iup_region_found_ids)

#Calculating the iupred scores for the sequences: calculating the difference of the two scores for a residual,then adding them together
##Calculating for DOMAINS predicted by IUPRED2A:
# m=0
# p=0

iup_domain_score=[]
iup_domain_redox_score=[]
list_sum_domain_diff=[]
dict_sum_domain_diff={}

for iup_acc,iup_seq in iup_domain_found_ids.items():
    #print(f"iup acc: {iup_acc}")
    #print(iup_seq)
    iup_score=iupred(iup_seq)[0]
    iup_redox_score=iupred_redox(iup_seq)[0]
    #print(f"iup score: {iup_score[0:3]}")
    #print(f"iup redox score: {iup_redox_score[0:3]}")
    list_diff=[] #creating a list for the calculated differences
    dict_sum_domain_diff[iup_acc]=0
    redox_score_list=[]
    for score in iup_score:
        #print(f"score: {score}")
        for redox_score in iup_redox_score:
            if redox_score in redox_score_list: #so that it just compares it with one redox score not all of them
                continue
            else:
                redox_score_list.append(redox_score)
                #print(f"redox score: {redox_score}")
                score_diff=abs(score-redox_score)
                #print(f"score difference :{score_diff}")
                list_diff.append(score_diff)
                #m += 1
                break            
        # if m >= 3:
        #     break
    #print(f"list of differences: {list_diff}")
    sum_diff=sum(list_diff)
    #print(f"sum of differences: {sum_diff}")
    dict_sum_domain_diff[iup_acc]=sum_diff
    list_sum_domain_diff.append(sum_diff)
    # p +=1
    # if p >=3:
    #     break
    # if m>=3:
    #     break
# print(f"sum of differences in domain sequences: {list_sum_domain_diff}")
# print(f"dictionary of domains and score differeneces: {dict_sum_domain_diff}")

##Calculating for REGIONS/NON-DOMAINS predicted by IUPRED2A:
# m=0
# p=0

iup_region_score=[]
iup_region_redox_score=[]
list_sum_region_diff=[]
dict_sum_region_diff={}

for iup_acc,iup_seq in iup_region_found_ids.items():
    #print(f"iup acc: {iup_acc}")
    #print(iup_seq)
    iup_score=iupred(iup_seq)[0]
    iup_redox_score=iupred_redox(iup_seq)[0]
    #print(f"iup score: {iup_score[0:3]}")
    #print(f"iup redox score: {iup_redox_score[0:3]}")
    list_diff=[] #creating a list for the calculated differences
    dict_sum_region_diff[iup_acc]=0
    redox_score_list=[]
    for score in iup_score:
        #print(f"score: {score}")
        for redox_score in iup_redox_score:
            if redox_score in redox_score_list: #so that it just compares it with one redox score not all of them
                continue
            else:
                redox_score_list.append(redox_score)
                #print(f"redox score: {redox_score}")
                score_diff=abs(score-redox_score)
                #print(f"score difference :{score_diff}")
                list_diff.append(score_diff)
                # m += 1
                break            
        # if m >= 3:
        #     break
    #print(f"list of differences: {list_diff}")
    sum_diff=sum(list_diff)
    #print(f"sum of differences: {sum_diff}")
    dict_sum_region_diff[iup_acc]=sum_diff
    list_sum_region_diff.append(sum_diff)
    # p +=1
    # if p >=2:
    #     break
    # if m>=3:
    #     break
#print(f"sum of differences in region sequences: {list_sum_region_diff}")
#print(f"dictionary of domains and score differeneces: {dict_sum_domain_diff}")

#Creat text files to store these lists:
#For the domains:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_iupred_domains_confidence_scores.txt", "w") as file:
    for numbers in list_sum_domain_diff:
        file.write(f"{numbers}\n")

#For the regions:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_iupred_regions_confidence_scores.txt", "w") as file:
    for numbers in list_sum_region_diff:
        file.write(f"{numbers}\n")

#Create tsv files for the accession numbers and the sum of the confidence scores belonging to them
#For the domains:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_iupred_domains_confidence_scores.tsv","w") as tsv:
    tsv.write("Accession_number\tSum_of_confidence_scores\n")
    for acc,scores in dict_sum_domain_diff.items():
        tsv.write(f"{acc}\t{scores}\n")
        
#For the non-domains:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_iupred_regions_confidence_scores.tsv","w") as tsv:
    tsv.write("Accession_number\tSum_of_confidence_scores\n")
    for acc,scores in dict_sum_region_diff.items():
        tsv.write(f"{acc}\t{scores}\n")

end_time=time.time()
elapsed_seconds = end_time - start_time
elapsed_hours = elapsed_seconds / 3600

print(f"Script completed in {elapsed_hours:.2f} hours.")