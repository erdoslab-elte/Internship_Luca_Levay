#Calculation of redox-state dependent disorder prediction scores for the domain and non-domain regions identified by Iupred2a
#Comparing the distribution of the prediction scores between domain and non-domain regions

from aiupred import AIUPred
import pandas as pd
import time

start_time=time.time()

#Creating functions to calculate the redox-state dependent disorder prediction scores:
predictor = AIUPred()

def aiupred_score(seq):
    prediction = predictor.predict_disorder(seq)
    return prediction

def aiupred_redox_score(seq):
    return aiupred_score(seq.replace("C", "S"))

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

#Read in the Aiupred predicte region files:
##Domains: regions that overlap more than 50% iwth annotated domains
aiup_domain_file=pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_aiupred_overlap_above_50.tsv",sep="\t")

##Non-domains/Regions:called here regions: regions that overlap less than 50% with annotated domains:
aiup_region_file=pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_aiupred_overlap_below_50.tsv",sep="\t")

#Create a list of the accesion numbers of the predicted redox-state dependent disordered regions by Aiupred:
aiup_domain_data = []
for _,aiup_row in aiup_domain_file.iterrows():
    aiup_acc=aiup_row["Accession_number"]
    if aiup_acc not in aiup_domain_data:
        aiup_domain_data.append(aiup_acc)
#print(f"aiup domain data: {aiup_domain_data[:3]}")

aiup_region_data = []
for _,aiup_row in aiup_region_file.iterrows():
    aiup_acc=aiup_row["Accession_number"]
    if aiup_acc not in aiup_region_data:
        aiup_region_data.append(aiup_acc)  
#print(f"aiup region data: {aiup_region_data[:3]}")

#Save the predicted redox sensitive regions' sequences in a dictionary
aiup_domain_found_ids={}
aiup_region_found_ids={}

#n=0
for id,seq in fasta_data.items():
    for aiup_dom_acc in aiup_domain_data:
        if id == aiup_dom_acc and aiup_dom_acc not in aiup_domain_found_ids:
            if len(seq) <= 1000:
                aiup_domain_found_ids[aiup_dom_acc]=""
                aiup_domain_found_ids[aiup_dom_acc]=seq
            else:
                continue
    for aiup_reg_acc in aiup_region_data:
        if id == aiup_reg_acc and aiup_reg_acc not in aiup_region_found_ids:
            if len(seq) <= 1000:
                aiup_region_found_ids[aiup_reg_acc]=""
                aiup_region_found_ids[aiup_reg_acc]=seq
            else:
                continue      
    # n += 1
    # if n >15:
    #     break

# print(f"found aiup domains: {aiup_domain_found_ids}")
# print(f"found aiup regions: {aiup_region_found_ids}")

#Calculating the aiupred scores for the sequences: calculating the difference of the two scores for a residual,then adding them together
#Calculating for DOMAINS predicted by AIUPRED:

# #p=0

aiup_domain_score=[]
aiup_domain_redox_score=[]
list_sum_domain_diff=[]
dict_sum_domain_diff={}

for aiup_acc,aiup_seq in aiup_domain_found_ids.items():
    #m=0
    #print(f"aiup acc: {aiup_acc}")
    #print(aiup_seq)
    aiup_score=aiupred_score(aiup_seq)
    aiup_redox_score=aiupred_redox_score(aiup_seq)
    #print(f"aiup score: {aiup_score[0:3]}")
    # print(f"aiup redox score: {aiup_redox_score[0:3]}")
    list_diff=[] #creating a list for the calculated differences
    dict_sum_domain_diff[aiup_acc]=0
    redox_score_list=[]
    for score in aiup_score:
        #print(f"score: {score}")
        for redox_score in aiup_redox_score:
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
        # if m >= 2:
        #     break
    #print(f"list of differences: {list_diff}")
    sum_diff=sum(list_diff)
    #print(f"sum of differences: {sum_diff}")
    dict_sum_domain_diff[aiup_acc]=sum_diff
    list_sum_domain_diff.append(sum_diff)
    # p +=1
    # if p >=2:
    #     break
    # if m>=2:
    #     break

#print(f"sum of differences in domain sequences: {list_sum_domain_diff}")
#print(f"dictionary of domains and score differeneces: {dict_sum_domain_diff}")

##Calculating for REGIONS predicted by AIUPRED:

#p=0

aiup_region_score=[]
aiup_region_redox_score=[]
list_sum_region_diff=[]
dict_sum_region_diff={}

for aiup_acc,aiup_seq in aiup_region_found_ids.items():
    #m=0
    #print(f"aiup acc: {aiup_acc}")
    #print(aiup_seq)
    aiup_score=aiupred_score(aiup_seq)
    aiup_redox_score=aiupred_redox_score(aiup_seq)
    #print(f"aiup score: {aiup_score[0:3]}")
    #print(f"aiup redox score: {aiup_redox_score[0:3]}")
    list_diff=[] #creating a list for the calculated differences
    dict_sum_region_diff[aiup_acc]=0
    redox_score_list=[]
    for score in aiup_score:
        #print(f"score: {score}")
        for redox_score in aiup_redox_score:
            if redox_score in redox_score_list: #so that it just compares it with one redox score not all of them
                continue
            else:
                redox_score_list.append(redox_score)
                print(f"redox score: {redox_score}")
                score_diff=abs(score-redox_score)
                #print(f"score difference :{score_diff}")
                list_diff.append(score_diff)
                #m += 1
                break            
        # if m >= 2:
        #     break
    #print(f"list of differences: {list_diff}")
    sum_diff=sum(list_diff)
    #print(f"sum of differences: {sum_diff}")
    dict_sum_region_diff[aiup_acc]=sum_diff
    list_sum_region_diff.append(sum_diff)
    # p +=1
    # if p >=2:
    #     break
    # if m>=2:
    #     break
# print(f"sum of differences in region sequences: {list_sum_region_diff}")
# print(f"dictionary of domains and score differeneces: {dict_sum_region_diff}")

#Creat text files to store these lists:
#For the domains:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_aiupred_domains_confidence_scores.txt", "w") as file:
    for numbers in list_sum_domain_diff:
        file.write(f"{numbers}\n")

#For the regions:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_aiupred_regions_confidence_scores.txt", "w") as file:
    for numbers in list_sum_region_diff:
        file.write(f"{numbers}\n")

#Create tsv files for the accession numbers and the sum of the confidence scores belonging to them
#For the domains:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_aiupred_domains_confidence_scores.tsv","w") as tsv:
    tsv.write("Accession_number\tSum_of_confidence_scores\n")
    for acc,scores in dict_sum_domain_diff.items():
        tsv.write(f"{acc}\t{scores}\n")
        
#For the non-domains:
with open("/home/guest/Internship/results/Interproscan_Pfam/14_pfam_aiupred_regions_confidence_scores.tsv","w") as tsv:
    tsv.write("Accession_number\tSum_of_confidence_scores\n")
    for acc,scores in dict_sum_region_diff.items():
        tsv.write(f"{acc}\t{scores}\n")

end_time=time.time()
elapsed_seconds = end_time - start_time

elapsed_hours = elapsed_seconds / 3600

print(f"Script completed in {elapsed_hours:.2f} hours.")