#Calculating overlapping
import pandas as pd
from datetime import datetime

start_time = datetime.now()

#Calculating the length of the overlapping parts of the regions:
def overlap(x, y):
    overlap_start = max(x[0], y[0])
    overlap_end = min(x[-1], y[-1])
    overlap_length = max(0, overlap_end - overlap_start + 1)
    overlap_region = [overlap_start, overlap_end]
    return overlap_region,overlap_length

#Read the Pfam tsv file
Pfam_tsv = pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/02_pfam_filtered_UP000005640_9606.fasta.tsv",sep="\t")

#Read the Iupred tsv file
Iupred_tsv=pd.read_csv("/home/guest/Internship/results/01_iupred_disordered_regions.tsv",sep="\t")

#Read the Iupred tsv file
AIUpred_tsv=pd.read_csv("/home/guest/Internship/results/01_aiupred_disordered_regions.tsv",sep="\t")

#n=0  #Counter for the foor loop

#Create a counter for the overlapping regions

iup_overlap_above= 0
iup_overlap_below= 0

aiup_overlap_above= 0
aiup_overlap_below= 0

#Create a tsv Iupred > 0.5 file and an <0.5

iup_overlap_above_file=open("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_iupred_overlap_above_50.tsv","w")
iup_overlap_below_file=open("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_iupred_overlap_below_50.tsv","w")

iup_overlap_above_file.write("Accession_number\tRegion\tDomain\n")
iup_overlap_below_file.write("Accession_number\tRegion\tDomain\n")

#Create a tsv AIUpred > 0.5 and an <0.5

aiup_overlap_above_file=open("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_aiupred_overlap_above_50.tsv","w")
aiup_overlap_below_file=open("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_aiupred_overlap_below_50.tsv","w")

aiup_overlap_above_file.write("Accession_number\tRegion\tDomain\n")
aiup_overlap_below_file.write("Accession_number\tRegion\tDomain\n")

#Creating two lists to keep track on the already found domains:
#iup_found_domains=[]
#aiup_found_domains=[]

#Counting the ratio of the overlap between the predicted disordered regions and the annotated domains:
for pfam_index,pfam_row in Pfam_tsv.iterrows():
    pfam_acc=pfam_row["accession"]
    pfam_region=[pfam_row["start"], pfam_row["end"]]
    pfam_domain=pfam_row["domain"]
    #Comparing with Iupred2a results
    for iup_index,iup_row in Iupred_tsv.iterrows():
        iup_acc=iup_row["accession"]
        iup_region=[iup_row["start"],iup_row["end"]]
        iup_region_length = iup_row["end"] - iup_row["start"] + 1
        if pfam_acc == iup_acc:
            overlap_result=overlap(pfam_region,iup_region)
            if overlap_result[1] != 0:
                overlap_ratio= overlap_result[1]/ iup_region_length
                #print("Overlapping regions:" , pfam_acc , overlap_result)
                #print(overlap_ratio)
                #Counting the overlapping regions,and collect them in separate files based on the overlap ratio
                if overlap_ratio >= 0.5: #and pfam_domain not in iup_found_domains:
                    #iup_found_domains.append(pfam_domain)
                    iup_overlap_above+=1
                    iup_overlap_above_file.write(f"{pfam_acc}\t{overlap_result[0]}\t{pfam_domain}\n")
                elif overlap_ratio < 0.5: #and pfam_domain not in iup_found_domains
                    #iup_found_domains.append(pfam_domain)
                    iup_overlap_below+=1
                    iup_overlap_below_file.write(f"{pfam_acc}\t{overlap_result[0]}\t{pfam_domain}\n")
    #Comparing with AIUpred results:
    for aiup_index,aiup_row in AIUpred_tsv.iterrows():
        aiup_acc=aiup_row["accession"]
        aiup_region=[aiup_row["start"],aiup_row["end"]]
        aiup_region_length = aiup_row["end"] - aiup_row["start"] + 1
        if pfam_acc == aiup_acc:
            overlap_result=overlap(pfam_region,aiup_region)
            if overlap_result[1] != 0:
                overlap_ratio= overlap_result[1]/ aiup_region_length
                #print("Overlapping regions:" , pfam_acc , overlap_result)
                #print(overlap_ratio)
                if overlap_ratio >= 0.5: #and pfam_domain not in aiup_found_domains:
                    #aiup_found_domains.append(pfam_domain)
                    aiup_overlap_above+=1
                    aiup_overlap_above_file.write(f"{pfam_acc}\t{overlap_result[0]}\t{pfam_domain}\n")
                elif overlap_ratio < 0.5: #and pfam_domain not in iup_found_domains:
                    #aiup_found_domains.append(pfam_domain)
                    aiup_overlap_below+=1
                    aiup_overlap_below_file.write(f"{pfam_acc}\t{overlap_result[0]}\t{pfam_domain}\n")                
    #n += 1
    #if n>20:
        #break

iup_overlap_above_file.close()
iup_overlap_below_file.close()
aiup_overlap_above_file.close()
aiup_overlap_below_file.close()

#Summarize the results in a txt file
with open("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_number_of_overlapping_regions.txt","w") as file:
    file.write("Overlapping equal to and more than 50% Pfam vs Iupred2a: {0}\n".format(iup_overlap_above))
    file.write("Overlapping less than 50% Pfam vs Iupred2a: {0}\n".format(iup_overlap_below))
    file.write("Overlapping equal to and more than 50% Pfam vs AIUpred: {0}\n".format(aiup_overlap_above))
    file.write("Overlapping less than 50% Pfam vs AIUpred: {0}\n".format(aiup_overlap_below))

end_time=datetime.now()
duration= end_time - start_time
print("Started: ",start_time)
print("Finished: ", end_time)
print("Duration: ", duration)