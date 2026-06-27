#Counting the ratio of the overlap between the predicted disordered regions and the annotated domains:
import pandas as pd

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
Aiupred_tsv=pd.read_csv("/home/guest/Internship/results/01_aiupred_disordered_regions.tsv",sep="\t")

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


#Create a dictionary for annotated regions: key: accession number, value: list of regions
pfam_data = {}

for pfam_index,pfam_row in Pfam_tsv.iterrows():
    pfam_acc=pfam_row["accession"]
    pfam_region=[pfam_row["start"], pfam_row["end"]]
    pfam_domain=pfam_row["domain"]
    if pfam_acc not in pfam_data:
        pfam_data[pfam_acc] = []
    pfam_data[pfam_acc].append({
        "region":pfam_region,
        "domain":pfam_domain})

#Create a dictionary for predicted regions by Iupred2a: key: accesion number, value:list of regions
iupred_data = {}
for iup_index,iup_row in Iupred_tsv.iterrows():
    iup_acc=iup_row["accession"]
    iup_region=[iup_row["start"],iup_row["end"]]
    #iup_region_length = iup_row["end"] - iup_row["start"] + 1
    if iup_acc not in iupred_data:
        iupred_data[iup_acc] = []
    iupred_data[iup_acc].append(iup_region)

#Create a dictionary for predicted regions by Aiupred: key: accesion number, value:list of regions
aiupred_data = {}
for aiup_index,aiup_row in Aiupred_tsv.iterrows():
    aiup_acc=aiup_row["accession"]
    aiup_region=[aiup_row["start"],aiup_row["end"]]
    #aiup_region_length = aiup_row["end"] - aiup_row["start"] + 1
    if aiup_acc not in aiupred_data:
        aiupred_data[aiup_acc] = []
    aiupred_data[aiup_acc].append(aiup_region)

##Iupred2a
#Create a counter for the overlapping regions

iup_overlap_above= 0
iup_overlap_below= 0

k= 0 #Counter for predicted regions that has no annotated PFAM domains

for acc, regions in iupred_data.items():
    if acc not in pfam_data:
        #print(f'{acc} does not have annotated PFAM domains')
        k += len(regions)
        continue
    for iupred_region in regions:
        overlap_above = False
        has_overlap= False
        found_overlap_region= None
        found_domain= None
        iup_region_length = iupred_region[1] - iupred_region[0] + 1
        for pfam_values in pfam_data[acc]:      #pfam_data[acc] refers to the list of regions that belongs to a certain accession number 
            pfam_region= pfam_values["region"]
            pfam_domain= pfam_values["domain"]
            overlap_result = overlap(pfam_region, iupred_region)  
            if overlap_result[1] != 0:
                has_overlap=True
                if found_overlap_region is None:
                    found_overlap_region = overlap_result[0]
                    found_domain= pfam_domain
                overlap_ratio= overlap_result[1]/ iup_region_length
#                 print("Overlapping regions:" acc , overlap_result)
#                 print(overlap_ratio)
#                 Counting the overlapping regions,and collect them in separate files based on the overlap ratio
                if overlap_ratio >= 0.5:
                    #print(f'{acc} PFAM: {pfam_region}, IUP: {iupred_region} overlaps more than 50%')
                    overlap_above= True
                    iup_overlap_above +=1
                    #print(acc,overlap_result[0],pfam_domain)    
                    iup_overlap_above_file.write(f"{acc}\t{overlap_result[0]}\t{pfam_domain}\n")
                    break
        if not overlap_above:
            iup_overlap_below +=1
            if has_overlap:
                iup_overlap_below_file.write(f"{acc}\t{found_overlap_region}\t{found_domain}\n")
            else:
                 iup_overlap_below_file.write(f"{acc}\t{iupred_region}\tNo_PFAM_overlap\n")
           
      
# print(f'{iup_overlap_above} regions overlap more than 50%')
# print(f'{iup_overlap_below} region overlap less than 50%')
# print(f'{k} region dont have PFAM at all')
# print(f'{sum([len(x) for x in iupred_data.values()])} regions total')        

iup_overlap_above_file.close()
iup_overlap_below_file.close()

##Aiup
#Create a counter for the overlapping regions

aiup_overlap_above= 0
aiup_overlap_below= 0

m= 0 #Counter for predicted regions that has no annotated PFAM domains

for acc, regions in aiupred_data.items():
    if acc not in pfam_data:
        #print(f'{acc} does not have annotated PFAM domains')
        m += len(regions)
        continue
    for aiupred_region in regions:
        overlap_above = False
        has_overlap= False
        found_overlap_region= None
        found_domain= None
        aiup_region_length = aiupred_region[1] - aiupred_region[0] + 1
        for pfam_values in pfam_data[acc]:      #pfam_data[acc] refers to the list of regions that belongs to a certain accession number 
            pfam_region= pfam_values["region"]
            pfam_domain= pfam_values["domain"]
            overlap_result = overlap(pfam_region, aiupred_region)
            if overlap_result[1] != 0:
                has_overlap=True
                if found_overlap_region is None:
                    found_overlap_region = overlap_result[0]
                    found_domain= pfam_domain
                overlap_ratio= overlap_result[1]/ aiup_region_length
#                 print("Overlapping regions:"acc , overlap_result)
#                 print(overlap_ratio)
#                 Counting the overlapping regions,and collect them in separate files based on the overlap ratio
                if overlap_ratio >= 0.5:
                    #print(f'{acc} PFAM: {pfam_region}, AIUP: {aiupred_region} overlaps more than 50%')
                    overlap_above= True
                    aiup_overlap_above +=1
                    #print(acc,overlap_result[0],pfam_domain)
                    aiup_overlap_above_file.write(f"{acc}\t{overlap_result[0]}\t{pfam_domain}\n")
                    break
        if not overlap_above:
            aiup_overlap_below +=1
            if has_overlap:
                print(acc, found_overlap_region, found_domain)
                aiup_overlap_below_file.write(f"{acc}\t{found_overlap_region}\t{found_domain}\n")
            else:
                 aiup_overlap_below_file.write(f"{acc}\t{aiupred_region}\tNo_PFAM_overlap\n")
                 
# print(f'{aiup_overlap_above} regions overlap more than 50%')
# print(f'{aiup_overlap_below} region overlap less than 50%')
# print(f'{m} region dont have PFAM at all')
# print(f'{sum([len(x) for x in aiupred_data.values()])} regions total')        

aiup_overlap_above_file.close()
aiup_overlap_below_file.close()

#Summarize the results in a txt file
with open("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_summary_number_of_overlapping_regions.txt","w") as file:
    file.write("Overlapping regions equal to and more than 50% Pfam vs Iupred2a: {0}\n".format(iup_overlap_above))
    file.write("Overlapping regions less than 50% Pfam vs Iupred2a: {0}\n".format(iup_overlap_below))
    file.write("Overlapping regions equal to and more than 50% Pfam vs AIUpred: {0}\n".format(aiup_overlap_above))
    file.write("Overlapping regions less than 50% Pfam vs AIUpred: {0}\n".format(aiup_overlap_below))
