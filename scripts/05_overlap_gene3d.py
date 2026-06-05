#Counting the ratio of the overlap between the predicted disordered regions and the annotated domains:
import pandas as pd

#Calculating the length of the overlapping parts of the regions:
def overlap(x, y):
    overlap_start = max(x[0], y[0])
    overlap_end = min(x[-1], y[-1])
    overlap_length = max(0, overlap_end - overlap_start + 1)
    overlap_region = [overlap_start, overlap_end]
    return overlap_region,overlap_length

#Read the Gene3D tsv file
Gene3D_tsv = pd.read_csv("/home/guest/Internship/results/Interproscan_Gene3D/03_gene3d_filtered_UP000005640_9606.fasta.tsv",sep="\t")

#Read the Iupred tsv file
Iupred_tsv=pd.read_csv("/home/guest/Internship/results/01_iupred_disordered_regions.tsv",sep="\t")

#Read the Iupred tsv file
Aiupred_tsv=pd.read_csv("/home/guest/Internship/results/01_aiupred_disordered_regions.tsv",sep="\t")

#Create a tsv Iupred > 0.5 file and an <0.5

iup_overlap_above_file=open("/home/guest/Internship/results/Interproscan_Gene3D/05_gene3d_iupred_overlap_above_50.tsv","w")
iup_overlap_below_file=open("/home/guest/Internship/results/Interproscan_Gene3D/05_gene3d_iupred_overlap_below_50.tsv","w")

iup_overlap_above_file.write("Accession_number\tRegion\tDomain\n")
iup_overlap_below_file.write("Accession_number\tRegion\tDomain\n")

#Create a tsv AIUpred > 0.5 and an <0.5

aiup_overlap_above_file=open("/home/guest/Internship/results/Interproscan_Gene3D/05_gene3d_aiupred_overlap_above_50.tsv","w")
aiup_overlap_below_file=open("/home/guest/Internship/results/Interproscan_Gene3D/05_gene3d_aiupred_overlap_below_50.tsv","w")

aiup_overlap_above_file.write("Accession_number\tRegion\tDomain\n")
aiup_overlap_below_file.write("Accession_number\tRegion\tDomain\n")


#Create a dictionary for annotated regions: key: accession number, value: list of regions
gene3d_data = {}

for gene3d_index,gene3d_row in Gene3D_tsv.iterrows():
    gene3d_acc=gene3d_row["accession"]
    gene3d_region=[gene3d_row["start"], gene3d_row["end"]]
    gene3d_domain=gene3d_row["domain"]
    if gene3d_acc not in gene3d_data:
        gene3d_data[gene3d_acc] = []
    gene3d_data[gene3d_acc].append({
        "region":gene3d_region,
        "domain":gene3d_domain})

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

#Create a counter for the overlapping regions

iup_overlap_above= 0
iup_overlap_below= 0

k= 0 #Counter for predicted regions that has no annotated PFAM domains

for acc, regions in iupred_data.items():
    if acc not in gene3d_data:
        #print(f'{acc} does not have annotated PFAM domains')
        k += len(regions)
        continue
    for iupred_region in regions:
        does_overlap = False
        iup_region_length = iupred_region[1] - iupred_region[0] + 1
        for gene3d_values in gene3d_data[acc]:      #gene3d_data[acc] refers to the list of regions that belongs to a certain accession number 
            gene3d_region= gene3d_values["region"]
            gene3d_domain= gene3d_values["domain"]
            overlap_result = overlap(gene3d_region, iupred_region)
            if overlap_result[1] != 0:
                overlap_ratio= overlap_result[1]/ iup_region_length
#                 print("Overlapping regions:" acc , overlap_result)
#                 print(overlap_ratio)
#                 Counting the overlapping regions,and collect them in separate files based on the overlap ratio
                if overlap_ratio >= 0.5:
                    #print(f'{acc} PFAM: {gene3d_region}, IUP: {iupred_region} overlaps more than 50%')
                    does_overlap= True
                    iup_overlap_above +=1
                    #print(acc,overlap_result[0],gene3d_domain)    
                    iup_overlap_above_file.write(f"{acc}\t{overlap_result[0]}\t{gene3d_domain}\n")
                    break
        if not does_overlap:
            iup_overlap_below +=1
            iup_overlap_below_file.write(f"{acc}\t{overlap_result[0]}\t{gene3d_domain}\n")
  
               
# print(f'{iup_overlap_above} regions overlap more than 50%')
# print(f'{iup_overlap_below} region overlap less than 50%')
# print(f'{k} region dont have PFAM at all')
# print(f'{sum([len(x) for x in iupred_data.values()])} regions total')        

iup_overlap_above_file.close()
iup_overlap_below_file.close()

#Create a counter for the overlapping regions

aiup_overlap_above= 0
aiup_overlap_below= 0

m= 0 #Counter for predicted regions that has no annotated PFAM domains

for acc, regions in aiupred_data.items():
    if acc not in gene3d_data:
        #print(f'{acc} does not have annotated PFAM domains')
        m += len(regions)
        continue
    for aiupred_region in regions:
        does_overlap = False
        aiup_region_length = aiupred_region[1] - aiupred_region[0] + 1
        for gene3d_values in gene3d_data[acc]:      #gene3d_data[acc] refers to the list of regions that belongs to a certain accession number 
            gene3d_region= gene3d_values["region"]
            gene3d_domain= gene3d_values["domain"]
            overlap_result = overlap(gene3d_region, aiupred_region)
            if overlap_result[1] != 0:
                overlap_ratio= overlap_result[1]/ aiup_region_length
#                 print("Overlapping regions:"acc , overlap_result)
#                 print(overlap_ratio)
#                 Counting the overlapping regions,and collect them in separate files based on the overlap ratio
                if overlap_ratio >= 0.5:
                    #print(f'{acc} PFAM: {gene3d_region}, AIUP: {aiupred_region} overlaps more than 50%')
                    does_overlap= True
                    aiup_overlap_above +=1
                    #print(acc,overlap_result[0],gene3d_domain)
                    aiup_overlap_above_file.write(f"{acc}\t{overlap_result[0]}\t{gene3d_domain}\n")
                    break
        if not does_overlap:
            aiup_overlap_below +=1
            aiup_overlap_below_file.write(f"{acc}\t{overlap_result[0]}\t{gene3d_domain}\n")
                 
# print(f'{aiup_overlap_above} regions overlap more than 50%')
# print(f'{aiup_overlap_below} region overlap less than 50%')
# print(f'{m} region dont have PFAM at all')
# print(f'{sum([len(x) for x in aiupred_data.values()])} regions total')        

aiup_overlap_above_file.close()
aiup_overlap_below_file.close()

#Summarize the results in a txt file
with open("/home/guest/Internship/results/Interproscan_Gene3D/05_gene3d_summary_number_of_overlapping_regions.txt","w") as file:
    file.write("Overlapping regions equal to and more than 50% Gene3D vs Iupred2a: {0}\n".format(iup_overlap_above))
    file.write("Overlapping regions less than 50% Gene3D vs Iupred2a: {0}\n".format(iup_overlap_below))
    file.write("Overlapping regions equal to and more than 50% Gene3D vs AIUpred: {0}\n".format(aiup_overlap_above))
    file.write("Overlapping regions less than 50% Gene3D vs AIUpred: {0}\n".format(aiup_overlap_below))
