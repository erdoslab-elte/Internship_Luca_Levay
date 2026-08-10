#FONTOS: 1000 AS SZEKVENCIA SZŰRÉSSEL MI LEGYEN?
#Normalizing the predicted scores for the regions identified using the modified Aiupred get_redox_region function

import pandas as pd

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

#Read in the file containing the iroginal aiupred predicted regions that overlap less than 50% with the annotated doamins (by Pfam)
#We will use this to filter on the non-domain regions
below_file = pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_aiupred_overlap_below_50.tsv",sep="\t", usecols=["Accession_number", "Region"])


#Read in all the modified files:



#Read in the multifasta file (with the defined multifasta function)
human_data = multi_fasta_reader('/home/guest/Internship/data/UP000005640_9606.fasta')

#Filter the aiupred_below_file so it contains only proteins that consist of less than 1000 aminoacids
#Because this is how it was done in the 14_pfam_aiup_predicted_score script

filtered_below_list=[]

for _,row in below_file.iterrows():
    acc=row["Accession_number"]
    region=row["Region"]
    for id,seq in human_data.items():
        if acc == id:
            if len(seq) <= 1000:
                filtered_below_list.append({
                    "Accession_number": acc,
                    "Region": region
                })

filtered_below_file=pd.DataFrame(filtered_below_list)

for _, m_row in modified_file.iterrows():
    m_acc=m_row["Accession_number"]
    m_region=[m_row["Start"],m_row["End"]]
    m_start=m_row
    for _, b_row in filtered_below_file.iterrows():
        b_acc=b_row["Accession_number"]
        b_region=b_row["Region"]
        if b_acc == m_acc and m_region == b_region:
            region
            aiupred_score
#Creating a data frame:


