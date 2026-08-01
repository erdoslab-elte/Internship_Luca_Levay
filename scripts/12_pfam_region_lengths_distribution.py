#Analyzing the distirbuiton of the lengths of the predicted redox sensitive disordered regions
#Pfam
import pandas as pd
import ast
import numpy as np
import matplotlib.pyplot as plt

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

def create_acc_list(file):
    list=[]
    for _,row in file.iterrows():
        acc=row["Accession_number"]
        if acc not in list:
            list.append(acc)
    return list

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

def length_calculation(found_ids_dict,tsv_file,result_file):
    list_lengths=[]
    dict_lengths={}
    for _,row in tsv_file.iterrows():
        acc1=row['Accession_number']
        region=ast.literal_eval(row["Region"])
        start,end=region
        length=end - start +1
        for acc2,_ in found_ids_dict.items():
            if acc1 == acc2:
                list_lengths.append(length)
                if acc1 not in dict_lengths:
                    dict_lengths[acc1]=[]
                dict_lengths[acc1].append(length)
                result_file.write(f"{acc1}\t{region}\t{length}\n")
    return list_lengths, dict_lengths

#Read in the multifasta file (with the defined multifasta function)
fasta_data = multi_fasta_reader('/home/guest/Internship/data/UP000005640_9606.fasta')

#Reading in the files: predicted by Iupred2a and Aiupred overlapping more than 50% and less than 50%
#Iupred2a
iup_domain_file = pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_iupred_overlap_above_50.tsv",sep="\t")
iup_segment_file = pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_iupred_overlap_below_50.tsv",sep="\t")

#Aiupred
aiup_domain_file = pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_aiupred_overlap_above_50.tsv",sep="\t")
aiup_segment_file = pd.read_csv("/home/guest/Internship/results/Interproscan_Pfam/04_pfam_aiupred_overlap_below_50.tsv",sep="\t")

#Create a list of the accesion numbers of the predicted domains and segments by Iupred2a:
iup_domain_list=create_acc_list(iup_domain_file)

iup_segment_list=create_acc_list(iup_segment_file)

#Create a list of the accesion numbers of the predicted domains and segments by Aiupred:
aiup_domain_list=create_acc_list(aiup_domain_file)

aiup_segment_list=create_acc_list(aiup_segment_file)

#Save the redox sensitive domains and segments sequences in a dictionary predicted by Iupred2a:
iup_domain_found_ids_dict=create_id_dict(fasta_data,iup_domain_list)
iup_segment_found_ids_dict=create_id_dict(fasta_data,iup_segment_list)


#Save the redox sensitive domains segments sequences in a dictionary predicted by Aiupred:
aiup_domain_found_ids_dict=create_id_dict(fasta_data,aiup_domain_list)
aiup_segment_found_ids_dict=create_id_dict(fasta_data,aiup_segment_list)

##Calculating segment lengths:
#IUPRED2A

#For domains:
#Create a tsv file,and save the lengths in a list:
iup_domain_lengths_file=open("/home/guest/Internship/results/Interproscan_Pfam/12_pfam_iupred_domain_lengths.tsv","w")
iup_domain_lengths_file.write("Accession_number\tRegion\tRegion_length\n")

list_iup_domain_lengths=length_calculation(iup_domain_found_ids_dict,iup_domain_file,iup_domain_lengths_file)[0]

iup_domain_lengths_file.close()

#For non-domains/segments:
#Create a tsv file,and save the lengths in a list:
iup_segment_lengths_file=open("/home/guest/Internship/results/Interproscan_Pfam/12_pfam_iupred_segment_lengths.tsv","w")
iup_segment_lengths_file.write("Accession_number\tRegion\tRegion_length\n")

list_iup_segment_lengths=length_calculation(iup_segment_found_ids_dict,iup_segment_file,iup_segment_lengths_file)[0]

iup_segment_lengths_file.close()

#Aiupred
#For domains:
#Create a tsv file,and save the lengths in a list:
aiup_domain_lengths_file=open("/home/guest/Internship/results/Interproscan_Pfam/12_pfam_aiupred_domain_lengths.tsv","w")
aiup_domain_lengths_file.write("Accession_number\tRegion\tRegion_length\n")

list_aiup_domain_lengths=length_calculation(aiup_domain_found_ids_dict,aiup_domain_file,aiup_domain_lengths_file)[0]

aiup_domain_lengths_file.close()

#For non-domains/segments:
#Create a tsv file,and save the lengths in a list:
aiup_segment_lengths_file=open("/home/guest/Internship/results/Interproscan_Pfam/12_pfam_aiupred_segment_lengths.tsv","w")
aiup_segment_lengths_file.write("Accession_number\tRegion\tRegion_length\n")

list_aiup_segment_lengths=length_calculation(aiup_segment_found_ids_dict,aiup_segment_file,aiup_segment_lengths_file)[0]

aiup_segment_lengths_file.close()

# #Filtering the outliers:IQR method (boxplot rule)
data = [list_iup_domain_lengths, list_iup_segment_lengths, list_aiup_domain_lengths, list_aiup_segment_lengths]

filtered_data=[]

for list in data: 
    q1 = np.percentile(list, 20)
    q3 = np.percentile(list, 80)
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    filtered_list = [x for x in list if lower_bound <= x <= upper_bound]
    filtered_data.append(filtered_list)

#Creating violin plots 
plt.figure(figsize=(8, 6))
vp=plt.violinplot(filtered_data,showmeans=True)
plt.xticks([1, 2, 3, 4],
           ["Iupred2a domains", "Iupred2a non-domains", "Aiupred domains", "Aiupred non-domains"])
plt.ylabel("Region lengths (AA number)")
plt.title("Comparison of distributions of region lengths")

colors=['blue','red','blue','red']

for i, body in enumerate(vp["bodies"]):
    body.set_facecolor(colors[i])
    body.set_edgecolor('black')
    body.set_alpha(0.5)
    body.set_linewidth(1.5)

quantilies = [np.percentile(d, [25, 50, 75]) for d in filtered_data]

for i, q in enumerate(quantilies):
    plt.scatter([i+1]*3, q, color="white", s=20,marker="s")

vp['cmins'].set_color('black')
vp['cmaxes'].set_color('black')
vp['cbars'].set_color('black')
vp['cmeans'].set_color('white')
plt.savefig("/home/guest/Internship/results/Interproscan_Pfam/12_pfam_region_lengths_violinplot.png", dpi=300, bbox_inches="tight")
plt.show()

#print(vp.keys())

# #Creating histograms for that:

# fig, axes = plt.subplots(2, 2, figsize=(10, 8))

# axes[0, 0].hist(filtered_data[0], bins=10)
# axes[0, 0].set_title("Iupred2a domains")

# axes[0, 1].hist(filtered_data[1], bins=10)
# axes[0, 1].set_title("Iupred2 a non-domains")

# axes[1, 0].hist(filtered_data[2], bins=10)
# axes[1, 0].set_title("Aiupred domains")

# axes[1, 1].hist(filtered_data[3], bins=10)
# axes[1, 1].set_title("Aiupred non-domains")

# plt.tight_layout()
# plt.show()