#Visualize the disorder tendency of each residue in a given protein with a graph
#Visualize found examples for redox-sensitive disordered proteins

import iupred2a_lib
import matplotlib.pyplot as plt
import requests
import matplotlib.pyplot as plt
import iupred2a_lib

def get_sequence(accession):
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    
    fasta_file = requests.get(url)
    fasta_file.raise_for_status()
    
    fasta_text = fasta_file.text
    fasta_data=fasta_text.splitlines()
    sequence = "".join(fasta_data[1:])
    header= fasta_data[0]
    protein_name=header.split(" ", 1)[1].split(" OS=", 1)[0]
    
    return sequence,protein_name

def visualize(seq,name,accession):
    pred1 = iupred2a_lib.iupred(seq)[0]
    plt.plot(pred1, label="ox")

    pred2 = iupred2a_lib.iupred_redox(seq)[0]
    plt.plot(pred2, label="redox")
    plt.title(f"Redox-dependent prediction of\n{name}")
    plt.margins(y=0.3)
    plt.legend()
    plt.figtext(
    0.05, 0.02,
    "Ox: oxidized-state disorder scores "
    "Redox: reduced-state disorder scores",
    ha="left",
    fontsize=9)
    plt.savefig(f"/home/guest/Internship/results/Interproscan_Pfam/19_aiupred_redox_prediction_graph_{accession}.png",dpi=300, bbox_inches="tight")
    plt.show()
    plt.close()

# UniProt accession ID of two examples:
# Copper chaperone for superoxide dismutase: 014618
#Mitochondrial intermembrane space import and assembly protein 40: Q8N4Q1
accessions = ["O14618","Q8N4Q1"]

# Download and extract sequence
for accession in accessions:
    seq,name = get_sequence(accession)
# Visualize the sequence
    visualize(seq,name,accession)



