**README**

**Project overview**



This repository contains the scripts, results, and documentation from my internship project on "Bioinformatic analysis and prediction of redox potential-sensitive disordered protein regions"



The project mainly uses IUPred2A and AIUPred to predict disordered and redox-state-dependent regions. The analysis was performed on the human proteome, predicted regions were compared with annotated protein domains using Pfam and Gene3D, followed by functional annotation and Gene Ontology (GO) analysis.



Scripts were written in Python 3.14.7 using Visual Studio Code.



Main Python packages:

pandas

matplotlib

numpy

time

gzip

ast

iupred.iupred2a\_lib

aiupred

goatools

goscripts

contextlib

requests



**Repository structure:**



**Internship\_Luca\_Levay/**

**│**

**├── data/**

**│   └── data\_list.txt**

**│**

**├── results/**

**│   ├── Interproscan\_Pfam/**

**│   └── Interproscan\_Gene3D/**

**│**

**└── scripts/**



**Repository description:**



**data/**



Contains information about the input data used in the project.



required\_files\_for\_the\_project.txt – list of the datasets required for the analysis.



The original datasets cannot be included in the repository because some of them are too large.



**results/**



Contains the results generated throughout the project.



Results independent of annotation databases

Interproscan\_Pfam/ – results based on the Pfam-annotated human proteome

Interproscan\_Gene3D/ – results based on the Gene3D-annotated human proteome



**scripts/**



Contains all Python scripts created during the internship.





**Software and databases**



The project used the following main resources:



IUPred2A – prediction of intrinsically disordered regions and redox-state-dependent disorder

AIUPred – prediction of intrinsically disordered and redox-sensitive regions

InterProScan6 – protein domain annotation

Pfam – protein domain database

Gene3D – protein domain database

PANTHER – functional annotation and GO analysis

Gene Ontology – functional annotation and enrichment analysis

UniProt – protein sequences and ID mapping

REVIGO – visualization of GO terms



**Important notes**

The original human proteome and annotation datasets are not included in the repository because of their large file sizes.

The analyses involving AIUPred were restricted to protein sequences of ≤1000 amino acids where indicated, because of computational time.

A 50% overlap threshold was used to distinguish domain and non-domain regions.

From the comparison of Pfam and Gene3D, the subsequent analysis focused on Pfam-annotated regions.

For the final selection of potential redox-sensitive regions, the AIUPred 0.7 cutoff was used.

