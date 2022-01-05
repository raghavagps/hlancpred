# hlancpred
## Introduction
HLAncPred is developed for predicting and scanning the peptides with the ability to bind non-classical class-I HLA alleles such as HLA-G\*01:01,HLA-G\*01:03,HLA-G\*01:04, HLA-E\*01:01, HLA-E\*01:03 . More information on HLAncPred is available from its web-server https://webs.iiitd.edu.in/raghava/hlancpred/stand.html. This page provides information about stnadalone version of HLAncPred. Please read/cite the content about the HLAncPred for complete information including algorithm behind HLAncPred.

## Models
In this program, five best models devoted to each allele, have been implemented for predicting the non-classical class-I HLA binder peptides. The models were trained on experimentally verified binders and randomly generated non-binders using Swiss-Prot database.

## Modules/Jobs
This program implements two modules (job types); i) Predict: for prediction of non-classical class-I HLA allele binder peptides, ii) Scan: for creating all possible overlapping peptides of given length (window) and computing binding  potential (score) of the overlapping peptides.

## Minimum USAGE
Minimum usage is "python3 hlancpred.py -i peptide.fa -a G0101," where peptide.fa is a input fasta file, and 'G0101' is the name for HLA-G\*01:01 against which the binders will be predicted. This will predict the binding ability of sequence in fasta format against HLA-G\*01:01. It will use other parameters by default. It will save output in "outfile.csv" in CSV (comma seperated variables).

-------------------------------------------------------------------------------------------------------------

## Full Usage
Following is complete list of all options, you may get these options by "python hlancpred.py -h"


## Usage
#### hlancpred.py [-h] -i INPUT -a {G0101,G0103,G0104,E0101,E0103} [-o OUTPUT] [-j {1,2}] [-w {8,9,10,11,12,13,14,15}] [-d {1,2}]


### Please provide following arguments
======================================
-------------------------------------------------------------------------

## Optional arguments:

  ##### -h, --help            
  >show this help message and exit
  
  ##### -i INPUT, --input INPUT
  >Input: protein or peptide sequence in FASTA format or single sequence per line in single letter code
                        
  ##### -a {G0101,G0103,G0104,E0101,E0103}, --allele {G0101,G0103,G0104,E0101,E0103}
  >Please provide the name of allele for the prediction of binder peptides
                        
  ##### -o OUTPUT, --output OUTPUT
  >Output: File for saving results by default outfile.csv
  ##### -j {1,2}, --job {1,2}
  >Job Type: 1:predict, and 2:scan, by default 1
  ##### -w {8,9,10,11,12,13,14,15}, --winleng {8,9,10,11,12,13,14,15}
  >Window Length: 8 to 35 (scan mode only), by default 9
  ##### -d {1,2}, --display {1,2}
  >Display: 1:Only binders, 2: All peptides, by default 
  
-------------------------------------------------------------------------


## **Input File:** 
It allow users to provide input in two format; i) FASTA format (standard) and ii) Simple Format. In case of simple format, file should have one peptide sequence in a single line in single letter code (eg. peptide.seq). 


## **Note:**
>1: In case of predict and design module (job), the length of peptide should be upto 15 amino acids. If a sequence with length more than 15 will be provided, the program will take first 15 residues, and ignore the rest. In case of scan module, minimum length of protein/peptide sequence should be equal to window length (pattern), see peptide.fa.
>2: Program will ignore peptides having length less than 8 residues (e.g., protein.fa).

## **Output File:** 
Program will save the results in the CSV format, in case user do not provide output file name, it will be stored in "outfile.csv".

## **Allele:** 
The program needs the name of the allele as shown in the usage, it could be G0101, G0103, G0104, E0101, and E0103.

## HLAncPred Package Files
It contantain following files, brief descript of these files given below

- INSTALLATION  			: Installations instructions

- LICENSE       			: License information

- README.md     			: This file provide information about this package

- Models           		: This folder comprises for five models devoted to five alleles

- hlancpred.py 			: Main python program 

- peptide.fa			: Example file contain peptide sequenaces in FASTA format

- peptide.seq			: Example file contain peptide sequenaces in simple format

- protein.fa			: Example file contain protein sequenaces in FASTA format 

- example_predict_output.csv	: Example output file for predict module

- example_scan_output.csv		: Example output file for scan module
-------------------------------------------------------------------
~Address for contact~
Prof. G. P. S. Raghava, Head Department of Computational Biology,            
Indraprastha Institute of Information Technology (IIIT), 
Okhla Phase III, New Delhi 110020 ; Phone:+91-11-26907444; 
Email: raghava@iiitd.ac.in  Web: http://webs.iiitd.edu.in/raghava/
--------------------------------------------------------------------
