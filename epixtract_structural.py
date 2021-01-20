#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 14:21:51 2020
@author: rocfarriolduran
"""
### EPITOPE EXTRACTOR
## GOAL: extract linear epitopes in tabulated data either from linear or structural prediction origin.
## Variant for structural Discotope predictions at: 
## https://services.healthtech.dtu.dk/service.php?DiscoTope-2.0
## create epitopes from continous residues that surpass a given quality threshold.


# IMPORTS
import sys
import csv
import os.path
from pathlib import Path
import pandas as pd

# HELP
h = '''
    To right usage of this script:
        $ python3 epitetons.py 
    The files in use have to be provided by stating "location/file_name.csv"
    to the input questions that appear in the console.
    <file_name> should be a .csv separated by ","
    OPTIONAL 
    The script returns a <file_name>_out.xlsx as output.
    You can also do the following (if using Mac/Linux OS):
        $ chmod +x epitetons.py 
        $ ./epitetons.py  <location/file_name.csv>
    You need python 3 installed in your computer!!!
    '''


# FILE EXISTANCE
def fileExist(file):
    if file!="":
        if Path(file).is_file():
            return True
        elif file=='-h':
            print(h)
        else:
            print("\n>>>>>>> "+file + " File not exist or is not accessible\n")
            return False
    else:
        return False


# FILE UPLOAD
data_file=""
num_args=len(sys.argv)


if num_args>=2:
    data_file = sys.argv[1]

while not fileExist(data_file):
    data_file = input(''' Provide location and name of the file dataset.
        For example: data/file_test.csv
            (-h for help)
        Here: ''')

print("Input file: " + data_file)

# NAME THE OUTPUT FILE
extension = os.path.splitext(data_file)[1]
lenextension=len(extension)
nameOutFile=data_file[:-lenextension]+"_Out"+extension
out_file=open(nameOutFile, "w")
out_file.writelines("Epitope_id,Epitope_seq\n")

# READ THE CSV INPUT & ADD HEADERS TO DISCOTOPE DATAFRAME
'''# chain_id	residue_id	residue_name	contact_number	propensity_score	discotope_score	status'''

prediction = pd.read_csv(data_file, 
			 sep='\t', 
			 names = ["chain_id", "residue_id", "residue_name", "contact_number", "propensity_score", "discotope_score", "status"])

#f=open(data_file, "r", encoding = 'utf-8-sig')
#inputFile = csv.reader(f, delimiter='\t')


# PIPELINE

## 3-LETTER TO 1-LETTER AA SCRIPT

# Create dictionary
aa_dict = {"ALA" : "A",
           "ARG" : "R",
           "ASN" : "N",
           "ASP" : "D",
           "CYS" : "C",
           "GLN" : "Q", 
           "GLU" : "E",
           "GLY" : "G",
           "HIS" : "H",
           "ILE" : "I",
           "LEU" : "L",
           "LYS" : "K",
           "MET" : "M",
           "PHE" : "F",
           "PRO" : "P",
           "SER" : "S",
           "THR" : "T",
           "TRP" : "W",
           "TYR" : "Y",
           "VAL" : "V"}

# Extract residues to convert
residues = prediction["residue_name"] #column from df as series
           
# Use dictionary to convert 3LAA to 1LAA    
prediction["aa"]=residues.map(aa_dict)
    
# Filter by discotope_score threshold at -3.7
scored = prediction[prediction.discotope_score >= -3.7]

# CONSECUTIVE RESIDUES EXTRACTOR
# PREDEFINE SCRIPT VARIABLES
y = 0 			# counter of registers read from in input data (iterator).
positionAnt = 0 	# position of the previous amino acid. 
AminoAcidAnt = ""	# 
Epitope_seq = ""	# epitope sequence to accumulate and extract
Epitope_id=0		# epitope identifier	
consecutius=1		# consecutive accumulator.


#for index,row in df.iterrows():
#    print("{0} has name: {1}".format(index,row["name"]))

'''# chain_id	residue_id	residue_name	contact_number	propensity_score	discotope_score	status'''
for index,reg in scored.iterrows():
    if y>=2:
        chain=reg["chain_id"]
        Position=reg["residue_id"]
        AminoAcid=reg["aa"]
        contact_num = reg["contact_number"]
        propensity_score = reg["propensity_score"]
        discotope_score = reg["discotope_score"]
        status = reg["status"]
        positionAnt = Position
        AminoAcidAnt = AminoAcid

        if int(Position)-int(positionAnt)==1:
            Epitope_seq=Epitope_seq+AminoAcid
            consecutius=consecutius+1
        else:
            if consecutius>=2:
                Epitope_id=Epitope_id + 1
                out_file.writelines(str(Epitope_id)+","+Epitope_seq+"\n")
            Epitope_seq = AminoAcid
            consecutius=1
    else:
        chain=reg["chain_id"]
        Position=reg["residue_id"]
        AminoAcid=reg["aa"]
        contact_num = reg["contact_number"]
        propensity_score = reg["propensity_score"]
        discotope_score = reg["discotope_score"]
        status = reg["status"]

        Epitope_seq = AminoAcid

    y=y+1
if consecutius>=2:
    Epitope_id=Epitope_id + 1
    out_file.writelines(str(Epitope_id)+","+Epitope_seq+"\n")



# OUTPUT FILE
out_file.close()
print("Output file: "+nameOutFile)
