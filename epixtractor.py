#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 14:21:51 2020

@author: rocfarriolduran
"""
### EPITOPE EXTRACTOR
## GOAL: extract linear epitopes in tabulated data either from linear or structural prediction origin.
## create epitopes from continous residues that surpass a given quality threshold.


# IMPORTS
import sys
import csv
import os.path
from pathlib import Path

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

extension = os.path.splitext(data_file)[1]
lenextension=len(extension)
nameOutFile=data_file[:-lenextension]+"_Out"+extension
out_file=open(nameOutFile, "w")
out_file.writelines("Epitope_id,Epitope_seq\n")

f=open(data_file, "r", encoding = 'utf-8-sig')
inputFile = csv.reader(f, delimiter=',')

# EPITOPE EXTRACTOR
y = 0
positionAnt = 0
AminoAcidAnt = ""
Epitope_seq = ""
Epitope_id=0
consecutius=1
for reg in inputFile :
    if y>=2:
        positionAnt = Position
        AminoAcidAnt = AminoAcid

        ID=reg[0]
        Entry=reg[1]
        Position=reg[2]
        AminoAcid=reg[3]

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
        ID=reg[0]
        Entry=reg[1]
        Position=reg[2]
        AminoAcid=reg[3]

        Epitope_seq = AminoAcid

    y=y+1
if consecutius>=2:
    Epitope_id=Epitope_id + 1
    out_file.writelines(str(Epitope_id)+","+Epitope_seq+"\n")

# OUTPUT FILE
out_file.close()
print("Output file: "+nameOutFile)
