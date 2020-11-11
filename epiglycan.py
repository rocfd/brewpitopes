#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 18:16:13 2020
Description: detection of glycosilated residues (positions) in a protein sequences
or given an epitope within a protein.
Input: epitope sequence, parental protein sequence & glycosilated positions.
Output: epitope glycosilation: YES or NO.
@author: rocfarriolduran
"""
# imports
import sys
import csv
import os.path
from pathlib import Path

h = '''
    To right usage of this script:
        $ python3 busca3030.py
    The files in use have to be provided by stating "location/file_name.csv"
    to the input questions that appear in the console.
    <file_name> should be a .csv separated by ";"
    OPTIONAL you can provide coma separated glycans
        For example: 0,16,32,64
        If you provide the golas....????????????
    The script returns a <file_name>_out.xlsx as output.
    You can also do the following (if using Mac/Linux OS):
        $ chmod +x busca3030.py
        $ ./busca3030.py <location/file_name.csv> <comaSeparateGoals>
    You need python 3 installed in your computer!!!
    '''



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



data_file=""
inputGlycans=""
glycans=""

num_args=len(sys.argv)


if num_args>=2:
    data_file = sys.argv[1]
if num_args==3:
    inputGlycans = sys.argv[2]


while not fileExist(data_file):
    data_file = input(''' Provide location and name of the file dataset.
        For example: data/file_test.csv
            (-h for help)
        Here: ''')
print("Input file: " + data_file)

#data_file= 'test_protein_parser.csv'

extension = os.path.splitext(data_file)[1]
lenextension=len(extension)
nameOutFile=data_file[:-lenextension]+"_Out"+extension

if num_args<3:
    inputGlycans = input(''' OPTIONAL Provide coma separated glycans.
        (zero is first position)
        For example: 0,16,32,64
        Here: ''')

#glycans=[0,16,32,64]
if len(inputGlycans)>0:
    glycans=inputGlycans.split(",")

if len(glycans)!=0:
    print("glycans: " + inputGlycans)
else:
    print("No glycans provided")


def isOnGoal(start,llpep):
    isOnGoal=-1
    if len(glycans) > 0:
        for goal in glycans:
            if start<=int(goal) and int(goal)<=start+llpep-1:
                if isOnGoal==-1:
                    isOnGoal=goal
                else:
                    isOnGoal = isOnGoal + "," + goal
    return isOnGoal

out_file=open(nameOutFile, "w")
if len(glycans) > 0:
    out_file.writelines("Id;Sequence;Peptide;Position;Pre;lenPre;Post;lenPost;onGoal\n")
else:
    out_file.writelines("Id;Sequence;Peptide;Position;Pre;lenPre;Post;lenPost\n")

f=open(data_file, "r", encoding = 'utf-8-sig')
inputFile = csv.reader(f, delimiter=';')
y=0
for reg in inputFile :
    if y!=0:
        id=reg[0]
        seq=reg[1]
        pep=reg[2]
        pre=""
        post=""
        llpep=len(pep)
        start=seq.find(pep)
        if start!=-1:
            if start<30 :
                if start==0:
                    pre=""
                else:
                    pre=seq[:start-1]
            else:
                pre=seq[start-31:start-1]
            post=seq[start+llpep+1:start+llpep+31]
            if len(glycans)>0:
                onGoal=isOnGoal(start,llpep)

        if len(glycans) > 0:
            out_file.writelines(id+";"+seq+";"+pep+";"+str(start)+";"+pre+";"+str(len(pre))+";"+post+";"+str(len(post))+";"+str(onGoal)+ "\n")
        else:
            out_file.writelines(id+";"+seq+";"+pep+";"+str(start)+";"+pre+";"+str(len(pre))+";"+post+";"+str(len(post))+"\n")
    y=y+1

out_file.close()
print("Output file: "+nameOutFile)
