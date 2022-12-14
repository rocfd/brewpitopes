#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 14:21:51 2020
@author: rocfarriolduran
"""
### EPITOPE EXTRACTOR
## GOAL: extract linear epitopes in tabulated data either from linear or structural prediction origin.
## Variant for linear epitope predictions at: 
## https://services.healthtech.dtu.dk/service.php?BepiPred-2.0
## SUBGOAL: Create epitopes from continous residues that surpass a given quality threshold in Linear Predictions predictions.


# IMPORTS
import sys
import csv
import os.path
from pathlib import Path
import pandas as pd
import more_itertools as mit

## HELP
h = '''
    To right usage of this script:
        $ python3 epixtractor_linear_bepipred.py
        The files in use have to be provided by stating
        "location/file_name.csv"
    to the input questions that appear in the console.
    <file_name> should be a .csv separated by ","
    OPTIONAL 
    The script returns a <file_name>_out.csv as output.
    You can also do the following (if using Mac/Linux OS):
        $ chmod +x epixtractor_linear.py 
        $ ./epixtractor_linear_bepipred.py  <location/file_name.csv>
    You need python 3 installed in your computer!!!
    '''


## FILE EXISTANCE
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


## FILE UPLOAD
num_args=len(sys.argv)

if num_args>=2:
    data_file = sys.argv[1]
if num_args==3:
    outpath = sys.argv[2]

## VARIABLES
data_file=""
outpath=""

while not fileExist(data_file):
    data_file = input(''' Provide location and name of the file dataset.
        For example: data/file_test.csv
            (-h for help)
        Here: ''')

print("Input file: " + data_file)

# OUTPUT PATH
while not os.path.exists(outpath):
    outpath = input(''' Provide the desired output path.
        For example: /brewpitopes/C_epixtractor
            (-h for help)
        Here: ''')
print("Input file:" + outpath)


# NAME THE OUTPUT FILE
extension = os.path.splitext(data_file)[1]
lenextension=len(extension)
nameOutFile=outpath+"/bepipred_results_extracted"+extension


## PIPELINE for EXTRACTION OF EPITOPES FROM LINEAR PREDICTIONS

# READ THE CSV INPUT & ADD HEADERS TO BEBIPRED DATAFRAME
prediction = pd.read_csv(data_file, 
			 sep=',', 
			 header = "infer")
    
# Filter by bebipred_score threshold at 0.55
scored = prediction[prediction.EpitopeProbability >= 0.55]

# Reset index of filtered DF to use indexes later on
scored = scored.reset_index(drop=True)

# Loop to append a position if the previous is contiguous
resid = []
for i in range(len(scored)):
    try:
        if scored.loc[i,'Position']-scored.loc[i-1,'Position'] == 1:
            resid.append(scored.loc[i,'Position'])
    
    ## misses first aminoacid of a stretch.
    except:
        resid.append(scored.loc[i,'Position'])
        continue
        
# Group the epitopes and extract positions
resid_grouped = [list(group) for group in mit.consecutive_groups(resid)]

# filter for groups larger than 2 elements
resid_grouped_filtered = [group for group in resid_grouped if len(group)>2]
resid_grouped_filtered

# EXTRACT THE CONTINOUS SEQUENCES
sequences = []
for group in resid_grouped_filtered:
    sequence = []
    for n in group:
        sequence.append(scored.AminoAcid[scored.Position==n].values[0])
    sequences.append(''.join(sequence))
sequences

# LOOP TO EXTRACT START AND END POSITION
start = []
end = []
for group in resid_grouped_filtered:
    start.append(group[0])
    end.append(group[-1])
    #print(start)
    #print(end)

# EXTRACT BEBIPRED SCORE
BebiScore = scored.EpitopeProbability
    
# CREATE LIST TO PREPARE OUTPUT DATAFRAME
out = list(zip(sequences,start,end,resid_grouped_filtered,BebiScore))
#print (out)

# CREATE OUTPUT DATAFRAME
out_file = pd.DataFrame(out, columns = ("Sequence", "Start", "End", "Positions", "BebiScore"))
#print(out_file)

# EXPORT OUTPUT FILE
out_file.to_csv(nameOutFile, sep = ";", index = True, index_label = "Rank")
print("Output file: "+nameOutFile)

if out_file.empty:
    print("Bepipred 2.0 could not predict any epitopes in your target sequence. You will get an empty dataframe and you can continue the pipeline with the other predictors.")
else:
    print("Bepipred 2.0 could predict one or more epitopes in your target sequence. Go ahead!")
