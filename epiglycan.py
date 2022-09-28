# IMPORTS
import sys
import csv
import os.path
from pathlib import Path
import pandas as pd

# HELP
h = '''
    To right usage of this script:
        $ python3 epiglycan.py
    The files in use have to be provided by stating "location/file_name.csv"
    to the input questions that appear in the console.
    <file_name> should be a .csv separated by ";"
    Additionally you have to provide coma separated glycans
        For example: 0,16,32,64
    The script returns a <file_name>_out.csv as output.
    You need python 3 installed in your computer!!!
    '''


# FILE CHECK
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

num_args=len(sys.argv)


if num_args>=2:
    data_file = sys.argv[1]
if num_args==3:
    inputGlycans = sys.argv[2]
if num_args==4:
    outpath = sys.argv[3]

# VARIABLES
data_file=""
inputGlycans=""
glycans=""
outpath=""

# FILE NOT EXIST
while not fileExist(data_file):
    data_file = input(''' Provide location and name of the EPITOPE dataset.
        For example: data/file_test.csv
            (-h for help)
        Here: ''')
print("Input file:" + data_file)

# FILE NOT EXIST
while not os.path.exists(outpath):
    outpath = input(''' Provide the desired output path.
        For example: brewpitopes/F_epiglycan
            (-h for help)
        Here: ''')
print("Input file:" + outpath)


# OPEN READ FILE
extension = os.path.splitext(data_file)[1]
lenextension=len(extension)
nameOutFile=outpath+"/glycan_extracted"+extension

# GLYCAN INPUT
while not fileExist(inputGlycans):
    inputGlycans = input(''' Provide location and name of the file dataset.
        For example: data/file_test.csv
            (-h for help)
        Here: ''')
print("Input file:" + inputGlycans)

## READ DATA
data = pd.read_csv(data_file, sep = ";")
data

## READ GLYCANS (single column file with "glyc_pos" as column title and positions in separate rows)
glycans_df = pd.read_csv(inputGlycans, sep = ",")
glycans_df

## DF TO LIST
glycans_list = glycans_df["glyc_pos"].values.tolist()
glycans_list
type(glycans_list)

## LOOP TO DETERMINE GLYCOSILATED EPITOPES
z = []
for index, row in data.iterrows():
    for glyc in glycans_list:
        y = "Non-glycosilated"
        if glyc >= row["Start"] and int(glyc) <= row["Start"] + row["Length"] - 1:
            y = "Glycosilated"
            z.append(y)
            break
    if y == "Non-glycosilated":
        z.append(y)
print(z)
len(z)

## APPEND GLYCOSILATION TO DF
data['Glycosilation'] = z
data

## EXPORT DATA
data.to_csv(path_or_buf= nameOutFile,
         index = True, index_label = "Rank")
