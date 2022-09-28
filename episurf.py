## EPISURF
## GOAL: LABEL THE EPITOPE BASED ON ACCESSIBILITY ON THE PARENTAL PROTEIN

# IMPORTS
import sys
import csv
import os.path
from pathlib import Path
import pandas as pd

# HELP
h = '''
    To right usage of this script:
        $ python3 episurf.py
    The files in use have to be provided by stating "location/file_name.csv"
    to the input questions that appear in the console.
    <file_name> should be a .csv separated by ";"
    You have to provide coma separated buried positions
        For example: 0,16,32,64
    The script returns a <file_name>_out.xlsx as output.
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
    inputaccess = sys.argv[2]
if num_args==4:
    outpath = sys.argv[3]

# VARIABLES
data_file=""
inputaccess=""
access=""
outpath=""

# FILE NOT EXIST
while not fileExist(data_file):
    data_file = input(''' Provide location and name of the file dataset.
        For example: brewpitopes/F_epiglycan/glycan_extracted.csv
            (-h for help)
        Here: ''')
print("Input file:" + data_file)

# FILE NOT EXIST
while not os.path.exists(outpath):
    outpath = input(''' Provide the desired output path.
        For example: brewpitopes/G_episurf
            (-h for help)
        Here: ''')
print("Input file:" + outpath)

# OPEN READ FILE
extension = os.path.splitext(data_file)[1]
lenextension=len(extension)
nameOutFile=outpath+"/access_extracted"+extension

# GLYCAN INPUT
while not fileExist(inputaccess):
    inputaccess = input(''' Provide location and name of the file dataset.
        For example: brewpitopes/G_episurft/buried_positions_list.csv
            (-h for help)
        Here: ''')
print("Input file:" + inputaccess)

## READ DATA
data = pd.read_csv(data_file, sep = ",")
data


## READ MOLSOFT EXTRACTED RESULTS (single column file titled "buried" and containing buried positions separated by rows")
access_df = pd.read_csv(inputaccess, sep = ",")
access_df

## DF TO LIST
access = access_df["buried"]
access

## LOOP TO EXTRACT access
z = []
for index, row in data.iterrows():
    for ac in access:
        y = "Accessible"
        if int(ac) >= row["Start"] and int(ac) <= row["Start"] + row["Length"] - 1:
            y = "Buried"
            z.append(y)
            break
    if y == "Accessible":
        z.append(y)
print(z)
len(z)

## APPEND GLYCOSILATION TO DF
data['accessibility_icm'] = z
data

## EXPORT DATA
data.to_csv(path_or_buf= nameOutFile,
         index = True, index_label = "Rank")
