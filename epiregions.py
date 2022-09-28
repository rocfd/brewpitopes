### EPIREGIONS

### EXTRACTS EPITOPE REGIONS

# IMPORTS
import sys
import csv
import os.path
from pathlib import Path
import pandas as pd

# HELP
h = '''
    To right usage of this script:
        $ python3 epiepregs.py
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
    outpath = sys.argv[2]

# VARIABLES
outpath=""
data_file = ""

# FILE NOT EXIST
while not fileExist(data_file):
    data_file = input(''' Provide location and name of the EPITOPE dataset.
        For example: /brewpitopes/I_final_candidates/brewpitopes_results_df.csv
            (-h for help)
        Here: ''')
print("Input file:" + data_file)

# FILE NOT EXIST
while not os.path.exists(outpath):
    outpath = input(''' Provide the desired output path.
        For example: brewpitopes/K_epitope_regions
            (-h for help)
        Here: ''')
print("Input file:" + outpath)

# OPEN READ FILE
extension = os.path.splitext(data_file)[1]
lenextension=len(extension)
nameOutFile=outpath+"/epitope_regions_extracted"+extension

### SCRIPT
## UPLOAD DATA
data = pd.read_csv(data_file, sep = ";")

## REMOVE RANK COLUMNS
data = data[data.columns.drop(list(data.filter(regex='Rank')))]

## SORT ASCENDINGLY BY START
data_sorted = data.sort_values(by ='Start', ascending = 1)
data_sorted

### epregs, SCORE AND POSITION EXTRACTOR
inici=0
iniAnt=0
fiAnt=0
iniAct=0
fiAct=0
epregs=""
epregs_vector=[]
epregs_vectors=[]
Score_ABCpred=0.0
Score_Bebipred_2_0=0.0
Score_Discotope_2_0=-10.0
#data = pd.read_csv(data_file, sep = ",")
#Rank,Sequence,Start,End,Positions,Score,Tool,Length,extracellular1,Glycosilation,accessibility
epregs_vector=["Sequence","Start","End","Score_ABCpred","Score_Bebipred_2_0","Score_Discotope_2_0"]
epregs_vectors.append(epregs_vector)
for index, row in data_sorted.iterrows():
    iniAct=row[1]
    fiAct=row[2]
    if iniAct>=iniAnt and iniAct<=fiAnt and fiAct> fiAnt:
        epregs=epregs + row[1][fiAnt-fiAct:]
    elif iniAct>iniAnt and fiAct>fiAnt:
        if epregs != "":
            epregs_vector=[]
            epregs_vector.append(epregs)
            epregs_vector.append(inici)
            epregs_vector.append(fiAnt)
            epregs_vector.append(Score_ABCpred)
            epregs_vector.append(Score_Bebipred_2_0)
            epregs_vector.append(Score_Discotope_2_0)
            epregs_vectors.append(epregs_vector)
            Score_ABCpred=0.0
            Score_Bebipred_2_0=0.0
            Score_Discotope_2_0=-10.0
        epregs= row[0][:]
        inici=iniAct
    if row[5] == "ABCpred" and row[4]>Score_ABCpred:
        Score_ABCpred=row[4]
    if row[5] == "Bebipred 2.0" and row[4]>Score_Bebipred_2_0:
        Score_Bebipred_2_0=row[4]
    if row[5] == "Discotope 2.0" and row[4]>Score_Discotope_2_0:
        Score_Discotope_2_0=row[4]
    iniAnt=iniAct
    if fiAct > fiAnt:
        fiAnt=fiAct
epregs_vector=[]
epregs_vector.append(epregs)
epregs_vector.append(inici)
epregs_vector.append(fiAnt)
epregs_vector.append(Score_ABCpred)
epregs_vector.append(Score_Bebipred_2_0)
epregs_vector.append(Score_Discotope_2_0)
epregs_vectors.append(epregs_vector)
    #for vector in epregs_vectors:
        #print(vector)

### SAVE AS DATAFRAME
epregs_df = pd.DataFrame(epregs_vectors)

### REMOVE FIRST ROW
epregs_df.columns = epregs_df.iloc[0]
epregs_df = epregs_df[1:]

## ORDER epregs BY MULTIPLE COLUMNS
epregs_df.sort_values(by = ["Score_Bebipred_2_0", "Score_Discotope_2_0", "Score_ABCpred"], ascending = [False, False, False])
epregs_df_sorted = epregs_df.sort_values(by = ["Score_Bebipred_2_0", "Score_Discotope_2_0", "Score_ABCpred"], ascending = [False, False, False])

## RANK epregs
epregs_df_sorted["epregs_rank"] = range(0,len(epregs_df_sorted["Sequence"]),1)

## ADD epregs LENGTH
epregs_df_sorted['epregs_length'] = epregs_df_sorted['Sequence'].apply(len)

## EXPORT DATA
epregs_df_sorted.to_csv(outpath+"/epitope_regions_extracted.csv", index = False)
