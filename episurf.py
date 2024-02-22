#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on XX XX XX 14:21:51 2022
@autor: rocfarriolduran

Modified on Tue Feb 5 14:21:51 2024
@autor: victormontal
"""
### EPITOPE Surface label
## GOAL: Identify if (lineal) epitopes are buried or no

# NOTES:
# Conformational epitopes will not be evaluated since we have already filter-out
# those residues that are NOT on the surface

# We will cluster the residues based on surface-based neighbours



# ---------------------------
# Import libraries
# ---------------------------
import os
import sys
from os.path import join, dirname, basename
from shutil import copyfile
import subprocess
import pandas as pd
import numpy as np
import mdtraj as md
import argparse
import ipdb

# ---------------------------
# Parser Options
# ---------------------------
HELPTEXT = f"""

episurf.py [dev version]

Parse the output from epiptm.R

Steps:
- input the brewpitope project folder
- load epitope table
- compute SASA and filter residues > XX
- iterate each residue and mark those with >80% accesible
- generate outputs [.csv]

Example:
--------
python3 episurf.pyy --path example/path/project/brewpitope

Author:
------
Victor Montal
victor.montal [at] protonmail [dot] com


"""

USAGE = r"""

"""

def options_parser():
    """
    Command Line Options Parser:
    initiate the option parser and return the parsed object
    """
    parser = argparse.ArgumentParser(description = HELPTEXT,usage = HELPTEXT)

    # help text
    h_projectpath = 'Path to project folder'

    # Parser
    parser.add_argument('--path',
                        dest = 'ipath', action = 'store',
                        help = h_projectpath, required = True)

    args = parser.parse_args()
    return args

# ---------------------------
# Extra functions
# ---------------------------
def run_cmd(cmd,err_msg):
    """
    execute the comand
    """
    print('#@# Command: ' + cmd+'\n')
    retcode = subprocess.Popen(cmd,shell=True, executable='/bin/bash').wait()
    if retcode != 0 :
        print('ERROR: '+err_msg)
        sys.exit(1)
    print('\n')

def overlap(set1, set2):
    """
    Compute the overlap of set1 within set2
    """
    inter = np.intersect1d(set1, set2)
    overlap = inter.size / len(set1)
    return overlap

# ---------------------------
# Main Code
# ---------------------------
def episurf(args):
    # Defaults
    ifile = join(args.ipath,"F_epiptm","ptm_extracted.csv")
    if not os.path.exists(ifile):
        print("File does not exist:", ifile)
        print(f"Remember to run epiptm.R before running episurf.py (!!) ")
        sys.exit(1)

    inpdb = join(args.ipath,"G_episurf","episurf.pdb")
    if not os.path.exists(inpdb):
        print("File does not exist:", inpdb)
        print(f"Remember that PDB output MUST be named episurf.pdb ")
        print(f"  and placed at folder 'G_episurf/'" )
        sys.exit(1)

    print("> Load .csv (after PTM labeling)")
    epi_df = pd.read_csv(ifile, sep=";")
    oepi_df = epi_df.copy()
    oepi_df["accessibility"] = "Buried"

    # Define MaxASA
    maxasa_dic = {
                    "ALA" : 129,
                    "ARG" : 274,
                    "ASN" : 195,
                    "ASP" : 193,
                    "CYS" : 167,
                    "GLU" : 223,
                    "GLN" : 225,
                    "GLY" : 104,
                    "HIS" : 224,
                    "ILE" : 197,
                    "LEU" : 201,
                    "LYS" : 236,
                    "MET" : 224,
                    "PHE" : 240,
                    "PRO" : 159,
                    "SER" : 155,
                    "THR" : 172,
                    "TRP" : 285,
                    "TYR" : 263,
                    "VAL" : 174
    }

    # Compute SASA
    print("> Compute residue-level SASA (and filter residues)")
    itraj = md.load(inpdb)
    resid = [xx.resSeq for xx in itraj.topology.residues]
    resname = [xx.name for xx in itraj.topology.residues]
    resid_sasa = md.shrake_rupley(itraj, mode='residue').flatten()
    resid_maxasa = [ maxasa_dic[xx] for xx in resname]
    sasa_pct = (resid_sasa*100) / np.array(resid_maxasa)   # scale to Amst
    pos_sasa = np.where(sasa_pct > 0.2)[0]
    resid_sasa = [resid[xx] for xx in pos_sasa]

    print("> Iterate over all epitope")
    for idx,row in epi_df.iterrows():
        cpos = row["Positions"]
        cpos = [int(xx) for xx in cpos.split(",")]
        coverlap = overlap(cpos,resid_sasa)

        if coverlap > 0.8:
            oepi_df.loc[idx,"accessibility"] = "Accessible"


    print("> Dump results")
    ocsv = join(args.ipath,"G_episurf","access_extracted.csv")
    oepi_df.to_csv(ocsv, index=False)




#  -Run Code -
# ------------
if __name__ == "__main__":
    args = options_parser()
    episurf(args)
