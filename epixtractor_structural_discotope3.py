#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 5 14:21:51 2024
@autor: victormontal
"""
### EPITOPE EXTRACTOR
## GOAL: Obtain conformational epitopes (3D groups of residues )
## from Discope3.0

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
import argparse
import ipdb

# ---------------------------
# Parser Options
# ---------------------------
HELPTEXT = f"""

epixtractor_structural_discotope3.py [dev version]

Parse the output from (conformational) epitope prediction from Discotope3.0

Steps:
- input the brewpitope project folder
- load pdb results from Discotope3.0
- generate binary/mask file (i.e residues > threshold)
- surface-based cluster based on surface neighbours
- generate outputs [.ply surface, .pdb, .csv]

Example:
--------
python3 epixtractor_structural_discotope3.py --path example/path/project/brewpitope

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

# ---------------------------
# Main Code
# ---------------------------
def parse_discotope(args):
    # Defaults
    ifile = join(args.ipath,"B_structural_predictions","discotope","discotope_pred.csv")
    if not os.path.exists(ifile):
        print("File does not exist:", ifile)
        print(f"Remember that Discotope output MUST be named discotope_pred.csv ")
        print(f"  and placed at folder 'B_structural_predictions/discotope/'" )
        sys.exit(1)

    inpdb = join(args.ipath,"B_structural_predictions","discotope","discotope_pred.pdb")
    if not os.path.exists(inpdb):
        print("File does not exist:", inpdb)
        print(f"Remember that Discotope output MUST be named discotope_pred.pdb ")
        print(f"  and placed at folder 'B_structural_predictions/discotope/'" )
        sys.exit(1)

    cwd = os.getcwd()
    clustf = join(cwd,"prot_surface_cluster.py")
    if not os.path.exists(clustf):
        print(f"File {clustf} NOT FOUND in current directory")
        print(f" It must be placed at the same directory where current script is run")
        print(f"     usually within brewpitopes repo.")
        sys.exit(1)

    # Load .csv and get epitope residues
    print("> Load .csv output from Discotope3.0")
    disco_df = pd.read_csv(ifile)
    last_residue = np.unique(disco_df["res_id"])[-1]
    group = np.zeros([last_residue])

    print("> Compute features mask")
    for cresidue in np.unique(disco_df["res_id"]):
        is_epi = disco_df[disco_df["res_id"] == cresidue]["epitope"]
        if is_epi.any():
            group[cresidue-1] = 1  #zero-coded
    omask = join(args.ipath, "B_structural_predictions","discotope", "features_aa_discotope.txt")
    np.savetxt(omask,group,fmt='%.1f')


    print("> Generate groups using a surface-based approach")
    oname = "discotope3"
    opath = join(args.ipath, "B_structural_predictions","discotope")
    cmd = ""
    cmd += "python prot_surface_cluster.py "
    cmd += f"--pdb {inpdb} --feature {omask} --outpath {opath} --outname {oname}"
    run_cmd(cmd,"Can not compute surface-based grouping of epitopes")


    print("> Dump results")
    cporigin = join(opath,"groups.discotope3.csv")
    df_clen = pd.read_csv(cporigin)
    # Compute average score per group of epitope
    for i,row in df_clen.iterrows():
        cpos = row["Positions"]
        cpos = [int(xx) for xx in cpos.split(",")]   # 0 vs 1 encoded
        residue_scores = [disco_df[disco_df["res_id"] == xx]["DiscoTope-3.0_score"].values for xx in cpos]
        cscore = np.average(residue_scores)
        df_clen.loc[i,"Score"] = cscore
        df_clen.loc[i,"Start"] = cpos[0]
        df_clen.loc[i,"End"] = cpos[-1]
    df_clen['Tool'] = "Discotope3"
    df_clen.to_csv(cporigin, index=False)
    cpdest = join(args.ipath, "C_epixtractor", "discotope_results_extracted.csv")
    copyfile(cporigin, cpdest)




#  -Run Code -
# ------------
if __name__ == "__main__":
    args = options_parser()
    parse_discotope(args)
