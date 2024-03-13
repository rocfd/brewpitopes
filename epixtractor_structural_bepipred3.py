#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 5 14:21:51 2024
@autor: victormontal
"""
### EPITOPE EXTRACTOR
## GOAL: Obtain conformational epitopes (3D groups of residues )
## from bepipred33.0

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

epixtractor_structural_bepipred33.py [dev version]

Parse the output from (conformational) epitope prediction from bepipred33.0

Steps:
- input the brewpitope project folder
- load pdb results from bepipred33.0
- generate binary/mask file (using fasta upper letters)
- surface-based cluster based on surface neighbours
- generate outputs [.ply surface, .pdb, .csv]

Example:
--------
python3 epixtractor_structural_bepipred33.py --path example/path/project/brewpitope

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
def parse_bepipred3(args):
    # Defaults
    ifile = join(args.ipath,"B_structural_predictions","bepipred3_conf","bepipred3_conf_output.csv")
    if not os.path.exists(ifile):
        print("File does not exist:", ifile)
        print(f"Remember that bepipred3 output MUST be named bepipred3_conf_output.csv")
        print(f"  and placed at folder 'B_structural_predictions/bepipred3_conf/'" )
        sys.exit(1)

    inpdb = join(args.ipath,"B_structural_predictions","bepipred3_conf","bepipred3_renumpdb.pdb")
    if not os.path.exists(inpdb):
        print("File does not exist:", inpdb)
        print(f"Remember that protein PDB MUST be named bepipred3_renumpdb.pdb ")
        print(f"  and placed at folder 'B_structural_predictions/bepipred3_conf/'" )
        sys.exit(1)

    cwd = os.getcwd()
    clustf = join(cwd,"prot_surface_cluster.py")
    if not os.path.exists(clustf):
        print(f"File {clustf} NOT FOUND in current directory")
        print(f" It must be placed at the same directory where current script is run")
        print(f"     usually within brewpitopes repo.")
        sys.exit(1)

    """
    # Load fasta output from bepipred33.0
    print("> Loading fasta file of bepipred3.0 results")
    with open(ifile, 'r') as f:
        for line in f:
            if line.startswith(">"):
                continue
            else:
                sequence = line.strip()
    """
    print("> Compute features mask")
    epi_csv = pd.read_csv(ifile)
    score = epi_csv["BepiPred-3.0 score"]
    pos_thr = np.where(score >= 0.1512)[0]
    mask = np.zeros([score.size])
    mask[pos_thr] += 1
    omask = join(args.ipath, "B_structural_predictions","bepipred3_conf", "features_aa_bepipred3.txt")
    np.savetxt(omask,mask,fmt='%.1f')

    print("> Generate groups using a surface-based approach")
    oname = "bepipred3"
    opath = join(args.ipath, "B_structural_predictions","bepipred3_conf")
    cmd = ""
    cmd += "python prot_surface_cluster.py "
    cmd += f"--pdb {inpdb} --feature {omask} --outpath {opath} --outname {oname}"
    run_cmd(cmd,"Can not compute surface-based grouping of epitopes")

    print("> Dump results")
    cporigin = join(opath,f"groups.{oname}.csv")
    df_clen = pd.read_csv(cporigin)
    df_clen['Tool'] = "Bepipred3.0_conf"
    # Compute average score per group of epitope
    for i,row in df_clen.iterrows():
        cpos = row["Positions"]
        cpos = [int(xx) -1 for xx in cpos.split(",")]   # 0 vs 1 encoded
        cscore = np.average(score[cpos])
        df_clen.loc[i,"Score"] = cscore
        df_clen.loc[i,"Start"] = cpos[0] +1
        df_clen.loc[i,"End"] = cpos[-1]
    df_clen.to_csv(cporigin, index=False)
    cpdest = join(args.ipath, "C_epixtractor", "bepipred3_conf_results_extracted.csv")
    copyfile(cporigin, cpdest)




#  -Run Code -
# ------------
if __name__ == "__main__":
    args = options_parser()
    parse_bepipred3(args)
