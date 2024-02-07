#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 5 14:21:51 2024
@autor: victormontal
"""
### EPITOPE EXTRACTOR
## GOAL: Obtain conformational epitopes (3D groups of residues )
## from SEPPAd3.0

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
from biopandas.pdb import PandasPdb
import numpy as np
import argparse
import ipdb

# ---------------------------
# Parser Options
# ---------------------------
HELPTEXT = f"""

epixtractor_structural_seppa3.py [dev version]

Parse the output from (conformational) epitope prediction fromSeppa3.0

Steps:
- input the brewpitope project folder
- load pdb results from SEPPA3.0
- generate binary/mask file (using prediction from seppa)
- surface-based cluster based on surface neighbours
- generate outputs [.ply surface, .pdb, .csv]

Example:
--------
python3 epixtractor_structural_seppa3.py --path example/path/project/brewpitope

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
def parse_seppa(args):
    # Defaults
    iseppa = join(args.ipath,"B_structural_predictions","seppa","seppa_predict.pdb")
    if not os.path.exists(iseppa):
        print("File does not exist:", iseppa)
        print(f"Remember that Seppa output MUST be named seppa_predict.pdb ")
        print(f"  and placed at folder 'B_structural_predictions/seppa/'" )
        sys.exit(1)

    cwd = os.getcwd()
    clustf = join(cwd,"prot_surface_cluster.py")
    if not os.path.exists(clustf):
        print(f"File {clustf} NOT FOUND in current directory")
        print(f" It must be placed at the same directory where current script is run")
        print(f"     usually within brewpitopes repo.")
        sys.exit(1)

    # Load output from Seppa3.0
    print("> Loading PDB file of Seppa3.0 results")
    # Get residue numbers
    ppdb = PandasPdb().read_pdb(iseppa)
    atoms = ppdb.df["ATOM"]
    all_resid = atoms["residue_number"].values

    # Get residue binary value per atom
    with open(iseppa, 'r') as f:
        isepi = []
        for line in f:
            line = line.strip()
            epitmp = line.split("\t")[1]  # get bin epitope or not
            if epitmp != "NonEpitope":
                isepi.append(1)
            else:
                isepi.append(0)

    print("> Compute features mask")
    resid = np.unique(all_resid)
    mask = np.zeros([resid[-1]])
    for cresidue in resid:
        pos = np.where(all_resid == cresidue)[0]
        value = isepi[pos[0]]
        mask[cresidue-1] = value   #zero-coded

    omask = join(args.ipath, "B_structural_predictions","seppa", "features_aa_seppa.txt")
    np.savetxt(omask,mask,fmt='%.1f')

    print("> Generate groups using a surface-based approach")
    oname = "seppa"
    opath = join(args.ipath, "B_structural_predictions","seppa")
    cmd = ""
    cmd += "python prot_surface_cluster.py "
    cmd += f"--pdb {iseppa} --feature {omask} --outpath {opath} --outname {oname}"
    run_cmd(cmd,"Can not compute surface-based grouping of epitopes")


    print("> Dump results")
    cporigin = join(opath,f"groups.{oname}.csv")
    df_clen = pd.read_csv(cporigin)
    df_clen['Tool'] = "Seppa3"
    df_clen.to_csv(cporigin, index=False)
    cpdest = join(args.ipath, "C_epixtractor", "seppa_results_extracted.csv")
    copyfile(cporigin, cpdest)




#  -Run Code -
# ------------
if __name__ == "__main__":
    args = options_parser()
    parse_seppa(args)
