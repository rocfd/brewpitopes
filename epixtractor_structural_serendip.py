#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 5 14:21:51 2024
@autor: victormontal
"""
### EPITOPE EXTRACTOR
## GOAL: Obtain conformational epitopes (3D groups of residues )
## from Serendip-CE

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
from biopandas.pdb import PandasPdb
import argparse
import ipdb

# ---------------------------
# Parser Options
# ---------------------------
HELPTEXT = f"""

epixtractor_structural_serendipce.py [dev version]

Parse the output from (conformational) epitope prediction from Serendip-CE

Steps:
- input the brewpitope project folder
- load pdb and csv results from Serendip-ce
- generate binary/mask file (using csv file)
- surface-based cluster based on surface neighbours
- generate outputs [.ply surface, .pdb, .csv]

Example:
--------
python3 epixtractor_structural_serendipce.py --path example/path/project/brewpitope

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
def parse_serendip(args):
    # Defaults
    ifile = join(args.ipath,"B_structural_predictions","serendipce","serendipce_predict.csv")
    if not os.path.exists(ifile):
        print("File does not exist:", ifile)
        print(f"Remember that Serendip-ce output MUST be named serendipce_predict.csv ")
        print(f"  and placed at folder 'B_structural_predictions/serendipce/'" )
        sys.exit(1)

    inpdb = join(args.ipath,"B_structural_predictions","serendipce","serendipce_renumpdb.pdb")
    if not os.path.exists(inpdb):
        print("File does not exist:", inpdb)
        print(f"Remember that Serendip-ce output MUST be named serendipce_renumpdb.pdb ")
        print(f"  and placed at folder 'B_structural_predictions/serendipce/'" )
        sys.exit(1)

    cwd = os.getcwd()
    clustf = join(cwd,"prot_surface_cluster.py")
    if not os.path.exists(clustf):
        print(f"File {clustf} NOT FOUND in current directory")
        print(f" It must be placed at the same directory where current script is run")
        print(f"     usually within brewpitopes repo.")
        sys.exit(1)

    # Load fasta output from Bepipred3.0
    print("> Loading csv and pdb file of SerenDIP-CE results")
    serendf = pd.read_csv(ifile, sep="\t")

    ppdb = PandasPdb().read_pdb(inpdb)
    atoms = ppdb.df["ATOM"]
    all_resid = atoms["residue_number"].values


    print("> Compute features mask")
    resid = np.unique(all_resid)
    mask = np.zeros([resid[-1]])

    for row in serendf.iterrows():
        resid = row[1]["SeqPos"]
        epitope = row[1]["prediction"]
        if epitope == "I":
            mask[resid] = 1
    omask = join(args.ipath, "B_structural_predictions","serendipce", "features_aa_serendipce.txt")
    np.savetxt(omask,mask,fmt='%.1f')


    print("> Generate groups using a surface-based approach")
    oname = "serendipce"
    opath = join(args.ipath, "B_structural_predictions","serendipce")
    cmd = ""
    cmd += "python prot_surface_cluster.py "
    cmd += f"--pdb {inpdb} --feature {omask} --outpath {opath} --outname {oname}"
    run_cmd(cmd,"Can not compute surface-based grouping of epitopes")


    print("> Dump results")
    cporigin = join(opath,f"groups.{oname}.csv")
    df_clen = pd.read_csv(cporigin)
    df_clen['Tool'] = "Serendip-CE"
    df_clen.to_csv(cporigin, index=False)
    cpdest = join(args.ipath, "C_epixtractor", "serendipce_results_extracted.csv")
    copyfile(cporigin, cpdest)




#  -Run Code -
# ------------
if __name__ == "__main__":
    args = options_parser()
    parse_serendip(args)
