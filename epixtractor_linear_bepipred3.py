#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 14:21:51 2020
@author: rocfarriolduran

Modified on Jan 31 2024
@autor: victormontal
"""
### EPITOPE EXTRACTOR
## GOAL: extract linear epitopes in tabulated data either from linear prediction origin.
## Variant for linear epitope predictions at:
## https://services.healthtech.dtu.dk/services/BepiPred-3.0/


# ---------------------------
# Import libraries
# ---------------------------
import os
import sys
from os.path import join, dirname, basename
import pandas as pd
import numpy as np
import argparse
import ipdb

# ---------------------------
# Parser Options
# ---------------------------
HELPTEXT = f"""

epixtractor_liner_bepipred3.py [dev version]

Parse the output from (linear) epitope prediction from Bepipred3.0

Steps:
- input the brewpitope project folder
- load fasta results from Bepipred3.0 server ()
      [use the default parameters: high-confindence, threshold:0.1512]
- output the linear results to ProjectPath/C_epixtractor

Example:
--------
python3 epixtractor_liner_bepipred3.py --path example/path/project/brewpitope

Author:
------
Victor Montal
victor.montal [at] protonmail [dot] com

Roc Farriol Duran

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
# Main Code
# ---------------------------
def parse_bepipred(args):
    # Dfaults
    ifile = join(args.ipath,"A_linear_predictions","bepipred","bepipred3_lin_output.csv")
    if not os.path.exists(ifile):
        print("File does not exist:", ifile)
        print(f"Remember that bepipred3 output MUST be named bepipred3_lin_output.csv")
        print(f"  and placed at folder 'A_linear_predictions/bepipred/'" )
        sys.exit(1)

    ofile = join(args.ipath, "C_epixtractor", "bepipred_results_extracted.csv")

    """
    # Load fasta output from Bepipred3.0
    print("> Loading fasta file of Bepipred3.0 results")
    with open(ifile, 'r') as f:
        for line in f:
            if line.startswith(">"):
                continue
            else:
                sequence = line.strip()
    """
    # Read input data
    epi_csv = pd.read_csv(ifile)
    size_df = epi_csv.shape[0]
    print(size_df)
    scores = epi_csv["BepiPred-3.0 linear epitope score"].values

    # Define threshold
    thr = 0.1512

    # Define output pandas dict
    odict = {}
    odict["Rank"] = []
    odict["Sequence"] = []
    odict["Start"] = []
    odict["End"] = []
    odict["Positions"] = []
    odict["Score"] = []
    odict["Length"] = []
    odict["Tool"] = "Bepipred3.0_lin"

    # Loop over fasta seq and identify linear peptides based on capital letters
    print("> Identify linear epitopes")
    start = 0
    stop = 0
    idx1 = 0
    while idx1 < size_df:
        cscore = epi_csv.loc[idx1,"BepiPred-3.0 linear epitope score"]
        cepitope = []
        cepitopepos = []

        if cscore >= thr:
            cepitope.append(epi_csv.loc[idx1,"Residue"])
            cepitopepos.append(idx1)
            start = idx1

            for idx2 in range(idx1+1,size_df):
                evalscore = epi_csv.loc[idx2,"BepiPred-3.0 linear epitope score"]
                if evalscore >= thr:
                    cepitope.append(epi_csv.loc[idx2,"Residue"])
                    cepitopepos.append(idx2)
                else:
                    idx1 = idx2
                    stop = idx2-1

                    if len(cepitope)> 2:
                        posi_str = [str(xx) for xx in cepitopepos]

                        odict["Rank"].append(0)
                        odict["Sequence"].append("".join(cepitope))
                        odict["Start"].append(start)
                        odict["End"].append(stop)
                        odict["Positions"].append(",".join(posi_str))
                        odict["Score"].append(np.average(scores[cepitopepos]))
                        odict["Length"].append(len(cepitope))
                    break
        else:
            idx1 += 1

    # Df from dict
    outdf = pd.DataFrame(odict)

    # Output to csv
    print(f"> Dump to csv at: {ofile}")
    outdf.to_csv(ofile, index = False)


#  -Run Code -
# ------------
if __name__ == "__main__":
    args = options_parser()
    parse_bepipred(args)
