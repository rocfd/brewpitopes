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
    ifile = join(args.ipath,"A_linear_predictions","bepipred","bepipred_linear_pred.fasta")
    if not os.path.exists(ifile):
        print("File does not exist:", ifile)
        print(f"Remember that Bepitope3 linear output MUST be named bepipred_results.fasta ")
        print(f"  and placed at folder 'A_linear_predictions/bepipred/'" )
        sys.exit(1)

    ofile = join(args.ipath, "C_epixtractor", "bepipred_results_extracted.csv")

    # Load fasta output from Bepipred3.0
    print("> Loading fasta file of Bepipred3.0 results")
    with open(ifile, 'r') as f:
        for line in f:
            if line.startswith(">"):
                continue
            else:
                sequence = line.strip()

    # Define output pandas dict
    odict = {}
    odict["Rank"] = []
    odict["Sequence"] = []
    odict["Start"] = []
    odict["End"] = []
    odict["Positions"] = []
    odict["BebiScore"] = []

    # Loop over fasta seq and identify linear peptides based on capital letters
    print("> Identify linear epitopes")
    start = 0
    stop = 0
    idx1 = 0
    while idx1 < len(sequence):
        cstring = sequence[idx1]
        cepitope = []
        if cstring.isupper():
            cepitope = [cstring]
            cepitopepos = [idx1]
            start = idx1

            for idx2 in range(idx1+1,len(sequence)):
                evalstring = sequence[idx2]
                if evalstring.isupper():
                    cepitope.append(evalstring)
                    cepitopepos.append(idx2)
                else:
                    idx1 = idx2
                    stop = idx2-1

                    if len(cepitope)> 2:
                        odict["Rank"].append(0)
                        odict["Sequence"].append("".join(cepitope))
                        odict["Start"].append(start)
                        odict["End"].append(stop)
                        odict["Positions"].append(cepitopepos)
                        odict["BebiScore"].append("NA")
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
