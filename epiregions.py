#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2021
@autor: racfarriol

Modified on Tue March 12 14:21:51 2024
@autor: victormontal
"""

### Merge epitopes based on their
#


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

epiregions.py [dev version]

Merge epitopes based on their overlapping pattern.


Example:
--------
python3 epiregions.pyy --path example/path/project/brewpitope

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
# Main Code
# ---------------------------
def epiregions(args):
    # Defaults
    ifile = join(args.ipath,"I_final_candidates","brewpitopes_results_df.csv")
    if not os.path.exists(ifile):
        print("File does not exist:", ifile)
        sys.exit(1)
    # Load fasta file
    ifasta = join(args.ipath,"Z_fasta","protein.fasta")
    if not os.path.exists(ifasta):
        print("File does not exist:", ifile)
        sys.exit(1)
    print("> Loading protein fasta file")
    with open(ifasta, 'r') as f:
        for line in f:
            if line.startswith(">"):
                continue
            else:
                sequence = line.strip()


    opath = join(args.ipath,"K_epitope_regions")
    ofile = join(opath, "epitope_regions_extracted.csv")

    # Load dataframe
    data = pd.read_csv(ifile, sep = ";")

    #
    # Identify overlapping epitopes
    #

    # Init dicts
    tools = ["Bepipred3.0_lin", "EpitopeVec", "Seppa3","Discotope3",
             "Bepipred3.0_conf","Serendip-CE"]
    epi_dict = {}
    keys = ["Sequence","Positions", "Length"]
    for ctool in tools:
        keys.append(f"{ctool}")
    for ckey in keys:
        epi_dict[ckey] = []

    # Loop all epitopes
    nepitopes = 0
    for i,row in data.iterrows():
        # Get current epitope info
        cpositions = row["Positions"]
        cpositions = [int(xx) for xx in cpositions.split(",")]
        ctool = row["Tool"]
        cscore = float(row["Score"])

        # Search for overlap
        found_ov = False
        for idx,fetched in enumerate(epi_dict["Positions"]):
            overlap = np.intersect1d(cpositions,fetched)

            if overlap.any():
                found_ov = True
                cpositions.extend(fetched)
                new_positions = np.unique(cpositions)
                epi_dict["Positions"][idx] = new_positions
                if cscore > epi_dict[ctool][idx]:
                    epi_dict[ctool][idx] = cscore
                break

        # new entry if no overlap
        if not found_ov:
            epi_dict["Positions"].append(cpositions)
            epi_dict["Sequence"].append("")
            epi_dict["Length"].append(0)
            for xx in tools:
                    epi_dict[xx].append(0)
            epi_dict[ctool][nepitopes] = cscore
            nepitopes += 1

    # Get sequence from positions
    for idx,cpos in enumerate(epi_dict["Positions"]):
        cseq = [sequence[xx] for xx in cpos]
        cseq = "".join(cseq)
        tmppos = [str(xx) for xx in cpos]
        tmppos = ",".join(tmppos)

        epi_dict["Sequence"][idx] = cseq
        epi_dict["Positions"][idx] = tmppos
        epi_dict["Length"][idx] = len(cseq)

    # Dump results
    epi_df = pd.DataFrame(epi_dict)
    epi_df.to_csv(ofile,index=False)
    print(f"Find your epitope regions at: {ofile}")


#  -Run Code -
# ------------
if __name__ == "__main__":
    args = options_parser()
    epiregions(args)
