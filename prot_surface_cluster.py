#!/usr/bin/env python
# prot_surface_cluster.py
"""
Given a PDB, and a residue-size vector of properties [0,1] this script will
  1) Compute the MSMS surface
  2) Map values from the residue to the surface
  3) Cluster
     3.1) Based on geodesic distance
     or
     3.2) Based on connected elements
 4) Map back to PDB and group
 5) Dump PDB with group info (b-facor)
 6) Dump .csv with epitope information by cluster

"""

# ---------------------------
# Import Libraries
# ---------------------------
import os
from os.path import join
import sys
import ipdb
import argparse
import numpy as np
import pandas as pd

import subprocess
import ipdb

import prody
import meshio

# ---------------------------
# Parser Options
# ---------------------------
HELPTEXT = f"""

prot_surface_cluster.py [dev version]

Given a PDB, and a residue-size vector of properties [0,1] this script will
  1) Compute the MSMS surface
  2) Map values from the residue to the surface
  3) Cluster
     3.1) Based on geodesic distance
     or
     3.2) Based on connected elements
 4) Map back to PDB and group
 5) Dump PDB with group info (b-facor)
 6) Dump .csv with epitope information by cluster

Author:
------
Victor Montal
victor.montal [at] protonmail [dot] com \n

NOTES:
------
(1) MSMS MUST be installed and included in PATH env var
It is needed to generate xyzrn and mesh files

export MSMS_PATH=/path/to/downloaded/MSMS

try:
    echo $MSMS_PATH


(2) PDB MUST be clean (i.e no HEATM)
You can use something similar to
    grep "^ATOM" 1fv1_renum.pdb >1fv1_renum_clean.pdb

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
    h_ipdb= 'Path to input PDB'
    h_iprop = "Path to txt with residue-level feature/info.\n"
    h_iprop += "One row. Same number residues than PDB"
    h_outpath = "Output path"
    h_outname = "Output prefix name for dumped files"

    # Parser
    parser.add_argument('--pdb',
                        dest = 'ipdb', action = 'store',
                        help = h_ipdb, required = True)
    parser.add_argument('--feature',
                        dest = 'ifeat', action = 'store',
                        help = h_iprop, required = True)
    parser.add_argument('--outpath',
                        dest = 'outpath', action = 'store',
                        help = h_outpath, required = True)
    parser.add_argument('--outname',
                        dest = 'oname', action = 'store',
                        help = h_outname, required = True)

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

def compute_msms(pdb, opath):
    """
    Prepare the pdb to rxyzn and then run the MSMS to generate the surface
    """

def read_msms(file_root):
    """
    This script is a COPY from Masif repo
    """
    # read the surface from the msms output. MSMS outputs two files: {file_root}.vert and {file_root}.face
    vertfile = open(file_root + ".vert")
    meshdata = (vertfile.read().rstrip()).split("\n")
    vertfile.close()

    # Read number of vertices.
    count = {}
    header = meshdata[2].split()
    count["vertices"] = int(header[0])
    ## Data Structures
    vertices = np.zeros((count["vertices"], 3))
    normalv = np.zeros((count["vertices"], 3))
    atom_id = [""] * count["vertices"]
    res_id = [""] * count["vertices"]
    for i in range(3, len(meshdata)):
        fields = meshdata[i].split()
        vi = i - 3
        vertices[vi][0] = float(fields[0])
        vertices[vi][1] = float(fields[1])
        vertices[vi][2] = float(fields[2])
        normalv[vi][0] = float(fields[3])
        normalv[vi][1] = float(fields[4])
        normalv[vi][2] = float(fields[5])
        atom_id[vi] = fields[7]
        res_id[vi] = fields[9]
        count["vertices"] -= 1

    # Read faces.
    facefile = open(file_root + ".face")
    meshdata = (facefile.read().rstrip()).split("\n")
    facefile.close()

    # Read number of vertices.
    header = meshdata[2].split()
    count["faces"] = int(header[0])
    faces = np.zeros((count["faces"], 3), dtype=int)
    normalf = np.zeros((count["faces"], 3))

    for i in range(3, len(meshdata)):
        fi = i - 3
        fields = meshdata[i].split()
        faces[fi][0] = int(fields[0]) - 1
        faces[fi][1] = int(fields[1]) - 1
        faces[fi][2] = int(fields[2]) - 1
        count["faces"] -= 1

    assert count["vertices"] == 0
    assert count["faces"] == 0

    return vertices, faces, normalv, res_id

def sample_pdb_vertex(resid, pdb_feat):
    """
    Assign to each vertex, a residue metric
    """
    # Loop over each vertex and assign
    vtx_metric = []
    for idx,cresid in enumerate(resid):
        cid = int(cresid.split("_")[2])
        if pdb_feat[cid-1] == 1:
            vtx_metric.append(255)
        else:
            vtx_metric.append(0)
    return vtx_metric

def define_neighbours_vertex(faces):
    """
    Given a set of faces, create a dictionary of all neighbours
    for each vertex
    """
    neig_dict = {}
    for cface in faces:

        for cvtx in cface:
            if cvtx not in neig_dict.keys():
                neig_dict[cvtx] = []
            neig_dict[cvtx].extend(cface)
    # Remove duplicates
    for ckey in neig_dict.keys():
        neig_dict[ckey] = np.unique(neig_dict[ckey])
        pos_ckey = np.where(neig_dict[ckey] == ckey)[0]
        neig_dict[ckey] = np.delete(neig_dict[ckey],pos_ckey)

    return neig_dict


def seg_mesh_vtx(vertices, faces, vtx_metric):
    """
    Generate groups of vertex based on two conditions:
      1) Vertex metric must have value of 1
      2) One of vertex neighbour has value of one

    If none of the conditions is met, we assign such vertex to a new group
    """
    # Calculate vertex neighbours
    neig_dict = define_neighbours_vertex(faces)

    # Spread-like segmentation
    pos_mask = np.where(np.array(vtx_metric) > 0)[0]
    n_groups = 0
    vtx_group = np.zeros_like(vtx_metric)

    for cvtx in pos_mask:
        # Not assigned
        if vtx_group[cvtx] == 0:
            n_groups += 1
            visited = [cvtx]
            queue = [cvtx]
            while queue:
                # Select current working node
                node = queue.pop(0)
                # Get neighbours and their metrics
                cneigh = neig_dict[node]
                # Append neighbours with value to queue
                # and keep track of visited
                for neigh in cneigh:
                    if (neigh not in visited) and (vtx_metric[neigh] > 0):
                        queue.append(neigh)
                        visited.append(neigh)

            # Update vertex groups
            vtx_group[visited] = n_groups

    return vtx_group

def group_colors(vtx_group):
    """
    Define RGB color scheme for k different groups
    """
    kgroups = np.unique(vtx_group).size - 1
    ranks = np.round(np.linspace(0,255,kgroups)).astype(int)

    combinations = np.array(np.meshgrid(ranks, ranks, ranks)).T.reshape(-1, 3)
    rgb_row = np.random.choice(len(combinations), size=kgroups, replace=False)
    rgb_row = rgb_row + 1  # avoid (0,0,0) for non-groups
    rgb = np.ones([len(vtx_group),3]) * 128
    for idx,cgroup in enumerate(np.unique(vtx_group)[1:]):
        pos = np.where(vtx_group == cgroup)[0]
        ccolor = combinations[rgb_row[idx]]
        rgb[pos,:] = ccolor

    return rgb

def groups2pdb(res_id, segment_groups, pdb):
    """
    Given a set of surface groups, project back to residue (b-factor/Beta in Prody)
    If more than one group to one residue, select the most occurance
    """
    newpdb = pdb.copy()

    all_surf_res = [int(xx.split("_")[2]) for xx in res_id]
    surf_res = np.unique(all_surf_res)

    pdb_residues = newpdb.getResnums()
    betas = np.zeros_like(newpdb.getBetas())

    for cres in np.unique(pdb_residues):

        if cres in surf_res:
            pos_res = np.where(all_surf_res == cres)[0]
            pdb

            groups = [segment_groups[xx] for xx in pos_res ]
            if any(groups):
                uniq,count = np.unique(groups, return_counts=True)

                if uniq[0] == 0:
                    uniq = uniq[1:]    # ignore zeros
                    count = count[1:]  # ignore zeros.

                max_pos = np.where(count == np.max(count))[0]
                max_group = uniq[max_pos]

                pdb_pos = np.where(pdb_residues == cres)[0]
                betas[pdb_pos] = max_group*10

    newpdb.setBetas(betas)
    return newpdb

def pdb2csv(pdb):
    newpdb = pdb.select('calpha')
    groups = newpdb.getBetas()
    all_resid = newpdb.getResnums()
    all_resname = newpdb.getResnames()

    # Re-encode residue namess
    aa_dict = {"ALA" : "A",
               "ARG" : "R",
               "ASN" : "N",
               "ASP" : "D",
               "CYS" : "C",
               "GLN" : "Q",
               "GLU" : "E",
               "GLY" : "G",
               "HIS" : "H",
               "ILE" : "I",
               "LEU" : "L",
               "LYS" : "K",
               "MET" : "M",
               "PHE" : "F",
               "PRO" : "P",
               "SER" : "S",
               "THR" : "T",
               "TRP" : "W",
               "TYR" : "Y",
               "VAL" : "V"}
    all_resname = [ aa_dict[xx] for xx in all_resname]

    csvdict = {}
    csvdict["Rank"] = []
    csvdict["Sequence"] = []
    csvdict["Start"] = []
    csvdict["End"] = []
    csvdict["Positions"] = []
    csvdict["Score"] = []

    for cgroup in np.unique(groups)[1:]:
        pos_group = np.where(groups == cgroup)[0]
        cresi_names = [all_resname[xx] for xx in pos_group]
        cresi_pos = [str(all_resid[xx]) for xx in pos_group]

        csvdict["Rank"].append(cgroup)
        csvdict["Sequence"].append("".join(cresi_names))
        csvdict["Start"].append("NA")
        csvdict["End"].append("NA")
        csvdict["Positions"].append(",".join(cresi_pos))
        csvdict["Score"] = 1

    dfcsv = pd.DataFrame(csvdict)
    return dfcsv

# ---------------------------
# Main Code
# ---------------------------
def cluster_pdb(args):
    print(f"Working with pdb: {args.ipdb}")
    # DEFAULTS
    msms_path = os.environ["MSMS_PATH"]
    opath = args.outpath


    print(f"    load files")
    pdb = prody.parsePDB(args.ipdb)
    feat = np.loadtxt(args.ifeat)
    # Check size residues pdb and feat agree


    print(f"    generate surface mesh")
    cwd = os.getcwd()    # to make pdb2xyzrn work
    os.chdir(msms_path)
    # pdb2xyzrn
    xyzrn = join(opath,f"{args.oname}_pdb2xyzrn.txt")
    cmd = f""
    cmd += f"{msms_path}/pdb_to_xyzrn {args.ipdb}"
    cmd += f" > {xyzrn}"
    run_cmd(cmd,"Could not convert PDB to XYZRN file")

    osurf = join(opath,f"{args.oname}")
    cmd = f""
    cmd += f"{msms_path}/msms.x86_64Linux2.2.6.1 "
    cmd += f"-if {xyzrn} -of {osurf} "
    print(cmd)
    run_cmd(cmd,"Could not convert XYZRN to surface files")
    os.chdir(cwd)


    print(f"    sampling form PDB to surface")
    vertices, faces, normalv, res_id = read_msms(osurf)
    vtx_metric = sample_pdb_vertex(res_id, feat)
    # Dump temporal surface
    points = list(vertices)
    cells = [ ("triangle", faces)]
    mesh = meshio.Mesh(points,
                       cells,
                       point_data={"red": np.array(vtx_metric).astype(np.uint8)})
    oply = join(opath,f"tmp.{args.oname}.ply")
    mesh.write(oply)
    cmd = f""
    cmd += f"meshio ascii {oply}"
    run_cmd(cmd, "Could not force .ply file into ascii")

    print(f"    segment mesh")
    segment_groups = seg_mesh_vtx(vertices, faces, vtx_metric)
    rgb_groups = group_colors(segment_groups)

    print(f"    dump mesh with groups")
    mesh = meshio.Mesh(points,
                       cells,
                       point_data={"red": np.array(rgb_groups[:,0]).astype(np.uint8),
                                   "green": np.array(rgb_groups[:,1]).astype(np.uint8),
                                   "blue": np.array(rgb_groups[:,2]).astype(np.uint8)}
                       )
    oply = join(opath,f"groups.{args.oname}.ply")
    mesh.write(oply)
    cmd = f""
    cmd += f"meshio ascii {oply}"
    run_cmd(cmd, "Could not force .ply file into ascii")


    print(f"    map surface segments to residues")
    newpdb = groups2pdb(res_id,segment_groups,pdb)


    print(f"    dump PDB")
    opdb = join(opath,f"groups.{args.oname}.pdb")
    prody.writePDB(opdb,newpdb)


    print(f"    dump csv")
    dfcsv = pdb2csv(newpdb)
    odf = join(opath,f"groups.{args.oname}.csv")
    dfcsv.to_csv(odf)


#  -Run Code -
# ------------
if __name__ == "__main__":
    args = options_parser()
    cluster_pdb(args)
