#!/usr/bin/env python

'''
> make_dres_in_file.py start.pdb fixed_idr_residues.txt dres.in
'''

import sys
import numpy as np
import matplotlib.pyplot as plt
from soursop.sstrajectory import SSTrajectory
from soursop.ssprotein import SSProtein
import argparse

parser = argparse.ArgumentParser(description='create dres.in file')
parser.add_argument('pdb', type=str, help='input PDB file')
parser.add_argument('fixed_res_file', type=str, help='single line IDR fixed residues file')
parser.add_argument('outfile', type=str, help='output file')
parser.add_argument('--force-constant', type=float, default=20.0, 
                    help='Energy constant for distance restraint (default=20.0)')
parser.add_argument('--idr-first', action='store_true', help='flag if IDR is first chain in PDB')
args = parser.parse_args()

input_pdb = args.pdb
fixed_res_file = args.fixed_res_file
outfile = args.outfile

# Read in PDB as 1 frame trajectory
traj = SSTrajectory(trajectory_filename=input_pdb, pdb_filename=input_pdb).traj
prot = SSProtein(traj)

# Get distance map and contact map
dmap = prot.get_distance_map(verbose=False)[0]
cmap = np.ones((len(dmap), len(dmap))) * (dmap < 8.0)

# Get binary IDR fixed residues
with open(fixed_res_file) as f:
    lines = [x.strip().split() for x in f]
    idr_fixed = lines[0][1]

len_idr = len(idr_fixed)
fixed_start = idr_fixed.find('1')

idr_fixed = idr_fixed + '0'
fixed_end = idr_fixed[fixed_start:].find('0') + fixed_start
len_fd = len(dmap) - len_idr

if not args.idr_first:
	# Get residues where the fixed region of the IDR contacts the FD
	fd_CA_contacts, idr_CA_contacts = np.nonzero(cmap[:len_fd, len_fd+fixed_start:len_fd+fixed_end])

	# Adjust CA numbering for the IDR
	idr_CA_contacts += 1 + len_fd + fixed_start

else: # IDR first
	# Get residues where the fixed region of the IDR contacts the FD
	idr_CA_contacts, fd_CA_contacts = np.nonzero(cmap[fixed_start:fixed_end, len_idr:])

	# Adjust CA numbering for the IDR and FD
	idr_CA_contacts += 1 + fixed_start
	fd_CA_contacts += 2 + len_idr

# Write dres.in file from these contacts
with open(outfile, 'w') as out:
    out.write(f'{len(fd_CA_contacts)}\n')
    for r1, r2 in zip(fd_CA_contacts, idr_CA_contacts):
        a1 = prot.get_CA_index(r1)
        a2 = prot.get_CA_index(r2)
        out.write(f'{a1+1} {a2+1} 1 {prot.get_inter_residue_atomic_distance(r1, r2)[0]} {args.force_constant}\n')
