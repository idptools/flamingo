#!/usr/bin/env python

'''
> make_psw_file.py FD.isf fixed_idr_residues.txt output.psw
'''

import sys
import argparse

parser = argparse.ArgumentParser(description='create seq.in file')
parser.add_argument('fd_seq', type=str, help='FD amino acid sequence')
parser.add_argument('fixed_res_file', type=str, help='single line IDR fixed residues file')
parser.add_argument('outfile', type=str, help='output file')
parser.add_argument('--idr-first', action='store_true', help='flag if IDR is first chain in PDB')
parser.add_argument('--idr-caps', action='store_true', help='flag if IDR has ACE/NME caps')
args = parser.parse_args()

fd_seq = args.fd_seq
fixed_res_file = args.fixed_res_file
outfile = args.outfile

# Get FD aa sequence from isf file
# with open(fd_seqfile) as f:
#     lines = [x.strip().split() for x in f]
#     fd_seq = lines[0][0]

# Get binary IDR fixed residues
with open(fixed_res_file) as f:
    lines = [x.strip().split() for x in f]
    idr_fixed = lines[0][1]

fd_fixed = '1'*len(fd_seq)

if not args.idr_first:
    # Additional zeros are for ACE and NME caps 
    if args.idr_caps:
       full_fixed = fd_fixed + '0' + idr_fixed + '0'
    else:
       full_fixed = fd_fixed + idr_fixed
else:
    if args.idr_caps:
       full_fixed = '0' + idr_fixed + '0' + fd_fixed
    else:
       full_fixed = idr_fixed + fd_fixed

# Write file
fixed_line =    '0.0 0.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0'
idr_line =      '1.0 1.0 1.0 0.0 1.0 0.0 1.0 1.0 0.0 0.0 0.0'
s = 'R\n'
for i, state in enumerate(list(full_fixed)):
    if state == '1':
        s += f'{i+1}\t{fixed_line}\n'
    elif state == '0':
        s += f'{i+1}\t{idr_line}\n'
    else:
        raise ValueError('ERROR: invalid state!')

with open(outfile, 'w') as out:
    out.write(s)
