#!/usr/bin/env python

import sys

seq = sys.argv[1]
fixed = sys.argv[2]

flexible_residues = []
for state, res in zip(fixed, seq):
    if state == '0' and res == 'P':
        sys.exit(0)


sys.exit(1)