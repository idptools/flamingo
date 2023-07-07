#!/usr/bin/env python

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--out", default='run_sims.py', required=True,
                    help="name of output python script")
parser.add_argument("--keyfile", metavar='KEYFILE', required=True, 
                    default="run.key", help="name of keyfile")
parser.add_argument("--preeq", metavar='N_PRE_EQ', required=True, 
                    default=500000, type=int, help="number of pre-equilibration steps")
parser.add_argument("--preeq-helix", action='store_true',
                    help="start pre-equilibration from helical structure")
parser.add_argument("--eq", metavar='N_EQ', required=True, 
                    default=2000000, type=int, help="number of equilibration steps")
parser.add_argument("--mc", metavar='N_MC', required=True, 
                    default=50000000, type=int, help="number of non-equilibration MC steps")
parser.add_argument("--mode", choices=['standard', 'ev', 'ts', 'hs'], required=True,
                    help="MC simulation mode: standard, EV, temperature sweep (TSMC), or Hamiltonian switch (HSMC)")
parser.add_argument("--aux-enter-prob", metavar='p', default=0.8, type=float, 
                    help="probability of entering auxiliary chain")
parser.add_argument("--aux-chain-freq", metavar='N', type=int, 
                    help="number of MC steps between auxiliary chain attempts for TSMC or HSMC")
parser.add_argument("--aux-chain-steps", metavar='N', type=int, 
                    help="number of steps for each auxiliary chain, multiplied by number of distinct temps for TSMC. 500 recommended for TSMC, 50 recommended for HSMC")
parser.add_argument("--temp", metavar='TEMP',
                    type=int, help="simulation temperature (only needed for TSMC)")

args = parser.parse_args()

# Check arguments
if args.aux_enter_prob is not None:
    if args.aux_enter_prob > 1 or args.aux_enter_prob < 0:
        raise Exception('auxiliary chain enter probability must be 0 <= p <= 1')

if args.mode == 'ts' or args.mode == 'hs':
    # Require auxiliary chain frequency and number of steps
    if args.aux_chain_freq is None:
        raise Exception(f'For {args.mode}, you must specify the auxiliary chain frequency (--aux-chain-freq)')
    if args.aux_chain_steps is None:
        raise Exception(f'For {args.mode}, you must specify the number of auxiliary chain steps (--aux-chain-steps)')
    
    if args.mode == 'ts':
        # Must provide simulation temp for TSMC
        if args.temp is None:
            raise Exception(f'For TSMC, you must specify the number of simulation temperature (--temp)')


# Begin writing output script
script_str = "from MonteCarlo import Protein,MonteCarlo\nimport random\n\n"

# Pre-Eq / Equilibration steps required regardless of simulation type
script_str += f"simulation = MonteCarlo.MonteCarlo('{args.keyfile}')\n"
script_str += f"simulation.preEquilMC(steps={args.preeq},helix={args.preeq_helix},debug=True)\n"
script_str += f"simulation.EquilMC(steps={args.eq},debug=True)\n"

if args.mode == 'standard' or args.mode == 'ev':
    script_str += f"simulation.StandardMC(steps={args.mc},debug=True)\n"

else: # HSMC or TSMC
    n_attempts = int(np.ceil(args.mc / args.aux_chain_freq))
    script_str += f"for i in range({n_attempts}):\n"
    script_str += f"\tif i == 0:\n"
    script_str += f"\t\tsimulation.StandardMC(steps={args.aux_chain_freq},debug=True)\n"
    script_str += f"\telse:\n"
    script_str += f"\t\tsimulation.StandardMC(steps={args.aux_chain_freq},debug=False)\n"
    script_str += f"\trandom.seed()\n"
    script_str += f"\tif {args.aux_enter_prob} > random.uniform(0,1):\n"

    if args.mode == 'hs':
        script_str += f"\t\tsimulation.HamiltonianSwitchMC(steps={args.aux_chain_steps})\n"

    else: # Temp sweep
        temp_list = [args.temp+10, args.temp+1000, args.temp+500, args.temp+250, 
                     args.temp+100, args.temp+80, args.temp+60, args.temp+40, 
                     args.temp+20, args.temp+10]
        str_temp_list = str(temp_list)
        script_str += f"\t\tsimulation.TempSweepMC(steps={args.aux_chain_steps},auxiliary_temperatures={str_temp_list})\n"

# print(script_str)

with open(args.out, 'w') as f:
    f.write(script_str)

