# flamingo
Repo for running and analyzing simulations of IDR motifs bound to folded domains

### For Jeff & Borna - 7/7/23

The requisite scripts for setup and analysis of FD+IDR simulations should all be in the repository. The most time-consuming step of the process is minimization of an AF2 structure in preparation for MC simulations. Fortunately, for our TAZ2 systems, we should not need to do any additional minization. All of the structures we need should be located in my `/work/degriffith/TADS/p300_coactivators` repo (I can help you find specific files you might need). 

For setting up and running these sims, you will need to copy the `setup/build_FD_IDR_sim_infrastructure_v1.sh` into a clean simulation directory. This script requires a few bits of information to run. 

1. The location of the locally cloned flamingo repo
2. The path to an "IDR sequence" file. This is a whitespace separated file, one line per distinct type of simulation. It has the following format:

```
IDR_seq_ID1    <AA_seq1>
IDR_seq_ID2    <AA_seq2>
 .
 .
 .
IDR_seq_IDN	   <AA_seqN>
```

3. The path to a "fixed IDR residues" file. This has the same format as the the IDR sequence file (one line per simulation), but instead of amino acid sequences, the IDR is represented as a string of 0's and 1's, where 1's denote which residues should have fixed backbones in the simulation (i.e., the 'motif'). I have found that it works well to define the fixed residues as the motif *minus two residues on both ends*. See some of the example simulation directories in my `/work` repo for inspiration here (or ask me).

4. The path to a "FD sequence" file. Same format as the IDR sequence file, however now there is only one line. Note that the current setup is only designed to run one distinct FD per `build_infrastructure` script, so if you want to compare distinct FDs, you will need to run them in separate directories with separate FD sequence files and structure files.

5. Finally, you need to provide a path to a "PDB structure list file". This also has a format of one-line-per-simulation. Like #2 and #3, you need the first column in this whitespace-separated file to be the sequence IDs. The second column contains paths to the input PDB structures containing the bound FD-IDR complex. **Important: these structures MUST contain exactly the residues indicated as "fixed" in your fixed IDR residues file, no more, no less. If you include ACE/NME caps, the simulation will not run properly (but annoyingly it will still run, just not with the chains bound).**


With those files in place, you can also cahnge the other main CAMPARI simulation parameters. The important things to consider here are:

 - MC_MODE : this is where you denote temperature sweep, standard, or EV
 - TEMPERATURE : the base temperature 
 - REPS : number of replicate simulations per simulation ID (reps per line in files #2,#3,#5)
 - SIM_MODE : starting state of simulations is either coil, helical, or both. I've just been running coil_start
 - AUXCHAIN_* : these are the temperature sweep params 

If you want to modify the temperature sweep parameters, you will need to manually modify the file `setup/create_tsmc_run_script.py`

That should mostly be it! Let me know if you have questions.