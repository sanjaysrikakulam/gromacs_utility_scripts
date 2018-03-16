#! python2
"""
Description: Python wrapper for computing RMSD (on backbone atoms only) using GROMACS and 
					automatically plots the resulting output from gmx rms 
							(tested on GROMACS version 5.1.*)
Author: Sanjay kumar Srikakulam (sanjaysrikakulam@gmail.com)
WARNING: PLEASE MAKE SURE YOU HAVE THE 'PLOT_XVG.PY' SCRIPT AVAILABLE IN THE SAME FOLDER AS THIS SCRIPT!!!
Init: 15.03.2018
"""

# Import necessary libraries
import os
import argparse
import subprocess
from plot_xvg import plotxvg

# Input parsing and processing
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", dest = 'structure_file', help = "reference file (e.g., *.tpr or *.pdb)", required = True )
    parser.add_argument("-f", dest = 'trajectory', help = "Trajectory file (*.xtc)", required = True)
    parser.add_argument("-d", dest = 'out_dir', nargs = '?', help = "Path to the output directory", default = os.getcwd())

    args = parser.parse_args()
    structure_file = args.structure_file
    trajectory = args.trajectory
    out_dir = args.out_dir
    
    # Computation and plotting of RMSD with least squares fit to backbone
    subprocess.call("echo 4 4 | gmx rms -s {0} -f {1} -o {2}'/backbone_rmsd.xvg' -tu ns".format(structure_file, trajectory, out_dir), shell = True)
    plotxvg("backbone_rmsd.xvg", out_dir)

# Function call
if __name__=='__main__':
    main()

