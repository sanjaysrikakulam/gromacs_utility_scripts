#! python2
"""
Description: A program that reads and converts the Secondary structure XPM file from GROMACS 'gmx do_dssp' to XVG file 
		(calculates the secondary structure propensity per residue) and also produces a plot of the same. 
								(tested on GROAMACS version 5.1.*)
Author: Sanjay kumar Srikakulam (sanjaysrikakulam@gmail.com)
Note: This script is based on Dr. Justin Lemkul's 'plot_ss_xmp.pl' script (http://www.bevanlab.biochem.vt.edu/Pages/Personal/justin/scripts.html)
WARNING: PLEASE MAKE SURE YOU HAVE ALL THESE IMPORTED LIBRARIES INSTALLED!!!
Init: 10.03.2018
"""

# Import necessary libraries
from __future__ import division
import os
import re
import subprocess
import argparse
import matplotlib
matplotlib.use('Agg')					#Unset matplotlib Xwindows default
import matplotlib.pyplot as plt
import os.path

# User input parsing and processing
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", dest = 'xpm_file', help = "Secondary structure XPM file from GROMACS (e.g., *.xpm)", required = True )
    parser.add_argument("-d", dest = 'out_dir', nargs = '?', help = "Path to the output directory", default = os.getcwd())

    args = parser.parse_args()
    xpm_file = args.xpm_file
    out_dir = args.out_dir
    xvg_out_file = "/".join([out_dir, xpm_file.replace("xpm", "xvg")])
    
    # Function call
    xpm_to_xvg(xpm_file, out_dir)
    plot_xvg(xvg_out_file)
   
# Parses and converts the GROAMACS (5.1.*) produced secondary structure XPM file to XVG file (calculates secondary structure propensity per residue)
def xpm_to_xvg(ss_xpm_file, out_dir): 
    with open(ss_xpm_file, 'r') as file:
        y_axis = list()
        number_of_frames = int()
        number_of_residues = int()
        xvg_out_file = "/".join([out_dir, ss_xpm_file.replace("xpm", "xvg")])
        
        # Parses the XPM file and keeps the necessary secondary structure information
        for line in file:
            if line.startswith("static char"):
                number_of_frames, number_of_residues = next(file).split()[0:2]
                number_of_frames = int(number_of_frames.strip('\"'))
                number_of_residues = int(number_of_residues)
                
            if len(line.split()) > 1:
                if re.match("y-axis", line.split()[1]):
                    y_axis = line.split()[2:-1]
                    
        # Keeps only the secondary structure data and discards the rest (comments, x- & y-axes information, etc)
        lines = file.readlines()
        lines = subprocess.check_output("cat {0} | sed '0,/y-axis/d'".format(ss_xpm_file), shell = True)
        lines = lines.replace('"', '').strip().split(",")
        lines = map(lambda s: s.strip(), lines)
        
        # Secondary structure types initialization
        alpha_helix = {key: 0 for key in range(1, number_of_residues+1)}                #(H)
        three_10_helix = {key: 0 for key in range(1, number_of_residues+1)}             #(G)
        pi_helix = {key: 0 for key in range(1, number_of_residues+1)}                   #(I)
        beeta_sheet = {key: 0 for key in range(1, number_of_residues+1)}                #(E)
        coil = {key: 0 for key in range(1, number_of_residues+1)}                       #(~)
        turn = {key: 0 for key in range(1, number_of_residues+1)}                       #(T)
        bridge = {key: 0 for key in range(1, number_of_residues+1)}                     #(B)
        bend = {key: 0 for key in range(1, number_of_residues+1)}                       #(S)
        
        # Parses the secondary structure data in the XPM file
        for i in range(number_of_residues, 0, -1):
            j = number_of_residues - i        
            for k in range(0, number_of_frames):
                if "H" in lines[j][k]:
                    alpha_helix[i] += 1
                elif "G" in lines[j][k]:
                    three_10_helix[i] += 1
                elif "I" in lines[j][k]:
                    pi_helix[i] += 1
                elif "E" in lines[j][k]:
                    beeta_sheet[i] += 1
                elif "~" in lines[j][k]:
                    coil[i] += 1
                elif "T" in lines[j][k]:
                    turn[i] += 1
                elif "B" in lines[j][k]:
                    bridge[i] += 1
                elif "S" in lines[j][k]:
                  bend[i] += 1

        # Creates the new XVG file containing secondary structure propensity per residue
        if not os.path.exists(xvg_out_file):
            with open(xvg_out_file, "a") as op:
                op.write("# Author: Sanjay kumar Srikakulam (sanjaysrikakulam@gmail.com)\n")
				op.write("# Information: XVG file produced from the script convert_ss_xpm_to_png.py (https://github.com/sanjaysrikakulam/gromacs_utility_scripts/gromacs_scripts/convert_ss_xpm_to_png.py) \n")
                op.write("# Description: Probability of various secondary structure elements, by residue.\n")
                op.write("@\ttitle\t\"Secondary Structure Content\"\n")
                op.write("@\txaxis\tlabel \"Residue\"\n")
                op.write("@\tyaxis\tlabel \"Probability\"\n")
                op.write("@TYPE xy\n")
                op.write("@ s0 legend \"\\f{Symbol}a\\f{Times}-Helix\"\n")
                op.write("@ s1 legend \"\\f{Symbol}p\\f{Times}-Helix\"\n")
                op.write("@ s2 legend \"3\\s10\\N-Helix\"\n")
                op.write("@ s3 legend \"\\f{Symbol}b\\f{Times}-Strand\"\n")
                op.write("@ s4 legend \"\\f{Symbol}b\\f{Times}-Bend\"\n")
                op.write("@ s5 legend \"\\f{Symbol}b\\f{Times}-Turn\"\n")
                op.write("@ s6 legend \"\\f{Symbol}b\\f{Times}-Bridge\"\n")
                op.write("@ s7 legend \"Coil\"\n")
        
                for z in range(1, number_of_residues + 1):
                    op.write("%s\t\t%.6f\t\t%.6f\t\t%.6f\t\t%.6f\t\t%.6f\t\t%.6f\t\t%.6f\t\t%.6f\t\t\n" % (y_axis[z-1], alpha_helix[z]/number_of_frames, pi_helix[z]/number_of_frames, three_10_helix[z]/number_of_frames, beeta_sheet[z]/number_of_frames, bend[z]/number_of_frames, turn[z]/number_of_frames, bridge[z]/number_of_frames, coil[z]/number_of_frames))
        else:
            print "The output XVG file {0} already exists. Aborting!".format(xvg_out_file)
            quit()

# Parses and plots the secondary structure output from the above function (xpm_to_xvg)
def plot_xvg(ss_xvg_file):
    with open(ss_xvg_file, 'r') as file:
        
        # Initialize and parse the input file
        data_c2 = list()                #(alpha_helix)
        data_c3 = list()                #(pi_helix)
        data_c4 = list()                #(three_10_helix)
        data_c5 = list()                #(beeta_sheet)
        data_c6 = list()                #(bend)
        data_c7 = list()                #(turn)
        data_c8 = list()                #(bridge)
        data_c9 = list()                #(coil)
        data_y = list()                 #(residue_list)
        
        for line in file:          
            if not line.startswith("@") and not line.startswith("#"):
                data_c2.append(float(line.split()[1]))
                data_c3.append(float(line.split()[2]))
                data_c4.append(float(line.split()[3]))
                data_c5.append(float(line.split()[4]))
                data_c6.append(float(line.split()[5]))
                data_c7.append(float(line.split()[6]))
                data_c8.append(float(line.split()[7]))
                data_c9.append(float(line.split()[8]))
                data_y.append(int(line.split()[0])) 

        # Plotting the parsed data
        fig, ax = plt.subplots()
        ax.plot(data_y, data_c2, 'b-', linewidth = 0.6, label = r'$\alpha$'"-helix")
        ax.plot(data_y, data_c3, 'm-', linewidth = 0.6, label = "pi-helix")
        ax.plot(data_y, data_c4, 'g--', linewidth = 0.6, label = r'$3_{10}$'"-helix")
        ax.plot(data_y, data_c5, 'r:', linewidth = 0.6, label = r'$\beta$'"-sheet")
        ax.plot(data_y, data_c6, 'brown', linewidth = 0.6, label = r'$\beta$'"-bend")
        ax.plot(data_y, data_c7, 'orange', linewidth = 0.6, label = r'$\beta$'"-turn")
        ax.plot(data_y, data_c8, 'yellow', linewidth = 0.6, label = r'$\beta$'"-bridge")
        ax.plot(data_y, data_c9, 'c-', linewidth = 0.6, label = "Coil")            
        plt.xlabel("Residue")
        plt.ylabel("Probability")
        plt.title("Secondary structure content")
        legends = ax.legend(bbox_to_anchor=(1, 1.01), loc=2)
        png_file_name = ss_xvg_file.replace("xvg", "png")
        plt.savefig(png_file_name, dpi = 300, bbox_extra_artists=(legends,), bbox_inches='tight')
        
# Function call
if __name__ == '__main__':
    main()
