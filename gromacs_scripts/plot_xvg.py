#! python2
"""
Description: A program that reads and plots the GROMACS XVG file (e.g., RMSD, RMSF) 
					(tested on GROMACS version 5.1.*)
Author: Sanjay kumar Srikakulam (sanjaysrikakulam@gmail.com)
WARNING: PLEASE MAKE SURE YOU HAVE ALL THESE IMPORTED LIBRARIES INSTALLED!!!
Init: 14.03.2018
"""

# Import necessary libraries
import os
import argparse
import re
import matplotlib
matplotlib.use('Agg')					#Unset matplotlib Xwindows default
import matplotlib.pyplot as plt

#Input parsing and processing
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", dest = 'xvg_file', help = "XVG file (e.g., *.xvg)", required = True )
    parser.add_argument("-d", dest = 'out_dir', nargs = '?', help = "Path to the output directory", default = os.getcwd())

    args = parser.parse_args()
    xvg_file = args.xvg_file
    out_dir = args.out_dir
    
    # Function call
    plotxvg(xvg_file, out_dir)
    
# GROMACS XVG file parsing and plotting (tested with the XVG files from gmx rms, and gmx rmsf)
def plotxvg(filename, out_dir):   
    with open(filename, 'r') as file:
        title = ""
        xyaxis_labels = list()
        data_x = list()
        data_y = list()
        for line in file:
            if line.find("@") != -1:

                # Extracts the title
                if re.match("title", line.split()[1]):
                    title = line.split()[2].strip('\"')

                # Extracts the legends for the axes
                if len(line.split()) > 3:
                    if re.match("label", line.split()[2]):
                        xyaxis_labels.append(line.split('"')[1])

            if not line.startswith("@") and not line.startswith("#"):
                                  data_x.append(float(line.split()[0]))
                                  data_y.append(float(line.split()[1]))

        # Creates and saves the plot
        plt.plot(data_x, data_y, 'k-', linewidth = 0.4)
        plt.xlabel(xyaxis_labels[0])
        plt.ylabel(xyaxis_labels[1])
        plt.title(title)
        plt.xticks(range(int(min(data_x)), int(max(data_x) + (int(max(data_x))-int(min(data_x))) / 10), (int(max(data_x))-int(min(data_x)))/10))
        out_file = "/".join([out_dir, filename.replace("xvg", "png")])
        plt.savefig(out_file, dpi = 300)
        
# Function call
if __name__ == '__main__':
    main()