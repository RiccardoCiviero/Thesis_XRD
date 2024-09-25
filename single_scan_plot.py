#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script to "rotate" and rsm dat file so that it's a list of omega scans instead of omega--2theta scans
#
# Author: Riccardo Civiero

import argparse
import math
import os.path
import re
import string
import sys
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Plot a single dat file in one command")
parser.add_argument('datfile', help=".dat rsm filename")

args = parser.parse_args()

datfile=args.datfile

datre = re.compile('.dat')

# load the .dat file

def contains_rotated(datfile): # to recognize if I am reading a "normal" omega2theta-scans dat file or a "rotated" omega-scans one 
    with open(datfile) as file:
        for _  in range(3):
            line = file.readline()
        return "rotated" in line

if os.path.isfile(datfile):
    if 'xrdml' in datfile:
        xrdmlfile=datfile
        xre = re.compile('.xrdml')
        datfile = xre.sub('.dat', xrdmlfile)
        if os.path.isfile(datfile):
            print('running on ' + datfile + ' associated with xrdml file ' + xrdmlfile)
            
        else:
            print(xrdmlfile + " exists but " + datfile + " does not exist, run xrd.py")
            exit(1)

    if os.path.getsize(datfile) > 0:
        # omega, ttheta, phi, chi, xpos, ypos, zpos, qx, qz, intens = np.loadtxt(datfile, dtype=float, skiprows=8, unpack=True)
        omega, ttheta, phi, chi, xpos, ypos, zpos, qx, qz, intens = np.genfromtxt(datfile, dtype=float, skip_header=7, unpack=True)
        # reloaded = pd.read_table(datfile, dtype=float, skiprows=8, skip_blank_lines=True, names=["omega", "ttheta", "phi", "chi", "xpos", "ypos", "zpos", "qx", "qz", "intens"])
        if contains_rotated(datfile):
            scantype = "omega"
        else:
            scantype = "omega2theta"
        
    else:
        print(datfile + " is empty")
        exit(1)
        
else:
    print(datfile + " does not exist")
    exit(1)


if scantype == "omega":
    angle = omega
    angleLabel = r"$\omega$"
else:
    angle = ttheta # ??? Not really sure will check if ever needed
    angleLabel = r"$\omega$ - $2\theta$"

plt.plot(angle, intens)
plt.xlabel(angleLabel)
plt.ylabel("Intensity")
plt.title(datfile)
plt.show()
