#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script to "rotate" and rsm dat file so that it's a list of omega scans instead of omega--2theta scans
#
# Author: Riccardo Civiero

import argparse
import os.path
import re
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

peak = np.argmax(intens)


# plt.plot(angle[0:np.size(angle)], intens[0:np.size(angle)], label = "{}".format(np.argmax(intens[0:np.size(angle)//2])))
# plt.legend()
# plt.show()

delta_angle = angle - angle[peak]


# cut_intens = intens[delta_angle > 0]
# delta_angle = delta_angle[delta_angle > 0]


# coefficients=np.polyfit(np.log(delta_angle)[9:21], np.log(cut_intens)[9:21], 1)
# polynomial=np.poly1d(coefficients)
# fit = lambda x: np.e**(coefficients[0]*x + coefficients[1])


# plt.loglog((delta_angle), (cut_intens), '+', label = "Scan")
# plt.loglog(delta_angle, fit(np.log(delta_angle)), label = "Fit slope = " + str(coefficients[0]))
# plt.xlabel(r"$\Delta$" + angleLabel)
# plt.ylabel("Intensity")
# plt.title(datfile)
# plt.legend()
# plt.show()

from scipy.optimize import curve_fit

i_peak = np.sum(intens)

f = lambda x, A, B: A * i_peak/(x**3)

popt, _ = curve_fit(f, delta_angle, intens)

print(popt)

plt.plot(delta_angle, intens, '+')
plt.plot(delta_angle, f(angle, *popt))
plt.show()