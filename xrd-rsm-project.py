#! /usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Script to extract some of the scans from an rsm dat file and collapse them into a single .dat file

import argparse
import os.path
import re
import numpy as np

parser = argparse.ArgumentParser(description="Extract some of the scans from an rsm .dat file and collapse them into a single  .dat file")
parser.add_argument("-s", "--start", nargs=1,
                  dest="startnum", default=[0], type=int,
                  help="which scan to start with")
parser.add_argument("-e", "--end", nargs=1,
                  dest="endnum", default=[200], type=int,
                  help="which scan to end with")
parser.add_argument("-d", "-destination", default = "", type = str, dest ="des", help="Destination folder")
parser.add_argument('datfile', help=".dat rsm filename")
args = parser.parse_args()
datfile=args.datfile
filename = datfile.split('/')[-1]

if args.des != "":
    if not os.path.isdir(args.des):
        print("Destination directory {} does not exist!".format(args.des))
        exit(1)

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
        omega, ttheta, phi, chi, xpos, ypos, zpos, qx, qz, intens = np.genfromtxt(datfile, dtype=float, skip_header=8, unpack=True)
        # reloaded = pd.read_table(datfile, dtype=float, skiprows=8, skip_blank_lines=True, names=["omega", "ttheta", "phi", "chi", "xpos", "ypos", "zpos", "qx", "qz", "intens"])
        if contains_rotated(datfile):
            scantype = "omega"
        else:
            scantype = "omega2theta"
        print(scantype)
        
    else:
        print(datfile + " is empty")
        exit(1)
        
else:
    print(datfile + " does not exist")
    exit(1)

# count the number of scans in the dat file and the length of each scan
i = 0
j = 0
maxj = 0
scans = 1
while i<len(intens):
    if(i+1<len(intens)):
        if(omega[i+1]<omega[i]):
            # print(omega[i],omega[i+1], j)
            scans=scans+1
            if j>maxj:
                maxj=j

            j=-1

    i=i+1
    j=j+1

points=maxj+1

print('There are {} scans in {}, each one has {} points'.format(scans,datfile,points))

# load all the data into an array
datarray=np.zeros((scans,10,points))
i = 0
j = 0
scan = 0
while i<len(intens):
    datarray[scan, :, j]=omega[i], ttheta[i], phi[i], chi[i], xpos[i], ypos[i], zpos[i], qx[i], qz[i], intens[i]
    # print(i,scan,j, omega[i], omega[i+1])
    if(i+1<len(intens)):
        if(omega[i+1]<omega[i]):
            scan=scan+1
            j=-1

    i=i+1
    j=j+1

#print(np.shape(datarray))

startnum=args.startnum[0]

if args.endnum[0]>scans:
    endnum=scans-1
else:
    endnum=args.endnum[0]

if startnum>endnum:
    print('starnum {} > endnum {}'.format(startnum,endnum))
    exit(1)

numscans=endnum-startnum

# select the scans to be averaged:
selarray=datarray[startnum:(endnum+1), :, :]
#print(np.shape(selarray))

# perform the average:
sumarray=np.sum(selarray,0)/float(numscans+1)
#print(np.shape(sumarray))
#print(sumarray)

zfiller=4
if scans < 1000:
    zfiller=3

if scans < 100:
    zfiller=2

if scans < 10:
    zfiller=1

out_string=('_{}-{}'.format(str(startnum).zfill(zfiller),str(endnum).zfill(zfiller)))
#print(out_string)
if scantype == "omega":
    outfile = datre.sub('_w'+out_string+'.dat', args.des + "/" + filename)
else:
    outfile = datre.sub('_w'+out_string+'.dat', args.des + "/" + filename)

print('saving to '+outfile)

with open(outfile, 'w') as datf:
    datf.write ('# filename: ' + filename + '\n')
    datf.write ('# reflection: (004)\n')
    if scantype == "omega":
        datf.write ('# type: Projection of {} {} scans\n'.format(endnum-startnum, scantype + " (rotated) "))
    if scantype == "omega2theta":
        datf.write ('# type: Projection of {} {} scans\n'.format(endnum-startnum, scantype))
    datf.write ('# x-axis: {} \n'.format(scantype))
    datf.write ('# \n')
    datf.write ('# Omega\t2Theta\tPhi  \tChi   \tX    \tY    \tZ    \tq_para\tq_perp\tIntensity\n')
    datf.write ('# [deg]\t[deg] \t[deg]\t[deg] \t[mm] \t[mm] \t[mm] \t[nm-1]\t[nm-1]\t[counts/s]\n')
    j = 0
    while j < points:
        i = 0
        while i < 10:
            if i < 2:
                sout_val = '{0:.4f}'.format(sumarray[i,j])
                
            elif i < 6:
                sout_val = '{0:.2f}'.format(sumarray[i,j])
                
            elif i == 6:
                sout_val = '{0:.3f}'.format(sumarray[i,j])
            
            elif i < 9:
                sout_val = '{0:.5f}'.format(sumarray[i,j])
                
            else:
                sout_val = '{0:.1f}'.format(sumarray[i,j])
            
            datf.write(sout_val+'\t')
            i = i + 1

        j = j + 1
        datf.write("\n")            

    
datf.close()




