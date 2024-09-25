#! /usr/bin/env python3
# -*- coding: utf-8 -*-

# Script to process .xrdml files from the PANalytical diffractometer
# This version 2.0 for Python 3
# 12/07/2017, -v flag added 15/03/2018, fixed the dat file number formats
# Danny Chrastina

# written in Python 2.6.5
# requires gnuplot 4.4 to plot the data in .ps files
# imagemagick mogrify and convert to create .jpegs
# and gzip to compress the .ps files

# version 1.1: added support for reciprocal space line scans,
# version 1.11: and -f and -n flags to choose whether to reprocess unchanged files
# version 2.0: converted to work with Python 3

# version 2.1: 22/08/2018, can't use Imagemagick due to https://www.kb.cert.org/vuls/id/332928
# call either the xrd-png or xrd-post bash scripts

# version 3.0: 29/10/2018, support for tilt correction
#     specify either -t/--tilt or -o/--omega and -q/--twotheta, tilt = omega-twotheta/2

# 11/06/2019, save gnuplot file with plt extension for windows compatibility
#
# version 3.1: 13/06/2019, an option to plot omega-2theta scans with
# theta or 2theta as the x-axis
#     specify either --w2t or --2tw respectively

# 26/11/2021, "set autoscale fix" for RSMs to allow map to fill more of the plot
# we also now use the MATLAB palette by default

# version 3.2: 09/01/2023, set the yrange bottom limit depending on whether
# * - it's a long triple-axis scan (use 0.1)
# * - it's a rocking-mode scan (autorange) except there's nothing in the xrdml
#     to indicate which kind of scan it is, if you can believe that
# * - otherwise 1.0
# * - so for now just autorange

import argparse
import math
import os.path
import re
import string
import sys
import time
from subprocess import call

def grep(string,mylist):
    """grep string in mylist
from http://casa.colorado.edu/~ginsbura/pygrep.htm"""
    expr = re.compile(string)
    return list(filter(expr.search,mylist))
    
def grepall(string,mylist):
    """grep all occurrences of a string in mylist"""
    expr = re.compile(string)
    return list(filter(expr.findall,mylist))

def grepn(string,mylist):
    """grepn string in mylist returns the index of the first occurrence of 
that string in the mylist"""
    for n in range(len(mylist)):
        if (re.search(string,mylist[n]) is not None):
            return n

def grepalln(string,mylist):
    """grepn string in mylist returns the indices of all the occurrences 
of that string in the mylist"""
    lines = []
    for n in range(len(mylist)):
        if (re.search(string,mylist[n]) is not None):
            lines = lines[:] + [n]
    return lines

def grepxml(tag,list):
    """grepxml tag in list returns the stuff between <tag> and </tag>"""
    nopen = grepn('<' + tag + '>',list)
    nclose = grepn('</'+ tag + '>',list)
    return list[nopen:nclose]

def tagstrip(tag,list):
    """tagstrip tag in list returns only what's between <tag> and </tag> 
on a line"""
    if (attrib(tag,list) is None):
        pre = re.compile('^.*<' + tag + '>')
        qre = re.compile('</' + tag +'>.*$')
        out = grep('<' + tag + '>',list)
        return qre.sub('',pre.sub('', out[0]))
    else:
        pre = re.compile('^.*<' + tag + ' ' + '[a-z]*="[a-z]*"' + '>')
        qre = re.compile('</' + tag +'>.*$')
        out = grep(tag,list)
        return qre.sub('',pre.sub('', out[0]))

def attrib(tag,list):
    """for <tag attrib="value"> will return 'attrib="value"'"""
    if grep('=',list):
        pre = re.compile('^.*<' + tag + ' ')
        qre = re.compile('>.*$')
        out = grep('<' + tag + ' ',list)
        if out:
            return qre.sub('',pre.sub('', out[0])).rstrip('\n')
    else:
        return None

def value(attrib,list):
    """for attrib="value" will return 'value'"""
    if grep('=',list):
        pat = re.compile(attrib + '="')
        qat = re.compile('".*$')
        out = pat.search(list)
        if out:
            return qat.sub('',pat.sub('', list))
    else:
        return None

def getpos(axis,list):
    """return a list containing position information of the given axis 
in the list"""
    axis_start = grepn(axis,data)
    for i in range(axis_start, len(data)):
        if (re.search('</positions>',data[i]) is not None):
            axis_end = i+1
            break
            
    return data[axis_start:axis_end]

# command line options

parser = argparse.ArgumentParser(description="Processes .xrdml files")
parser.add_argument("-f", "--force",
                  action="store_true", dest="force", default="False",
                  help="force reprocessing of the .xrdml file even if an up-to-date .ps or .ps.gz or already exists")
parser.add_argument("-n", "--noclobber",
                  action="store_true", dest="noclobber", default="False",
                  help="don't reprocesses the .xrdml file if the .ps or .ps.gz file already exists even if it is older")
parser.add_argument("-v", "--verbose",
                  action="store_true", dest="verbose", default="False",
                  help="list the .xrdml files as they are processed")
parser.add_argument("--w2t",
                  action="store_true", dest="omegatwotheta", default="False",
                  help="plot omega-2theta scans with theta as the y-axis")
parser.add_argument("--2tw",
                  action="store_true", dest="twothetaomega", default="False",
                  help="plot omega-2theta scans with 2theta as the y-axis")
parser.add_argument("-t", "-tau", "--tilt", nargs=1,
                  dest="taucorr", default=[0.0000], type=float,
                  help="tilt (in degrees) for correction of q values")
parser.add_argument("-o", "-om", "--omega", nargs=1,
                  dest="omcorr", default=[0.0000], type=float,
                  help="omega (in degrees) for correction of q values, requires twotheta to also be specified")
parser.add_argument("-q", "-tt", "--twotheta", nargs=1,
                  dest="ttcorr", default=[0.0000], type=float,
                  help="twotheta (in degrees) for correction of q values, requires omega to also be specified")
parser.add_argument('xrdfile', help=".xrdml filename")

# default is to reprocess the .xrdml file if it is newer than the .ps or .ps.gz file
# if both options are specified, -n wins

args = parser.parse_args()

xrdfile = args.xrdfile

tilt = args.taucorr[0]
if args.omcorr[0] != 0.0000:
    if args.ttcorr[0] != 0.0000:
        tilt = args.omcorr[0]-0.5*args.ttcorr[0]

if re.search('.xrdml',xrdfile) is None:
    print('Not an .xrdml file')
    sys.exit(2)

#don't process a file if the name says "fail"

if re.search('fail',xrdfile) is not None:
    print('Not processing ' + xrdfile)
    sys.exit(2)

#define output filenames
tilt_string=('')
if tilt != 0.0000:
    tilt_string=('_tc')

# print(tilt_string)

xrdmlre = re.compile('.xrdml')
datfile = xrdmlre.sub(tilt_string+'.dat', xrdfile)
gpfile = xrdmlre.sub(tilt_string+'.plt', xrdfile)
psfile = xrdmlre.sub(tilt_string+'.ps', xrdfile)
psgzfile = xrdmlre.sub(tilt_string+'.ps.gz', xrdfile)
#jpegfile = xrdmlre.sub(tilt_string+'.jpg', xrdfile)
#sjpegfile = xrdmlre.sub(tilt_string+'_s.jpg', xrdfile)

# -n
if args.noclobber is True:
    # print('-n set')
    # don't process .xrdml file if .ps or .ps.gz exists at all
    if os.path.isfile(psgzfile):
        print('Not reprocessing ' + xrdfile + ' since ' + psgzfile + ' already exists')
        sys.exit(0)
    elif os.path.isfile(psfile):
        print('Not reprocessing ' + xrdfile + ' since ' + psfile + ' already exists')
        sys.exit(0)

else:
    # -f
    if args.force is True:
        # print '-n not set, -f set'
        # reprocesses anyway
        print('Reprocessing ' + xrdfile)
        
    else:
        # print('-n not set, -f not set')
        # don't process .xrdml file is .ps.gz is newer
        if os.path.isfile(psgzfile):
            if os.path.getmtime(xrdfile) < os.path.getmtime(psgzfile):
                print('Not reprocessing ' + xrdfile + ' since ' + psgzfile + ' is already up-to-date')
                sys.exit(0)
        elif os.path.isfile(psfile):
            if os.path.getmtime(xrdfile) < os.path.getmtime(psfile):
                print('Not reprocessing ' + xrdfile + ' since ' + psfile + ' is already up-to-date')
                sys.exit(0)

# -v
if args.verbose is True:
    print('Processing ' + xrdfile)
    if tilt != 0.0000:
        print('Correcting tilt of ' + '{0:.4f}'.format(tilt) + '°')

wavelength = 0.1540562 #nm
k_len = 1/wavelength
aSi = 0.543102088 #nm

with open(xrdfile, 'r') as f:
    xrd = f.readlines()
    with open(datfile, 'w') as datf:
        datf.write ('# filename: ' + xrdfile + '\n')
        
        # what reflection is this?

        hkl = grepxml('reflection',xrd)
        hkl = grepxml('hkl',hkl)
        h = tagstrip('h',hkl)
        k = tagstrip('k',hkl)
        l = tagstrip('l',hkl)
        hkl = h[0], k[0], l[0]

        datf.write ('# reflection: ' + str(hkl) + '\n')

        scan_type = grep('measurementType',xrd)
        if grep('Scan',scan_type):
            ylow = ''
            yhigh = ''
            datf.write ('# type: Scan\n')
            data = grepxml('dataPoints',xrd)
            count_units = attrib('commonCountingTime',data)
            if (grep('"seconds"', count_units) is None):
                print('# Unrecognised counting time units')
                print('# ' + count_units)
                sys.exit(2)

            count_time = tagstrip('commonCountingTime',data)
            intensities_units = attrib('intensities',data)
            if (grep('"counts"', intensities_units) is None):
                print('# Unrecognised intensities units')
                print('# ' + intensities_units)
                sys.exit(2)
                
            counts = tagstrip('intensities',data)
            # atten = tagstrip('beamAttenuationFactors',data)
            # counts seem to already be multiplied by the attenuation factor
            counts = counts.split(" ")
            # atten = atten.split(" ")
            no_points = float(len(counts)-1)
            twotheta = getpos('2Theta',data)
            omega = getpos('Omega',data)
            phi = getpos('Phi',data)
            chi = getpos('Chi',data)
            x = getpos('X',data)
            y = getpos('Y',data)
            z = getpos('Z',data)
            twotheta = twotheta[1:-1]
            omega = omega[1:-1]
            phi = phi[1:-1]
            chi = chi[1:-1]
            xunits = '°'
            xaxis = '1'
            xlabel = '{/Symbol w}'
            if grep('start',twotheta):
                scan_type = '2Theta'
                twotheta_start = float(tagstrip('startPosition',twotheta))
                twotheta_end = float(tagstrip('endPosition',twotheta))
                xlabel = '2{/Symbol q}'
                xaxis = '2'
                xstart = str(twotheta_start)
                xend = str(twotheta_end)
            
            if grep('common',twotheta):
                twotheta_start = float(tagstrip('commonPosition',twotheta))
                twotheta_end = twotheta_start

            if grep('list',twotheta):
                scan_type = 'q'
                twotheta_list = tagstrip('listPositions',twotheta)
                twotheta_list = twotheta_list.split(" ")

            if grep('start',omega):
                omega_start = float(tagstrip('startPosition',omega))
                omega_end = float(tagstrip('endPosition',omega))
                xaxis = '1'
                xstart = str(omega_start)
                xend = str(omega_end)
                if (scan_type == '2Theta'):
                    scan_type = 'Omega_2Theta'
                    xlabel = '{/Symbol w}-2{/Symbol q}'
                    if args.omegatwotheta is True:
                        xlabel = '{/Symbol q}'
                        xaxis = '($2/2)'
                        xstart = str(twotheta_start/2)
                        xend = str(twotheta_end/2)
                        
                    if args.twothetaomega is True:
                        xlabel = '2{/Symbol q}'
                        xaxis = '2'
                        xstart = str(twotheta_start)
                        xend = str(twotheta_end)

                else:
                    scan_type = 'Omega'
                    xlabel = '{/Symbol w}'

            if grep('common',omega):
                omega_start = float(tagstrip('commonPosition',omega))
                omega_end = omega_start

            if grep('list',omega):
                scan_type = 'q'
                omega_list = tagstrip('listPositions',omega)
                omega_list = omega_list.split(" ")
                xlabel = '|q|'
                xunits = 'nm^{-1}'
                xaxis = '($8**2+$9**2)**0.5'

            if grep('start',phi):
                scan_type = 'Phi'
                phi_start = float(tagstrip('startPosition',phi))
                phi_end = float(tagstrip('endPosition',phi))
                xlabel = '{/Symbol f}'
                xaxis = '3'
                xstart = str(phi_start)
                xend = str(phi_end)

            if grep('common',phi):
                phi_start = float(tagstrip('commonPosition',phi))
                phi_end = phi_start

            if grep('start',chi):
                scan_type = 'Chi'
                chi_start = float(tagstrip('startPosition',chi))
                chi_end = float(tagstrip('endPosition',chi))
                xlabel = '{/Symbol c}'
                xaxis = '4'
                xstart = str(chi_start)
                xend = str(chi_end)
        
            if grep('common',chi):
                chi_start = float(tagstrip('commonPosition',chi))
                chi_end = chi_start

            if grep('start',x):
                scan_type = 'X'
                x_start = float(tagstrip('startPosition',x))
                x_end = float(tagstrip('endPosition',x))
                xlabel = 'x'
                xunits = 'mm'
                xaxis = '5'
                xstart = str(x_start)
                xend = str(x_end)
        
            if grep('common',x):
                x_start = float(tagstrip('commonPosition',x))
                x_end = x_start

            if grep('start',y):
                scan_type = 'Y'
                y_start = float(tagstrip('startPosition',y))
                y_end = float(tagstrip('endPosition',y))
                xlabel = 'y'
                xunits = 'mm'
                xaxis = '6'
                xstart = str(y_start)
                xend = str(y_end)
        
            if grep('common',y):
                y_start = float(tagstrip('commonPosition',y))
                y_end = y_start
                
            if grep('start',z):
                scan_type = 'Z'
                z_start = float(tagstrip('startPosition',z))
                z_end = float(tagstrip('endPosition',z))
                xlabel = 'z'
                xunits = 'mm'
                xaxis = '7'
                xstart = str(z_start)
                xend = str(z_end)
        
            if grep('common',z):
                z_start = float(tagstrip('commonPosition',z))
                z_end = z_start
                
            datf.write ('# x-axis: ' + scan_type + '\n')
            if tilt == 0.0000:
                datf.write ('# \n')

            else:
                datf.write ('# Correcting tilt of ' + '{0:.4f}'.format(tilt) + '° \n')

            datf.write ('# Omega\t2Theta\tPhi  \tChi   \tX    \tY    \tZ    \tq_para\tq_perp\tIntensity\n')
            datf.write ('# [deg]\t[deg] \t[deg]\t[deg] \t[mm] \t[mm] \t[mm] \t[nm-1]\t[nm-1]\t[counts/s]\n')
                   
            for index,count in enumerate(counts):
                if (scan_type == 'q'):
                    twotheta_val = float(twotheta_list[index])
                    omega_val = float(omega_list[index])
                  
                else:
                    omega_val = omega_start + (omega_end - omega_start)*index/no_points
                    twotheta_val = twotheta_start + (twotheta_end - twotheta_start)*index/no_points
                    
                phi_val = phi_start + (phi_end - phi_start)*index/no_points
                chi_val = chi_start + (chi_end - chi_start)*index/no_points
                x_val = x_start + (x_end - x_start)*index/no_points
                y_val = y_start + (y_end - y_start)*index/no_points
                z_val = z_start + (z_end - z_start)*index/no_points
                cps = float(count) / float(count_time)
                theta_rad = math.pi * twotheta_val / 360.0
                omega_rad = math.pi * omega_val / 180.0
                tilt_rad = math.pi * tilt / 180.0
                tau_rad = omega_rad - theta_rad - tilt_rad
                q_len = 2*k_len*math.sin(theta_rad)
                if (scan_type == 'q'):
                    if (index == 0):
                        xstart = '{0:.5f}'.format(q_len)
                                        
                q_para = '{0:.5f}'.format(q_len*math.sin(-tau_rad))
                q_perp = '{0:.5f}'.format(q_len*math.cos(-tau_rad))
                #out_line = omega_val, '\t', twotheta_val, '\t', phi_val, \
                #'\t', chi_val, '\t', x_val, '\t', y_val, '\t', z_val, \
                #'\t', q_para, '\t', q_perp, '\t', cps
                #sout_line = str(out_line)
                sout_line = '{0:.4f}'.format(omega_val) + '\t' + '{0:.4f}'.format(twotheta_val) + '\t'\
                + '{0:.2f}'.format(phi_val) + '\t' + '{0:.2f}'.format(chi_val) + '\t'\
                + '{0:.2f}'.format(x_val) + '\t' + '{0:.2f}'.format(y_val) + '\t' + '{0:.3f}'.format(z_val) + '\t'\
                + str(q_para) + '\t' + str(q_perp) + '\t'\
                + '{0:.1f}'.format(cps) + '\n'
                datf.write(sout_line)

            if (scan_type == 'q'):
                xend = '{0:.5f}'.format(q_len)
                # print(xstart + ' ' + xend)

            with open(gpfile, 'w') as gpf:
                gpf.write("""set encoding utf8
set xtics out
set mxtics 2
set logscale y
set ylabel "Intensity [ counts/s ]" offset -0.5,0
set format y "10^{%L}"
""")
                gpf.write('set yrange [' + ylow + ':' + yhigh + ']\n')
                gpf.write("""set encoding utf8
set style data lines
set linetype 1 lw 2 lc rgbcolor "red"
set linetype 2 lw 2 lc rgbcolor "blue"
set linetype 3 lw 2 lc rgbcolor "green"
set term postscript colour enhanced size 15cm,12cm
set palette rgb 30,31,32
unset key
""")
                underre = re.compile('_')
                filename = os.path.basename(xrdfile)
                title = underre.sub(' ', filename)
                gpf.write('set xrange [' + xstart + ':' + xend + ']\n')
                if (float(xend)-float(xstart) > 1.0):
                    if (float(xend)-float(xstart) < 4.0):
                        gpf.write('set xtics 0.5\n')
                        gpf.write('set format x "%.1f"\n')
                
                gpf.write('set title "' + title + '"\n')
                gpf.write('set out "' + psfile + '"\n')
                gpf.write('set xlabel "' + xlabel + ' [ ' + xunits + ' ]"\n')
                gpf.write('plot "' + datfile + '" using ' + xaxis + ':10\n')

            gpf.close()
        
        elif grep('Area measurement',scan_type):
            datf.write ('# type: Area measurement\n')
            y_scan_type = 'Omega'
            yaxis = '1'
            ylabel = '{/Symbol w}'
            yunits = '°'
            ycentre = grepxml('measurementStepAxisCenter',xrd)
            y_scan_type = attrib('position',ycentre)
            y_scan_type = value('axis',y_scan_type)
            if (y_scan_type == 'Phi'):
                yaxis = '3'
                ylabel = '{/Symbol f}'
                
            elif (y_scan_type == 'Chi'):
                yaxis = '4'
                ylabel = '{/Symbol c}'
                
            elif (y_scan_type == 'X'):
                yaxis = '5'
                ylabel = 'x'
                yunits = 'mm'
                
            elif (y_scan_type == 'Y'):
                yaxis = '6'
                ylabel = 'y'
                yunits = 'mm'
                
            elif (y_scan_type == 'Z'):
                yaxis = '7'
                ylabel = 'z'
                yunits = 'mm'
                
            data = grepxml('dataPoints',xrd)
            # assume that the first scan is representative
            count_units = attrib('commonCountingTime',data)
            if (grep('"seconds"', count_units) is None):
                print('# Unrecognised counting time units')
                print('# ' + count_units)
                sys.exit(2)
        
            count_time = tagstrip('commonCountingTime',data)
            intensities_units = attrib('intensities',data)
            if (grep('"counts"', intensities_units) is None):
                print('# Unrecognised intensities units')
                print('# ' + intensities_units)
                sys.exit(2)
                
            twotheta = getpos('2Theta',data)
            omega = getpos('Omega',data)
            phi = getpos('Phi',data)
            chi = getpos('Chi',data)
            x = getpos('X',data)
            y = getpos('Y',data)
            z = getpos('Z',data)
            twotheta = twotheta[1:-1]
            omega = omega[1:-1]
            phi = phi[1:-1]
            chi = chi[1:-1]
            xunits = '°'
            xaxis = '1'
            xlabel = '{/Symbol w}'
            if grep('start',twotheta):
                scan_type = '2Theta'

            if grep('start',omega):
                if (scan_type == '2Theta'):
                    scan_type = 'Omega_2Theta'

                else:
                    scan_type = 'Omega'

            if grep('start',phi):
                scan_type = 'Phi'

            if grep('start',chi):
                scan_type = 'Chi'
            
            if grep('start',x):
                scan_type = 'X'
            
            if grep('start',y):
                scan_type = 'Y'
            
            if grep('start',z):
                scan_type = 'Z'

            datf.write ('# x-axis: ' + scan_type + '\n')
            datf.write ('# y-axis: ' + y_scan_type + '\n')
            if tilt == 0.0000:
                datf.write ('# \n')

            else:
                datf.write ('# Correcting tilt of ' + '{0:.4f}'.format(tilt) + '° \n')

            datf.write ('# Omega\t2Theta\tPhi  \tChi   \tX    \tY    \tZ    \tq_para\tq_perp\tIntensity\n')
            datf.write ('# [deg]\t[deg] \t[deg]\t[deg] \t[mm] \t[mm] \t[mm] \t[nm-1]\t[nm-1]\t[counts/s]\n')
            scan_lines = grepalln('<dataPoints>',xrd)
            for this_scan in scan_lines:
                data = grepxml('dataPoints',xrd[this_scan:])
                counts = tagstrip('intensities',data)
                counts = counts.split(" ")
                no_points = float(len(counts)-1)
                twotheta = getpos('2Theta',data)
                omega = getpos('Omega',data)
                phi = getpos('Phi',data)
                chi = getpos('Chi',data)
                x = getpos('X',data)
                y = getpos('Y',data)
                z = getpos('Z',data)
                twotheta = twotheta[1:-1]
                omega = omega[1:-1]
                phi = phi[1:-1]
                chi = chi[1:-1]
                xformat = '#'
                xtics = '# '
                yformat = '#'
                if grep('start',twotheta):
                    scan_type = '2Theta'
                    twotheta_start = float(tagstrip('startPosition',twotheta))
                    twotheta_end = float(tagstrip('endPosition',twotheta))
                    xlabel = '2{/Symbol q}'
                    xaxis = '2'
                    xstart = str(twotheta_start)
                    xend = str(twotheta_end)
                
                if grep('common',twotheta):
                    twotheta_start = float(tagstrip('commonPosition',twotheta))
                    twotheta_end = twotheta_start

                if grep('start',omega):
                    omega_start = float(tagstrip('startPosition',omega))
                    omega_end = float(tagstrip('endPosition',omega))
                    if (scan_type == '2Theta'):
                        scan_type = 'Omega_2Theta'
                        xlabel = '{/Symbol w}-2{/Symbol q}'
                        if (y_scan_type == 'Omega'):
                            # reciprocal space map
                            xaxis = '8'
                            xlabel = 'q_{||}'
                            xunits = 'nm^{-1}'
                            xstart = ''
                            xend = ''
                            xformat = 'set format x "%.1f"'
                            if (h[0] == '0'):
                                # symmetrical
                                xtics = 'set xtics 0.1'
                                
                            yaxis = '9'
                            ylabel = r'q_{/Symbol \\^}' # use a raw string to avoid that "A backslash-character pair that is not a valid escape sequence now generates a SyntaxWarning, instead of DeprecationWarning."
                            yunits = 'nm^{-1}'
                            yformat = 'set format y "%.1f"'
                        
                    else:
                        scan_type = 'Omega'
                        xlabel = '{/Symbol w}'                  

                if grep('common',omega):
                    omega_start = float(tagstrip('commonPosition',omega))
                    omega_end = omega_start

                if grep('start',phi):
                    scan_type = 'Phi'
                    phi_start = float(tagstrip('startPosition',phi))
                    phi_end = float(tagstrip('endPosition',phi))
                    xlabel = '{/Symbol f}'
                    xaxis = '3'
                    xstart = str(phi_start)
                    xend = str(phi_end)

                if grep('common',phi):
                    phi_start = float(tagstrip('commonPosition',phi))
                    phi_end = phi_start

                if grep('start',chi):
                    scan_type = 'Chi'
                    chi_start = float(tagstrip('startPosition',chi))
                    chi_end = float(tagstrip('endPosition',chi))
                    xlabel = '{/Symbol c}'
                    xaxis = '4'
                    xstart = str(chi_start)
                    xend = str(chi_end)
            
                if grep('common',chi):
                    chi_start = float(tagstrip('commonPosition',chi))
                    chi_end = chi_start

                if grep('start',x):
                    scan_type = 'X'
                    x_start = float(tagstrip('startPosition',x))
                    x_end = float(tagstrip('endPosition',x))
                    xlabel = 'x'
                    xunits = 'mm'
                    xaxis = '5'
                    xstart = str(x_start)
                    xend = str(x_end)
            
                if grep('common',x):
                    x_start = float(tagstrip('commonPosition',x))
                    x_end = x_start

                if grep('start',y):
                    scan_type = 'Y'
                    y_start = float(tagstrip('startPosition',y))
                    y_end = float(tagstrip('endPosition',y))
                    xlabel = 'y'
                    xunits = 'mm'
                    xaxis = '6'
                    xstart = str(y_start)
                    xend = str(y_end)
            
                if grep('common',y):
                    y_start = float(tagstrip('commonPosition',y))
                    y_end = y_start
                    
                if grep('start',z):
                    scan_type = 'Z'
                    z_start = float(tagstrip('startPosition',z))
                    z_end = float(tagstrip('endPosition',z))
                    xlabel = 'z'
                    xunits = 'mm'
                    xaxis = '7'
                    xstart = str(z_start)
                    xend = str(z_end)
            
                if grep('common',z):
                    z_start = float(tagstrip('commonPosition',z))
                    z_end = z_start

                for index,count in enumerate(counts):
                    omega_val = omega_start + (omega_end - omega_start)*index/no_points
                    twotheta_val = twotheta_start + (twotheta_end - twotheta_start)*index/no_points
                    phi_val = phi_start + (phi_end - phi_start)*index/no_points
                    chi_val = chi_start + (chi_end - chi_start)*index/no_points
                    x_val = x_start + (x_end - x_start)*index/no_points
                    y_val = y_start + (y_end - y_start)*index/no_points
                    z_val = z_start + (z_end - z_start)*index/no_points
                    cps = float(count) / float(count_time)
                    theta_rad = math.pi * twotheta_val / 360.0
                    omega_rad = math.pi * omega_val / 180.0
                    tilt_rad = math.pi * tilt / 180.0
                    tau_rad = omega_rad - theta_rad - tilt_rad
                    q_len = 2*k_len*math.sin(theta_rad)
                    q_para = '{0:.5f}'.format(q_len*math.sin(-tau_rad))
                    q_perp = '{0:.5f}'.format(q_len*math.cos(-tau_rad))
                    #out_line = omega_val, '\t', twotheta_val, '\t', phi_val, \
                    #'\t', chi_val, '\t', x_val, '\t', y_val, '\t', z_val, \
                    #'\t', q_para, '\t', q_perp, '\t', cps
                    #sout_line = str(out_line)
                    sout_line = '{0:.4f}'.format(omega_val) + '\t' + '{0:.4f}'.format(twotheta_val) + '\t'\
                    + '{0:.2f}'.format(phi_val) + '\t' + '{0:.2f}'.format(chi_val) + '\t'\
                    + '{0:.2f}'.format(x_val) + '\t' + '{0:.2f}'.format(y_val) + '\t' + '{0:.3f}'.format(z_val) + '\t'\
                    + str(q_para) + '\t' + str(q_perp) + '\t'\
                    + '{0:.1f}'.format(cps) + '\n'
                    datf.write(sout_line)
                    
                datf.write('\n')
                
            with open(gpfile, 'w') as gpf:
                gpf.write("""set encoding utf8
set pm3d map
set size square
set xtics out
set ytics out
set mxtics 2
set mytics 2
set logscale z
set log cb
set format cb "10^{%L}"
set zrange [1:]
set cbrange [1:]
set autoscale fix
set term postscript colour enhanced size 20cm,16cm
# set palette rgb 7,5,15
# set palette rgb 30,31,32
set palette @MATLAB
# palette MATLAB macro must be defined in ~/.gnuplot
unset key
""")
                underre = re.compile('_')
                filename = os.path.basename(xrdfile)
                title = underre.sub(' ', filename)
                # gpf.write('set xrange [' + xstart + ':' + xend + ']\n')
                #if (float(xend)-float(xstart) > 1.0):
                #    if (float(xend)-float(xstart) < 4.0):
                #        gpf.write('set xtics 0.5\n')
                #        gpf.write('set format x "%.1f"\n')
                gpf.write(xformat + '\n')
                gpf.write(xtics + '\n')
                gpf.write(yformat + '\n')
                gpf.write('set title "' + title + '"\n')
                gpf.write('set out "' + psfile + '"\n')
                gpf.write('set xlabel "' + xlabel + ' [ ' + xunits + ' ]"\n')
                gpf.write('set ylabel "' + ylabel + ' [ ' + yunits + ' ]" offset -1.0,0.0\n')
                gpf.write('splot "' + datfile + '" using ' + xaxis + ':' + yaxis + ':10\n')

            gpf.close()
                
        else:
            print('Unrecognised scan type')
            print(scan_type)
            sys.exit(2)

    datf.close()

f.close()

call(['gnuplot', gpfile]) # if gnuplot was added to PATH this will also work on Windows
#call(['rm', gpfile])
#
# requires Ghostscript and ImageMagick:
call(['xrd-png.bat', psfile]) # output is a good-quality image, Windows version
#call(['xrd-png', psfile]) # output is a good-quality image, Linux and OSX version
#call(['xrd-post', psfile]) # output is compatible with lepecvd and lg2 databases
