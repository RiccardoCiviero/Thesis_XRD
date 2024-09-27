import os
import subprocess
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Batch rotate then extract omega scan profile in all files in a path. It will create a Omega and a Projected folder")
parser.add_argument('path', help="root path containing data")
args = parser.parse_args()

# print("Extracting omega2theta profiles using xrd-rsm-project and finding peaks!")
# for root, dirs, files in os.walk(args.path + "/"):
#     for name in files:
#             f = root + "/" + name
#             subprocess.run(['python', 'xrd-rsm-project.py', f, '-d', root + "Omega2Theta"], stdout=open(os.devnull, 'wb'))

peaklist = []
print("Extracting the peak position!")
for root, dirs, files in os.walk(args.path + "/Omega2Theta"):
    for name in files:
            f = root + "/" + name
            omega, ttheta, phi, chi, xpos, ypos, zpos, qx, qz, intens = np.genfromtxt(f, dtype=float, skip_header=7, unpack=True)
            peaklist.append(np.argmax(intens[0:omega.size()//2]))

print("And projecting on it!")
for root, dirs, files in os.walk(args.path + "/Omega"):
    for name in files:
        if "w-scans" in name:
            f = root + "/" + name
            subprocess.run(['python', 'xrd-rsm-project.py', f, '-s', str(args.startnum[0]), '-e', str(args.endnum[0]), '-d', args.path + ""], stdout=open(os.devnull, 'wb'))
