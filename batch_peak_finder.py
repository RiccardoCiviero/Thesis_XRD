import os
import subprocess
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Batch rotate then extract omega scan profile in all files in a path. It will create a Omega and a Projected folder")
parser.add_argument('path', help="root path containing data")
args = parser.parse_args()

print("Extracting omega2theta profiles using xrd-rsm-project and finding peaks!")
for root, dirs, files in os.walk(args.path + "/113xx"):
    for name in files:
            f = root + "/" + name
            subprocess.run(['python', 'xrd-rsm-project.py', f, '-d', args.path + "/Omega2Theta"], stdout=open(os.devnull, 'wb'))

for root, dirs, files in os.walk(args.path + "/113xx"):
    for name in files:
        f = root + "/" + name
        subprocess.run(['python', 'xrd_rsm_omega.py', f, "-d", args.path + "/Omega"], stdout=open(os.devnull, 'wb'))

peaklist = { "11352-annealed_004" : 0, "11352-as-grown_004" : 0, "11354-annealed_004" : 0, "11354-as-grown_004" : 0, "11355-annealed_004" : 0, "11355-as-grown_004" : 0, "11358-annealed_004" : 0, "11358-as-grown_004" : 0, "11371-annealed_004" : 0, "11371-as-grown_004" : 0 , "11352-annealed_224" : 0, "11352-as-grown_224" : 0, "11354-annealed_224" : 0, "11354-as-grown_224" : 0, "11355-annealed_224" : 0, "11355-as-grown_224" : 0, "11358-annealed_224" : 0, "11358-as-grown_224" : 0, "11371-annealed_224" : 0, "11371-as-grown_224" : 0 }

print("Extracting the peak position!")
for root, dirs, files in os.walk(args.path + "/Omega2Theta"):
    for name in files:
            f = root + "/" + name
            omega, ttheta, phi, chi, xpos, ypos, zpos, qx, qz, intens = np.genfromtxt(f, dtype=float, skip_header=7, unpack=True)
            n = name[27:26+19]
            print(n)
            peaklist[n] = (np.argmax(intens[0:np.size(ttheta)//2]))

i = 0
print("And projecting on it!")
for root, dirs, files in os.walk(args.path + "/Omega"):
    for name in files:
        f = root + "/" + name
        n = name[27:26+19]
        subprocess.run(['python', 'xrd-rsm-project.py', f, '-s', str(peaklist[n]), '-e', str(peaklist[n]), '-d', args.path + "/Projection_SiGePeak"], stdout=open(os.devnull, 'wb'))
        i = i+1

