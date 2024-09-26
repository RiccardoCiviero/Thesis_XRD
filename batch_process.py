import os
import subprocess
import argparse

parser = argparse.ArgumentParser(description="Batch rotate then extract omega scan profile in all files in a path. It will create a Omega and a Projected folder")
parser.add_argument("-s", "--start", nargs=1,
                  dest="startnum", default= [0], type=int,
                  help="which scan to start with")
parser.add_argument("-e", "--end", nargs=1,
                  dest="endnum", default= [200], type=int,
                  help="which scan to end with")
parser.add_argument('path', help="root path containing data")
args = parser.parse_args()

print("Executing xrd_rsm_omega.py in path {}!".format(args.path))
if not os.path.isdir(str(args.path) + "/Omega"):
    os.makedirs(str(args.path) + "/Omega")

for root, dirs, files in os.walk(args.path):
    for name in files:
        if not "w-scans" in name:
            f = root + "/" + name
            subprocess.run(['python', 'xrd_rsm_omega.py', f, "-d", str(args.path) + "/Omega"], stdout=open(os.devnull, 'wb'))

outpath = str(args.path) + "/Projected_{}-{}".format(str(args.startnum[0]), str(args.endnum[0]))
if not os.path.isdir(outpath):
    os.makedirs(outpath)

print("Extracting omega profiles from {} to {} using xrd-rsm-project!".format(args.startnum, args.endnum))
for root, dirs, files in os.walk(args.path + "/Omega"):
    for name in files:
        if "w-scans" in name:
            f = root + "/" + name
            subprocess.run(['python', 'xrd-rsm-project.py', f, '-s', str(args.startnum[0]), '-e', str(args.endnum[0]), '-d', outpath], stdout=open(os.devnull, 'wb'))

print("It was a smooth ride!")