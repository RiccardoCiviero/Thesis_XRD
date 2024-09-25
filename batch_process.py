import os
import subprocess
import argparse

parser = argparse.ArgumentParser(description="Batch rotate then extract omega scan profile in all files in a path")
parser.add_argument("-s", "--start", nargs=1,
                  dest="startnum", default= 0, type=int,
                  help="which scan to start with")
parser.add_argument("-e", "--end", nargs=1,
                  dest="endnum", default= 200, type=int,
                  help="which scan to end with")
parser.add_argument('path', help="root path containing data")
args = parser.parse_args()

print("Executing xrd_rsm_omega.py in path {}!\n\n\n".format(args.path))
for root, dirs, files in os.walk(args.path):
    for name in files:
        if not "w-scans" in name:
            f = root + "/" + name
            subprocess.run(['python', 'xrd_rsm_omega.py', f])

print("Extracting omega profiles from {} to {} using xrd-rsm-project!".format(args.startnum, args.endnum))
for root, dirs, files in os.walk(args.path):
    for name in files:
        if "w-scans" in name:
            f = root + "/" + name
            print(f, args.startnum, args.endnum)
            subprocess.run(['python', 'xrd-rsm-project.py', f, '-s', str(args.startnum), '-e', str(args.endnum)])

print("\n\n\nIt was a smooth ride!")