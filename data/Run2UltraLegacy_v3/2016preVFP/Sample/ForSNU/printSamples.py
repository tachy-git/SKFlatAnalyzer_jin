import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--baseDir", required=True, type=str, help="base directory")
args = parser.parse_args()

for f in os.listdir(args.baseDir):
    print(f"{args.baseDir}/{f}")
