import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--mother", "-m", required=True, type=str, help="mother directory")
parser.add_argument("--output_filename", "-o", required=True, type=str, help="output filename")
args = parser.parse_args()
mother = args.mother
output_filename = args.output_filename

main_dirs = os.listdir(mother)
with open(f"{output_filename}.txt", "w" ) as f:
		for d in main_dirs:
				sub_dirs = os.listdir(f"{mother}/{d}")
				for s in sub_dirs:
						f.write(f"{mother}/{d}/{s}\n")

