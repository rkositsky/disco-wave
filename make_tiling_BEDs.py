# Goal: Make windows of 20-bp in captured FISH regions

# Step 1: Cut FISH regions into x-bp tiling regions
# Step 2: Intersect with captured probes to get captured windows. Keep complete windows?

# Rachel Kositsky
# Created: 2020-01-16
# Updated: 2020-01-28

import argparse
import os
import subprocess
import sys


def parse_args(args=None):
	"""Parse command line arguments and return constants
	Get input BAM and output directory"""
	parser = argparse.ArgumentParser(
		description="Create tiling windows for genes listed in tiling file")

	parser.add_argument("input_dir",
		help="Input directory with BED files per each region in which to call" \
		" translocations.", type=str)

	parser.add_argument("tmp_dir",
		help="Output directory for temporary files", type=str)

	parser.add_argument("output_file",
		help="Output path of final BED file with tiled regions of each gene " \
		"region given in input_dir", type=str)

	parser.add_argument("-s", "--tiling_size", default=20,
		help="Input size of tiled windows. Default is 20bp. Smaller windows " \
		"provide greater resoultion at the cost of computation time.",
		type=int)

	parser.add_argument("-p", "--panel_bed",
		help="Input BED file of capture panel regions. Do not specify for " \
		"whole genome sequencing.", type=str)

	results = parser.parse_args(args)

	return results


def create_windows(in_bed, out_bed, tiling_size):
	"""Prints tiling window to bed file"""

	# Read in informations
	with open(in_bed, "r") as in_f:
		bed_line = in_f.readline().strip().split("\t")
		chrom = bed_line[0]
		in_start = int(bed_line[1])
		in_end = int(bed_line[2])
		annot = bed_line[3]

	# Write out tiled regions
	with open(out_bed, "w") as out_f:
		for start in range(in_start, in_end, tiling_size):
			out_f.write("{0}\t{1}\t{2}\t{3}\n".format(chrom, 
				start, start + tiling_size, annot))


def intersect_bed_with_baits(tiled_bed, panel_bed, output_bed, tmp_dir):
	"""Intersects windows with captured bed file, then merges duplicate lines"""
	
	# Instead of messing with subprocess.PIPE and closing processes, just
	# write to a temporary file and run two separate commands
	tmp_bed = os.path.join(tmp_dir, "intersect.tmp.bed")
	with open(tmp_bed, "w") as out_f:
		subprocess.call(["bedtools", "intersect", "-a", tiled_bed, 
			"-b", panel_bed, "-wa"], stdout=out_f)

	with open(output_bed, "w") as out_f:
		subprocess.call(["uniq", tmp_bed], stdout=out_f)

	os.remove(tmp_bed)


def main(input_dir, tmp_dir, output_file, tiling_size, panel_bed):

	# Make temporary output folder
	if not os.path.exists(tmp_dir):
		os.mkdir(tmp_dir)

	# Create captured tiled BEDs and merge all together
	with open(output_file, "a") as out_f:
		input_files = [f for f in os.listdir(input_dir)
					   if os.path.isfile(os.path.join(input_dir, f))]
		if input_files == []:
			raise Exception(f"No input files found for tiling in {input_dir}")

		for in_file in input_files:
			full_in_file = os.path.join(input_dir, in_file)
			tmp_out1 = os.path.join(tmp_dir, f"{in_file}.tiled.all.bed")
			tmp_out2 = os.path.join(tmp_dir, f"{in_file}.tiled.panel.bed")

			create_windows(full_in_file, tmp_out1, tiling_size)
			if panel_bed == None:
				tmp_out2 = tmp_out1
			else:
				intersect_bed_with_baits(tmp_out1, panel_bed, tmp_out2, tmp_dir)

			# Append tiled BED to final merged version
			with open(tmp_out2, "r") as in_f:
				out_f.write(in_f.read())

	# Make unique version
	tmp_bed = os.path.join(tmp_dir, "uniq.merged.bed")

	# Copy intersect_file to tmp_bed
	subprocess.call(["cp", output_file, tmp_bed])

	#### TODO: add BED file sorting?? #######

	# Copy unique version of tmp_bed to output_file
	with open(output_file, "w") as out_f:
		subprocess.call(["uniq", tmp_bed], stdout=out_f)

	os.remove(tmp_bed)


if __name__ == "__main__":
	args = parse_args(sys.argv[1:])
	main(args.input_dir, args.tmp_dir, args.output_file, args.tiling_size, 
		args.panel_bed)

