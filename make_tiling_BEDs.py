# Goal: Make windows of 20-bp in captured FISH regions

# Step 1: Cut FISH regions into x-bp tiling regions
# Step 2: Intersect with captured probes to get captured windows. Keep complete windows?

# Rachel Kositsky
# Created: 2020-01-16
# Updated: 2020-01-28

import argparse
import os
import subprocess


def parse_args(args=None):
	"""Parse command line arguments and return constants
	Get input BAM and output directory"""
	parser = argparse.ArgumentParser(
		description="Create tiling windows for genes listed in tiling file")

	parser.add_argument("input_dir",
		help="Input directory with BED files per each region in which to call" \
		" translocations.")

	parser.add_argument("tmp_dir",
		help="Output directory for temporary files")

	parser.add_argument("output_file",
		help="Output path of final BED file with tiled regions of each gene " \
		"region given in input_dir")

	parser.add_argument("tiling_size", default = 20,
		help="Input size of tiled windows. Default is 20bp. Smaller windows " \
		"provide greater resoultion at the cost of computation time.")

	parser.add_argument("panel_bed", default = None,
		help="Input BED file of capture panel regions. Do not specify for " \
		"whole genome sequencing.")

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


def intersect_bed_with_baits(tiled_bed, bait_bed, output_bed):
	"""Intersects windows with captured bed file, then merges duplicate lines"""
	
	# Instead of messing with subprocess.PIPE and closing processes, just
	# write to a temporary file and run two separate commands
	tmp_bed = "intersect.tmp.bed"
	with open(tmp_bed, "w") as out_f:
		subprocess.call(["bedtools", "intersect", "-a", tiled_bed, 
			"-b", bait_bed, "-wa"], stdout=out_f)

	with open(output_bed, "w") as out_f:
		subprocess.call(["uniq", tmp_bed], stdout=out_f)

	os.remove(tmp_bed)

def main(tiling_size, bait_bed, in_dir, tmp_dir, out_file):

	#bait_bed = os.path.expanduser("~/Dropbox/Rachel/Projects/Gene_panel/Twist_design/2019-10-21/annotation_beds/Twist_8MB_panel_with_ERCCs.gencode.sorted.bed")
	bait_bed = os.path.expanduser("~/Dropbox/Rachel/Projects/Gene_panel/Twist_design/2019-10-21/annotation_beds/Twist_8MB_panel_with_ERCCs.non-repetitive.bed")
	in_dir = os.path.expanduser("~/Dropbox/Rachel/Projects/WHO/Translocations/annotation/Vysis_FISH")
	#out_dir = os.path.expanduser("~/Dropbox/Rachel/Projects/WHO/Translocations/bed/tiling_windows/{0}bp_tiling".format(tiling_size))
	out_dir = os.path.expanduser("~/Dropbox/Rachel/Projects/WHO/Translocations/bed/tiling_windows/{0}bp_tiling_no_repeats".format(tiling_size))

	if not os.path.exists(out_dir):
		os.mkdir(out_dir)

	# Create captured tiled BEDs and merge all three together
	final_bed = os.path.join(out_dir, "BCL6_MYC_BCL2_FISH.tiled.bed")
	with open(final_bed, "a") as out_f:
		for gene in ["BCL6", "MYC", "BCL2"]:
			in_file = os.path.join(in_dir, gene + "_FISH.bed")
			out_file = os.path.join(out_dir, gene + "_FISH.tiled.all.bed")
			intersect_file = os.path.join(out_dir, gene + "_FISH.tiled.bed")

			create_windows(in_file, out_file, tiling_size)
			intersect_bed_with_baits(out_file, bait_bed, intersect_file)

			# Append tiled BED to final merged version
			with open(intersect_file, "r") as in_f:
				out_f.write(in_f.read())

	# Make unique version
	tmp_bed = "uniq.merged.bed"

	# Copy intersect_file to tmp_bed
	subprocess.call(["cp", final_bed, tmp_bed])

	# Copy unique version of tmp_bed to final_bed
	with open(final_bed, "w") as out_f:
		subprocess.call(["uniq", tmp_bed], stdout=out_f)

	os.remove(tmp_bed)


if __name__ == "__main__":
	args = parse_args(sys.argv[1:])
	main(args.tiling_size, args.panel_bed, args.input_dir, args.output_dir)

