#!/usr/bin/python3

# Get discordant reads and identify translocations using tiling windows.
#
# Rachel Kositsky
# Updated: 2020-06-17

import argparse
import os
import sort_BAM_file_by_mate_chromosome
import subprocess
import sys


def parse_args(args=None):
	"""Parse command line arguments"""
	parser = argparse.ArgumentParser(
		description="Find translocations using discordant reads from tiling windows "
		"from short-read sequencing aligned files")

	parser.add_argument("in_bam", help="Input BAM file")

	parser.add_argument("--out_translocation_table",
		help="File path to write tab-separated table of candidate translocations",
		default="candidate_translocations.tsv")

	parser.add_argument("--out_discordant_bam",
		help="File path to write subsetted BAM with only discordantly mapped reads",
		default="discordant_reads.diff_chrom.bam")
	
	parser.add_argument("--figure_output_dir",
		help="Directory where output figure will be written",
		default="figure_output_dir")

	parser.add_argument("--sample_name", default="",
		help="Sample name")

	parser.add_argument("--tiling_BED", type=str,
		help="Output from make_tiling_BEDs.py. BED file with tiled windows where "
		"translocations can be called. Default: MYC, BCL2, BCL6 FISH capture regions.",
		default="")

	parser.add_argument("--min_mapping_quality", default=0, type=int,
		help="Minimum read quality for read that get considered.")

	parser.add_argument("--merge_distance", default=0, type=int,
		help="Maximum distance between reads where they're still considered"
		" part of the same translocation. Default: 0 bp")

	parser.add_argument("--min_read_pairs", default=10, type=int,
		help="Minimum supporting reads for a translocation to be reported."
		" Default: 10 read pairs")

	parser.add_argument("--nr_cpus", default=1, type=int,
		help="Number of CPUs for parallelizing discordant read selection. "
		"Default: 1")

	# TODO: add option to define genome. Hardcoded for GENCODE right now.
	results = parser.parse_args(args)

	return results


def run_chrom_coverages(tiling_BED, input_dir, coverage_dir, diff_chrom_bam,
	resource_dir):
	"""Replaces run_chrom_coverages.bash"""

	# Set up some hard-coded paths
	if tiling_BED == "":
		tiling_BED = os.path.join(resource_dir, "BCL6_MYC_BCL2_FISH.tiled.bed")
	else:
		tiling_BED = tiling_BED
	genome_file = os.path.join(resource_dir, "hg38.sizes.gencode")

	if not os.path.exists(coverage_dir):
		os.mkdir(coverage_dir)

	# Run coverages for the total BAM
	out_bed = os.path.join(coverage_dir, "all.bed")
	with open(out_bed, "w") as out_f:
		subprocess.call(["bedtools", "coverage", "-sorted", "-g", genome_file,
			"-wa", "-a", tiling_BED, "-b", diff_chrom_bam], stdout=out_f)

	# Run it on each of the different mate chromosome BAMs
	chroms = ["chr" + i for i in list(map(str, range(1, 23))) + ["X", "Y", "M"]]

	for chrom in chroms:
		out_bed = os.path.join(coverage_dir, "{0}.bed".format(chrom))
		in_bam = os.path.join(input_dir, "{0}.bam".format(chrom))
		with open(out_bed, "w") as out_f:
			subprocess.call(["bedtools", "coverage", "-sorted", 
				"-g", genome_file, "-wa", "-a", tiling_BED, "-b", in_bam], 
				stdout=out_f)
		

def main(args):
	"""Find translocations using discordant reads from tiled windows from 
	short-read sequencing aligned files"""

	# Get folder of current executed file to get the repository path
	tool_dir = os.path.dirname(os.path.realpath(__file__))
	resource_dir = os.path.join(tool_dir, "resources")

	# Create temporary work directory
	work_dir = "tmp_dir"
	if not os.path.exists(work_dir):
		os.mkdir(work_dir)

	# 1. Extract discordant reads
	print("Extracting discordant reads...")
	extraction_script = os.path.join(tool_dir, "select_discordant_reads.bash")
	subprocess.call([extraction_script, args.in_bam, args.out_discordant_bam, 
		str(args.nr_cpus), str(args.min_mapping_quality)])

	# 2. Sort out the other chromosomes by where their mate maps to
	print("Sorting reads by their mate chromosome...")
	bam_by_chrom_dir = os.path.join(work_dir, "bams_by_chrom")
	sort_BAM_file_by_mate_chromosome.sort_BAM(args.out_discordant_bam, 
		bam_by_chrom_dir)

	# 3. Run bedtools coverage on tiled windows
	print("Getting tiling coverages...")
	tiling_coverage_dir = os.path.join(work_dir, "tiled_coverage_by_chrom")
	# TODO: include compatibility with other genome versions
	run_chrom_coverages(tiling_BED = args.tiling_BED, 
		input_dir = bam_by_chrom_dir, coverage_dir = tiling_coverage_dir,
		diff_chrom_bam = args.out_discordant_bam, resource_dir = resource_dir)

	# 4. Run R script to get the pretty pictures and the big output table
	print("Running tiling coverage parsing script and generating output table...")
	tiling_script = os.path.join(tool_dir, "tiling_discordant_reads.R")
	tiling_script_helper_functions = os.path.join(tool_dir, 
		"tiling_discordant_reads_helper_functions.R")
	
	subprocess.call(["Rscript", "--vanilla", tiling_script, 
		"--helper_functions", tiling_script_helper_functions, 
		"--input_directory", tiling_coverage_dir,
		"--output_table", args.out_translocation_table,
		"--figure_output_directory", args.figure_output_dir,
		"--sample_name", args.sample_name,
		"--min_reads", str(args.min_read_pairs),
		"--merge_distance", str(args.merge_distance)])

	print("Complete!")


if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))

