# Goal: Sort input BAM into different chromosomes by its mate
# Result: output directory has 25 BAM files, chr1.bam, chr2.bam, ..., chrX.bam,
# chrY.bam, chrM.bam. Each BAM file is the reads with MATES mapping to that
# chromosome.

import argparse
import os
import pysam
import sys


def parse_args(args=None):
	"""Parse command line arguments and return constants
	Get input BAM and output directory"""
	parser = argparse.ArgumentParser(
		description="Sort input BAM into different files by mate chromosome. " \
		"Intended for discordantly mapped input BAMs.")

	parser.add_argument("input_bam",
		help="Input BAM of discordantly mapped reads")

	parser.add_argument("output_dir",
		help="Directory for output files, e.g. my_dir/AAA343_A/")

	results = parser.parse_args(args)

	return (results.input_bam, results.output_dir)


def parse_UCSC_canonical_chromosome(ucsc_chrom):
	"""Used for extracting the canonical reference chromosome from a UCSC
	contig name. Meant to be temporary until a better conversion is found for 
	all GENCODE loci."""
	return(ucsc_chrom.split("_")[0])


def sort_BAM(input_bam, output_dir):
	# Chromosomes to sort into
	chroms = ["chr" + i for i in list(map(str, range(1, 23))) + ["X", "Y", "M"]]

	# Parse argument
	if not os.path.exists(output_dir):
		os.mkdir(output_dir)

	# Hardcoded/needs to come with the script:
	# GENCODE contig to chromosome mapping
	gencode_mapping_file = os.path.join(os.path.dirname(os.path.realpath(__file__)),
		"resources/GRCh38_gencode2UCSC.txt")

	# Open input BAM file for reading by pysam
	input_bam = pysam.AlignmentFile(input_bam, "rb")

	## Populate canonical chromosome conversion dictionary

	# chrom_dict: keys = all GENCODE contigs; values = canonical chromosomes or ""
	chrom_dict = {}
	# chrom_bamfile_dict: keys = chosen canonical chromosomes; values = python files
	chrom_bamfile_dict = {}
	
	with open(gencode_mapping_file, "r") as in_f:
		for line in in_f.readlines():
			gencode_contig, orig_chrom = line.strip("\n").split("\t")
			orig_chrom = parse_UCSC_canonical_chromosome(orig_chrom)

			chrom_dict[gencode_contig] = orig_chrom

			#if orig_chrom not in chroms:
				#pass
				#print("Warning: ignoring all reads with mate at {0}/{1}".format(
				#	gencode_contig, orig_chrom))
	
	# Parse in the rest of the header references - they're not all in the
	# GENCODE file
	for ref_name in input_bam.references:
		if ref_name not in chrom_dict.keys():
			chrom_dict[ref_name] = ""

	# Open all chromosome files and copy over the SAM header
	for chrom in chroms:
		chrom_file = os.path.join(output_dir, "{0}.bam".format(chrom))
		chrom_bamfile_dict[chrom] = pysam.AlignmentFile(chrom_file, "wb", 
			template=input_bam)


	# Read through file and sort out discordant reads
	n_reads = 0
	n_read_pairs_on_same_chrom_after_conversion = 0
	for read in input_bam.fetch():
		n_reads += 1

		# Convert from GENCODE contig to reference contig
		mate_chrom = chrom_dict[read.next_reference_name]
		read_chrom = chrom_dict[read.reference_name]

		# If, after converting, you find that it's the same reference chromsome
		# - skip it and move on. 
		# Here we assume that contigs and original chromsome would align well, 
		# which is not always true. Since we're using this for detecting
		# interchromsomal translocations, though, it's a useful assumption.
		if read_chrom == mate_chrom:
			n_read_pairs_on_same_chrom_after_conversion += 1
			continue

		# Append read to output BAM file if mate chromosome is one of our chosen chromosomes
		if mate_chrom in chroms:
			chrom_bamfile_dict[mate_chrom].write(read)

	# Print how many non-discordant reads you saw
	print("{0} reads of {1} reads parsed were aligned to the same primary " \
		"chromosome as their mate after contig conversion".format(
		n_read_pairs_on_same_chrom_after_conversion, n_reads))

	# Close all chromosome BAM files
	for chrom in chroms:
		chrom_bamfile_dict[chrom].close()

	# Close input BAM
	input_bam.close()


if __name__ == '__main__':
	arguments = parse_args(sys.argv[1:])
	sort_BAM(arguments.input_bam, arguments.output_dir)

