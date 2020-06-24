#!/usr/local/bin/Rscript


# Created: 2020-01-27
# Updated: 2020-06-17

#### Load argument-parsing libraries (rest are loaded with helper functions) ####
library("optparse")


#### List of command-line arguments ####
options_list <- list(
  make_option(c("-f", "--helper_functions"), dest = "helper_functions_file", 
              type = "character", help = "path to helper functions file"),
  make_option(c("-i", "--input_directory"), dest = "input_directory", 
              type = "character", 
              help = "path to input directory with coverage over tiled windows for 25 chromosomes"),
  make_option(c("-o", "--output_directory"), dest = "output_directory", 
              type = "character", 
              help = "path to output directory where output table and figures are printed"),
  make_option(c("-n", "--sample_name"), dest = "sample_name", 
              type = "character", 
              help = "sample name"),
  make_option(c("-r", "--min_reads"), dest = "min_reads", 
              type = "integer", 
              help = "minimum reads for a translocation tiling window to pass",
              default = 10),
  make_option(c("-d", "--merge_distance"), dest = "merge_distance", 
              type = "integer", 
              help = "maximum distance that windows can be apart when they're merged",
              default = 0)
)

#### Parse arguments ####
opt_parser <- OptionParser(option_list = options_list)
opt <- parse_args(opt_parser)

source(opt$helper_functions_file)

in_dir <- paste0(opt$input_directory, "/")
out_dir <- paste0(opt$output_directory, "/")

#### Read things in ####
regions_df <- parse_bedtools_coverages(in_dir = in_dir)

#### Plot figures of pct reads going to each chromosome vs. # reads ####
plot_pct_chrom_for_gene(sample_df = regions_df,
                        out_dir = out_dir,
                        sample_nickname = opt$sample_name)

# Write output
write_candidate_translocations(regions_df, out_dir, opt$sample_name, 
                               min_total_reads = opt$min_reads,
                               merge_distance = opt$merge_distance)
