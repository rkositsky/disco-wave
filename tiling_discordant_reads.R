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
  make_option(c("-o", "--output_table"), dest = "output_table", 
              type = "character", 
              help = "file path where output table with candidate translocations should be printed"),
  make_option(c("-g", "--figure_output_directory"), dest = "figure_output_directory", 
              type = "character", 
              help = "directory path where output figures should be printed"),
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
fig_out_dir <- paste0(opt$figure_output_directory, "/")

# Make figure output directory if it doesn't exist
if (!dir.exists(fig_out_dir)) {
  dir.create(fig_out_dir)  
}


#### Read things in ####
regions_df <- parse_bedtools_coverages(in_dir = in_dir)

#### Plot figures of pct reads going to each chromosome vs. # reads ####
# Get gene names from the 4th column of the tiling BED files
plot_pct_chrom_for_gene(sample_df = regions_df,
                        out_dir = fig_out_dir,
                        sample_nickname = opt$sample_name,
                        genes = unique(regions_df$gene))

# Write output
write_candidate_translocations(regions_df = regions_df, 
                               out_file = opt$output_table, 
                               sample_name= opt$sample_name,
                               min_total_reads = opt$min_reads,
                               merge_distance = opt$merge_distance)
