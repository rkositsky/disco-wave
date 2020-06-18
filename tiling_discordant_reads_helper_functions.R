# Functions to process single-sample discordant reads detection

library(tidyverse)
library(gridExtra)
library(viridis)


#### Input parsing ####
parse_bedtools_coverages <- function(in_dir){
  # Arguments:
  #    in_dir: directory with all chromosome coverages. See run_chrom_coverages.bash for how to produce these.
  # Output:
  #    df: dataframe with 55 columns. BED information about window region + gene region + 
  #        number discordant reads going to all chromosomes + for each chromosome, # disc. reads
  #        going to that chromosome and what percent it makes up of all discordant reads.
  
  # Start dataframe with all regions
  df <- read_delim(paste0(in_dir, "/all.bed"), delim="\t",
                   col_names = c("chrom", "start", "stop", "gene", 
                                 "all", "bases_covered", "region_length", "fraction_covered")) %>%
    select(c("chrom", "start", "stop", "gene", "all"))
  
  chromosomes <- c(unlist(map(1:22, function(x) {sprintf("chr%d", x)})), c("chrX", "chrY", "chrM"))
  
  for(chrom_name in chromosomes) {
    regions <- read_delim(paste0(in_dir, "/", chrom_name, ".bed"), delim="\t",
                          col_names = c("chrom", "start", "stop", "gene", 
                                        chrom_name, "bases_covered", "region_length", "fraction_covered")) %>%
      select(c(chrom_name))
    df[, chrom_name] <- regions %>% pull(chrom_name)
    df[, paste0("pct_", chrom_name)] <- df[, chrom_name] / df[, "all"] * 100
  }
  return(df)
}

#### Plotting functions ####
plot_pct_chrom <- function(sample_df, chrom, title=NA) {
  # Arguments:
  #     df: dataframe filtered for the gene of interest. Expects format from parse_bedtools_coverages.
  #     chrom: chromosome to plot
  #     min_total_reads: threshold to filter on for plotting minimum reads
  # Returns:
  #     ggplot of pct_chrom vs. all_reads, colored by start position (?)
  col_name <- paste0("pct_", chrom)
  sample_df$pct_chrom <- sample_df %>% pull(col_name)
  
  p <- ggplot(data = sample_df, aes(x = all, y = pct_chrom, color=start)) +
    geom_point(size = 0.5) +
    labs(x = "All", y = paste0("% on ", chrom)) +
    scale_y_continuous(limits = c(0, 100)) +
    scale_color_viridis(discrete = F) +
    # remove the legend; it takes up so much space!
    theme(legend.position = "none")
  if(!is.na(title)) {
    p <- p + labs(title=title)
  }
  return(p)
}

plot_pct_chrom_for_gene <- function(sample_df, out_dir, sample_nickname = "Sample",
                                    min_total_reads = 5, 
                                    genes=c("MYC_FISH", "BCL2_FISH", "BCL6_FISH")) {
  # For given sample regions and given gene(s), plot a 5x5 grid of 
  # pct_chr vs. total discordant reads for all chromosomes
  # Arguments:
  #     df: dataframe for a given sample
  #     out_dir: where to write file(s)
  #     sample_nickname: nickname to print at the top of the plot
  #     min_total_reads: threshold to filter on for plotting minimum reads
  #     genes: which genes to plot
  
  chromosomes <- c(unlist(map(1:22, function(x) {sprintf("chr%d", x)})), c("chrX", "chrY", "chrM"))
  
  # Filter for the gene and for minimium total reads
  for(gene_name in genes) {
    gene_df <- sample_df %>% filter(gene==gene_name) %>% filter(all >= min_total_reads)
    
    plot_list <- lapply(chromosomes, function(chrom) {plot_pct_chrom(gene_df, chrom)})
    names(plot_list) <- chromosomes
    
    png(paste0(out_dir, gene_name, "_pct_chroms_vs_all_reads.png"),
        width=18, height = 12, units = "in", res = 300)
    grid.arrange(grobs=plot_list, ncol = 5,
                 left = "% discordant reads mapping to mate chromosome",
                 bottom="Total discordant reads in sample",
                 top = paste(sample_nickname, "-", gene_name, 
                             "region - % discordant reads per chromosome vs. total discordant reads"))
    dev.off()
  }
}


#### Output functions ####
collapse_adjacent_windows <- function(df, merge_distance) {
  # Merge adjacent tiling windows within merge_distance if they have the same 
  # partner chromosome
  # When aggregated, keep information from highest-covered window
  # Arguments:
  #     df: tibble from write_candidate_translocations
  #     merge_distance: maximum distance between two adjacent windows that
  #                     can be merged
  
  # Return immediately if there's no passing rows in dataframe
  if(nrow(df) == 0) {
    return(df)
  }
  
  # Now we can assume there's at least one window. Copy over first row.
  output_df <- df[1,]
  
  # We assume that the input windows are sorted by sample and chromosomal position.
  # If this assumption does not hold, we will not optimally collapse.
  for(i in 2:nrow(df)) {
    j <- nrow(output_df)
    
    # If you can merge, update last row in output df. Otherwise, write new row as-is.
    if (((df[i, "start"] - output_df[j, "stop"]) <= merge_distance) && 
        (df[i, "chrom"] == output_df[j, "chrom"]) && 
        (df[i, "partner_chrom"] == output_df[j, "partner_chrom"]) &&
        (df[i, "sample"] == output_df[j, "sample"])) {
      # Let's merge!
      
      # 1) Extend window end to new stop position
      output_df[j, "stop"] <- df[i, "stop"]
      
      # 2) Update fields (all, max_pct, ... pct_chrM) IF the new window has 
      # more support (more total reads, or equal reads w/ higher %).
      # This means the output has stats from the best-supported window.
      if ((df[i, "all"] > output_df[j, "all"]) || 
          ((df[i, "all"] == output_df[j, "all"]) && 
           (df[i, "pct_reads_to_partner"] > output_df[j, "pct_reads_to_partner"]))) {
        output_df[j, which(colnames(df)=="all"):ncol(df)] <- df[i, which(colnames(df)=="all"):ncol(df)]
      }
    } else {
      # Otherwise add a new row
      output_df[j+1,] <- df[i,]
    }
  }
  
  return(output_df)
}

write_candidate_translocations <- function(
  regions_df, out_dir, sample_name, min_total_reads = 10, min_pct = 75,
  merge_distance = 0) {
  # Outputs candidate translocations to a tab-separated file
  # Arguments:
  #     sample_regions_list: list with dataframe for each sample
  #     out_dir: where to write file
  #     samples: vector with sample names to plot
  #     sample_nicknames: rownames are sample IDs, has a column named "nickname" for plotting
  #     min_total_reads: minimum reads for a window to be kept
  #     min_pct: minimum percent for a window to be kept
  #     merge_distance: maximum distance that windows can be apart when they're merged
  
  out_file <- paste0(out_dir, "detected_translocation_regions.min", 
                     min_total_reads,"reads.min", min_pct, "pct.merge", 
                     merge_distance, "bp.tsv")
  
  # For each sample, filter translocations, merge them, and then append to output dataframe
  # Create new columns for filtering
  df <- regions_df
  df$pct_reads_to_partner <- apply(df[,grep("pct_", colnames(df))], FUN=max, MARGIN=1)
  df$sample <- sample_name

  # Filter by total number of discordant reads and an overabundance of one chromosome
  df <- df %>% filter(all >= min_total_reads) 
  df <- df %>% filter(pct_reads_to_partner >= min_pct)
  
  ## Get the partner chromosome now that you have non-NAN values for things after filtering
  # First get the column names that are pct_chr#.
  pct_colnames <- grep("pct_", colnames(df), value=T)
  # Next, get the first instance of the chromosome with the greatest %. Then
  # get the name of that chromosome  by looking up in pct_colnames. Then
  # convert from pct_chr# to just chr# using substring.
  df$partner_chrom <- substring(pct_colnames[apply(df[,pct_colnames], MARGIN=1, FUN=which.max)], first=5)
  
  # Reorder columns
  df <- df %>% select("sample", 1:5, "pct_reads_to_partner", "partner_chrom", everything())
  
  # Merge adjacent windows if they share the same chromosome
  df <- collapse_adjacent_windows(df, merge_distance)
  
  # Write to output file, overwriting and initializing the column names for 
  # the first sample and appending after that
  write_tsv(df, path=out_file, na = ".")
}


