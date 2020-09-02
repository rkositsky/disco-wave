# DiscoWave: A bioinformatics tool to identify translocations using discordant reads 

DiscoWave identifies interchromosomal translocations from paired-end aligned
sequencing data. It focuses on just the discordant read pairs that map betwen 
different chromosomes. It accounts for variability in alignment quality at 
different loci by normalizing for the number of discordant reads at each tiling
 window examined.
 
One of the outputs is a scatterplot of stats for each tiling window. For each 
possible partner chromsome, the percentage of discordant reads mapping to that 
partner chromsome vs. the total reads is plotted. The sign of a "hit" is the
shape of a wave - locations in the genome where all discordantly mapping reads 
are going to just one other chromosome, with good coverage.

Currently this tool is geared towards displaying MYC, BCL2, and BCL6 
translocations. Further generalization is in progress.

## Installation

1. Clone this repository and install missing [dependencies](#dependencies)
2. Use Docker to access this repository with dependencies pre-installed

## Usage

### Source usage

```
usage: Main.py [-h] [--out_translocation_table OUT_TRANSLOCATION_TABLE]
               [--out_discordant_bam OUT_DISCORDANT_BAM]
               [--figure_output_dir FIGURE_OUTPUT_DIR]
               [--sample_name SAMPLE_NAME] [--tiling_BED TILING_BED]
               [--min_mapping_quality MIN_MAPPING_QUALITY]
               [--merge_distance MERGE_DISTANCE]
               [--min_read_pairs MIN_READ_PAIRS] [--nr_cpus NR_CPUS]
               in_bam

Find translocations using discordant reads from tiling windows from short-read
sequencing aligned files

positional arguments:
  in_bam                Input BAM file

optional arguments:
  -h, --help            show this help message and exit
  --out_translocation_table OUT_TRANSLOCATION_TABLE
                        File path to write tab-separated table of candidate
                        translocations
  --out_discordant_bam OUT_DISCORDANT_BAM
                        File path to write subsetted BAM with only
                        discordantly mapped reads
  --figure_output_dir FIGURE_OUTPUT_DIR
                        Directory where output figure will be written
  --sample_name SAMPLE_NAME
                        Sample name
  --tiling_BED TILING_BED
                        Output from make_tiling_BEDs.py. BED file with tiled
                        windows where translocations can be called. Default:
                        MYC, BCL2, BCL6 FISH capture regions.
  --min_mapping_quality MIN_MAPPING_QUALITY
                        Minimum read quality for read that get considered.
  --merge_distance MERGE_DISTANCE
                        Maximum distance between reads where they're still
                        considered part of the same translocation. Default: 0
                        bp
  --min_read_pairs MIN_READ_PAIRS
                        Minimum supporting reads for a translocation to be
                        reported. Default: 10 read pairs
  --nr_cpus NR_CPUS     Number of CPUs for parallelizing discordant read
                        selection. Default: 1
```

Example command, given my_sample.bam and my_sample.bam.bai in the current directory:
`python3 Main.py my_sample.bam --sample_name my_sample --nr_cpus 6`

### Docker usage

A Docker image is available at davelabhub:disco-wave.

Example Docker command, given my_sample.bam and my_sample.bam.bai in the current directory:
`sudo docker run --rm --user root -v ${PWD}:/data davelabhub/disco-wave:latest /bin/bash -c "python3 /disco-wave/Main.py /data/my_sample.bam --out_translocation_table /data/candidate_translocations.tsv --out_discordant_bam /data/discordant_reads.diff_chrom.bam --figure_output_dir /data/supporting_figures --sample_name my_sample --nr_cpus 7"`
When using Docker, it's helpful to define the outputs explicitly so that they remain after the Docker command finishes.


## Dependencies

Where available, versions are noted. Other versions may work, but have not been tested.

* python (used 3.6.5)
* python packages:
  - argparse
  - pysam
* samtools (used v1.9)
* awk
* bedtools (used v2.27.1)
* Rscript (used R v3.5.2)
* R packages: 
	- optparse (used v1.6.6)
	- tidyverse
	- gridExtra
	- viridis
  - microbiome 
