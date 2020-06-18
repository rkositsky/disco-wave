#!/bin/bash

# Rachel Kositsky
# Created: 2019-04-11
# Updated: 2020-06-20

# Goal: given a BAM file, produce another BAM file that selects discordant reads

# Usage if don't have inputs
if [[ $# -ne 4 ]] ; then
    echo "Usage: select_discordant_reads.bash [in.bam] [out.bam] [nr_cpus] [min_map_quality]"
    echo "Given a BAM file, produce another BAM file that selects discordant reads."
    echo "  in.bam, out.bam: input and output file paths"
    echo "  nr_cpus: number of CPUs over which to parallelize"
    echo "  min_map_quality: minimum mapping quality for selected reads"
    exit 0
fi

input=$1
output=$2
nr_cpus=$3
min_mq=$4

# Flag 1294: only keep read pairs where both read and mate are mapped,
# it's the read's primary alignment, and it's not a PCR/optical duplicate
samtools view -h -@ ${nr_cpus} -F 1294 ${input} -q ${min_mq} | awk '$7!="="' | samtools view -b -@ ${nr_cpus} -S - -o ${output}
samtools index ${output}
