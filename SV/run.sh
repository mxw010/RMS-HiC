#!/bin/bash
#SBATCH --partition=himem
#SBATCH --cpus-per-task=10

ml NeoLoopFinder/0.2.4.post2
calculate-cnv -H ../RH4.mcool::/resolutions/10000 --cachefolder cache --output RH4_cnv_norm.txt
segment-cnv --cnv-file RH4_cnv_norm.txt --binsize 10000 --output RH4_segment.txt --nproc 10
plot-cnv --cnv-profile RH4_cnv_norm.txt --cnv-segment RH4_segment.txt --output-figure-name cnv.pdf
