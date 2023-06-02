#!/bin/bash

#SBATCH --cpus-per-task=6
#SBATCH --partition=himem
### $1 is read1
### $2 is read2
### $3 is reference genome
### $4 is number of threads
### $5 is the output name

# if [ $# != 3 ]; then
#     echo "Usage: ./hic_align_pipeline.sh <read1> <read2> <reference genome> <number of threads> <output file name - will be *.bam>"
# else
ml HiNT/2.2.8; ml SAMtools/1.10

cell=$1

datadir="/home/gdstantonlab/mxw010/Data/SPICE-C/data"
scriptdir="/home/gdstantonlab/mxw010/Data/SPICE-C/hic_breakfinder/script"
ref_genome="/gpfs0/home2/gdstantonlab/mxw010/Data/SPICE-C/tmp/ENCFF643CGH/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"


jobID1=`sbatch --cpus-per-task=6 --partition=himem --export=cell=$cell --job-name=align_${cell}_R1 --output=$cell.r1.out --wrap="ml HiNT/2.2.8; ml SAMtools/1.10; ${scriptdir}/bwa_mem_hic_aligner.pl ${datadir}/${cell}_R1.fastq.gz $ref_genome 6 |\
samtools view -bS -o ${cell}/${cell}\_read1.bam -"`
jobID1=`echo $jobID1 | awk '{print $NF}'`

jobID2=`sbatch --cpus-per-task=6 --partition=himem --export=cell=$cell --job-name=align_${cell}_R1 --output=$cell.r2.out --wrap="ml HiNT/2.2.8; ml SAMtools/1.10; ${scriptdir}/bwa_mem_hic_aligner.pl ${datadir}/${cell}_R2.fastq.gz $ref_genome 6 |\
samtools view -bS -o ${cell}/${cell}\_read2.bam -"`
jobID2=`echo  $jobID2 | awk '{print $NF}'`

jobID3=`sbatch --dependency=afterok:$jobID1:$jobID2 --export=cell=$cell --job-name=Bam_${cell}  --output=$cell.makebam.out --wrap="ml HiNT/2.2.8; ml SAMtools/1.10; ${scriptdir}/two_read_bam_combiner.pl --qual 30 --single --file1 ${cell}/${cell}\_read1.bam --file2 ${cell}/${cell}\_read2.bam |\
samtools view -u -o - - |\
samtools sort -T 6 -o ${cell}/${cell}.bam -"`
jobID3=`echo  $jobID3 | awk '{print $NF}'`

exp_folder="/home/gdstantonlab/mxw010/Data/SPICE-C/expect_files"
sbatch --dependency=afterok:$jobID3 --export=cell=$cell --output=$cell.findbreaks.out --job-name=findBreaks_${cell} --wrap="ml HiNT/2.2.8; ml SAMtools/1.10; hic_breakfinder --bam-file ${cell}/${cell}.bam  --exp-file-inter ${exp_folder}/inter_expect_1Mb.hg38.txt \
--exp-file-intra ${exp_folder}/intra_expect_100kb.hg38.txt  \
--name  ${cell}/${cell}"


# fi
