#!/bin/bash

#SBATCH --cpus-per-task=10
#SBATCH --partition=himem

datadir=$1
outdir=$2
index_dir=$3


#ml Python/3.8.2
# ml GCCcore/10.3.0
# ml load GCC/10.3.0
# module load BWA/0.7.17
# ml GCC/9.3.0
# ml HiNT/2.2.8
# ml SAMtools/1.10
# source ~/python-virtual-environments/cooltools/bin/activate

#temp solution for using the same reference genome as ENCODE

indx_prefix=`realpath $index_dir/*.fna`
chrom_sizes=`realpath $index_dir/main_chrom.size`

#index=/home/gdlessnicklab/mxw010/GenRef/${genome}/Sequence/BWAIndex/version0.6.0/genome.fa
#chorm_sizes="/home/gdlessnicklab/mxw010/GenRef/${genome}/${genome}.chrom.sizes.main"

n_count=`echo $SLURM_ARRAY_TASK_ID`
prefix=`ls ${datadir}/*R1* | xargs -i basename {} | sed 's/.R1.fastq.gz//g' | head -n${n_count} | tail -1`
SORTED_PAIRS_PATH=${outdir}/${prefix}.sorted.pairs.gz

read1=`ls ${datadir}/${prefix}.R1*`
read2=`ls ${datadir}/${prefix}.R2*`

echo "file to align is $read1, $read2"
#align and conver to pairs format
#bwa mem -SP5M  -t 10 $index $read1 $read2 | samtools view -@ 9 -Shb - |
#bwa mem -SP5M  -t 10 $indx_prefix $read1 $read2 | samtools view -@ 10 -Shb -o ${prefix}_temp.bam

ml GCC/9.3.0 SAMtools/1.10
ml load Miniconda3/4.9.2
eval "$(conda shell.bash hook)"
conda activate hic
#from here on, drop chrom other than 1 - XY
pairtools parse --add-columns mapq --drop-sam --chroms-path ${chrom_sizes}  ${prefix}_temp.bam |
pairtools sort --nproc 4 -o ${SORTED_PAIRS_PATH}
