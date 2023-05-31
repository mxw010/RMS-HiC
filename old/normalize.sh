#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --partition=himem

cell=$1
indir=$2
index_dir=$3

#main chrom size
chrom_sizes=`realpath ${index_dir}/main_chrom.size`

echo "indir is: $indir"
echo "chrom size is: $chrom_sizes"
#chrom_sizes="/home/gdlessnicklab/mxw010/GenRef/${genome}/${genome}.chrom.sizes.main"
#chrom_sizes=`echo "/home/gdlessnicklab/mxw010/GenRef/${genome}/${genome}.chrom.sizes"`
res=1000

ml purge
ml GCC/9.3.0 SAMtools/1.10
ml load Miniconda3/4.9.2
eval "$(conda shell.bash hook)"
conda activate hic

pairtools merge --max-nmerge 20 --memory 10000G ${indir}/pairs/*.sorted.pairs.gz |
pairtools sort --memory 10000G --nproc 4 |
pairtools dedup  --nproc-in 4 --output-stats ${indir}/QC/${cell}.dedup.stats |
pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type == "RU")' |
pairtools select --output ${indir}/pairs/${cell}.filtered.pairs.gz '((chrom1==chrom2) and (abs(pos1-pos2))>300) or (chrom1 != chrom2)'
pairtools select 'mapq1 > 30 and mapq2 > 30' ${indir}/pairs/${cell}.filtered.pairs.gz |
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly hg38 ${chrom_sizes}:${res} - ${indir}/cools/${cell}.cool
cooler zoomify -n 10 -r 1000,5000,10000,100000,500000,1000000,5000000 --balance --out ${indir}/cools/${cell}.mcool ${indir}/cools/${cell}.cool
