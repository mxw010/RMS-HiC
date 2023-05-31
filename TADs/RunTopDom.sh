#!/bin/bash

#SBATCH --partition=himem

cell=$1
resolution=$2

if [[ $resolution -le 1000000 ]]; then
    res=`echo "$resolution/1000" | bc`
    res=`echo "${res}kb"`
else
    res=`echo "$resolution/1000000" | bc`
    res=`echo "${res}mb"`
fi

#chr="chr1"
chr=`cooler dump -t chroms /home/gdstantonlab/mxw010/Data/SPICE-C/${cell}/cools/${cell}.mcool::/resolutions/${resolution} | awk '{print $1}' | head -n${SLURM_ARRAY_TASK_ID} | tail -1`
echo "running cell line $cell at resolution $res for $chr"
# ml Miniconda3/4.7.10
# eval "$(conda shell.bash hook)"
# conda activate hic
python3 convertToTopDom.py --cell $cell --chr $chr --res $resolution
sed -ie 's/nan/0/g' data_${res}/freq_${chr}
paste -d "\t" data_${res}/chrom_${chr} data_${res}/freq_${chr} > data_${res}/${chr}
window=(2,3,5,7,10,12,15,16,18,20)
Rscript --vanilla TopDom.R $chr $res ${window[*]} $cell


