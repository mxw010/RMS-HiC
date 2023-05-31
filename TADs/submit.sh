#!/bin/bash

cell=$1
resolution=$2

if [[ $resolution -le 1000000 ]]; then
    res=`echo "$resolution/1000" | bc`
    res=`echo "${res}kb"`
else
    res=`echo "$resolution/1000000" | bc`
    res=`echo "${res}mb"`
fi

# rm -rf data* result*
mkdir -p data_$res
mkdir -p result_$res

n=`cooler dump -t chroms /home/gdstantonlab/mxw010/Data/SPICE-C/${cell}/cools/${cell}.mcool::/resolutions/${resolution} | awk '{print $1}' | wc -l`
jobID=`sbatch --job-name=$cell --array=1-${n} RunTopDom.sh $cell $resolution`
jobID=`echo $jobID | awk '{print $NF}'`

# sbatch --dependency=afterok:$jobID --wrap="ls result_${res}/*.bed  | xargs -i cat {} >> result;
# sort -k1,1 -k2,2n result | grep boundary  > ${cell}_${res}_boundary.bed;
# sort -k1,1 -k2,2n result | grep domain > ${cell}_${res}_domain.bed;"

