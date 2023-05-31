#!/bin/bash

SAMPLENAME=$1
DATADIR=$2
OUTDIR=$3
index_dir=$4

scriptdir="/home/gdstantonlab/mxw010/Data/SPICE-C/scripts"

#create output dir if it doesn't exist
#delete all output files
#create folders
if [ ! -e $OUTDIR ]; then
    mkdir $OUTDIR
fi

cd $OUTDIR
#rm -rf *

mkdir -p logs
mkdir -p pairs
mkdir -p cools
mkdir -p QC
#get absolute path
OUTDIR=`pwd`

#calculate the number of job arrays for aligning reads
n=`ls $DATADIR/*R1*  | wc -l`
echo "aligning $n dataset for $SAMPLENAME"
#alignment
#set up dependency
# jobID=`sbatch --export=$sample=$SAMPLENAME,n=$n --job-name=align-${sample} --array=1-$n \
#  --output=$OUTDIR/logs/align-%a.out ${scriptdir}/half.sh $DATADIR $OUTDIR/pairs $index_dir`
# jobID=`echo $jobID | awk '{print $NF}'`
# echo "Start aligning, jobID=$jobID"

#dedup, normalize and make cool files
# jobID1=`sbatch --dependency=afterok:$jobID --export=sample=$SAMPLENAME --job-name=normalize-$sample \
#     --output=$OUTDIR/logs/normalize.out ${scriptdir}/normalize.sh $SAMPLENAME $OUTDIR $index_dir`
jobID1=`sbatch  --export=sample=$SAMPLENAME --job-name=normalize-$sample \
    --output=$OUTDIR/logs/normalize.out ${scriptdir}/normalize.sh $SAMPLENAME $OUTDIR $index_dir`
jobID1=`echo $jobID | awk '{print $NF}'`
echo "start normalization, jobId=$jobID1"
#sbatch --dependency=afterok:$jobID1 --export=SAMPLENAME=$SAMPLENAME --array=1-4 --job-name=${SAMPLENAME}-compartments ${scriptdir}/call-compartment.sh $SAMPLENAME
