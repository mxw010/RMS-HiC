#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --partition=himem

cell=$1

ml GCCcore/9.3.0
#to generate bigwigs for compartments
ml Python/3.8.2
ml ucsc/default
ml Miniconda3/4.7.10
eval "$(conda shell.bash hook)"
conda activate hic

ref_dir="/home/gdstantonlab/mxw010/Data/SPICE-C/tmp/ENCFF643CGH"

# mkdir -p Compartments
cd Compartments
#resolutions=(100000 500000 1000000 5000000)
resolution=100000
#n_count=`echo $SLURM_ARRAY_TASK_ID`
#res=${resolutions[$n_count-1]}

if [[ $resolution -lt 1000000 ]]; then
    res=`echo "$resolution/1000" | bc`
    res=`echo ${res}kb`
else 
    res=`echo "$resolution/1000000" | bc`
    res=`echo ${res}mb`
fi 
#normalization
hicConvertFormat --matrices ../cools/${cell}.mcool::resolutions/$resolution  --outFileName ${cell}_${res}.h5 --inputFormat cool \
    --outputFormat h5 --load_raw_values --resolutions $resolution
hicCorrectMatrix diagnostic_plot --matrix ${cell}_${res}.h5 -o ${cell}_correction_plot.png 2> ${cell}_${res}_output.txt
min=`grep "mad threshold" ${cell}_${res}_output.txt | awk '{print $NF}'`
echo "min is $min"
hicCorrectMatrix correct --correctionMethod ICE --matrix ${cell}_${res}.h5 -o ${cell}_${res}_corrected.h5  --filterThreshold $min 3
rm ${cell}_${res}_output.txt ${cell}_correction_plot.png
#saddle 
hicConvertFormat --matrices ${cell}_${res}_corrected.h5  --outFileName ${cell}_${res}_hicexplorer_corrected.cool \
                  --inputFormat h5 --outputFormat cool
cooltools call-compartments --reference-track ${ref_dir}/gc_${resolution}.txt --bigwig -o ${cell}_${res} ${cell}_${res}_hicexplorer_corrected.cool
cooltools compute-expected --nproc 4 -o ${cell}_${res}_expected_tracks.tsv ${cell}_${res}_hicexplorer_corrected.cool 
cooltools compute-saddle -o ${cell}_saddle_${res} --qrange 0.02 0.98 --fig pdf ${cell}_${res}_hicexplorer_corrected.cool \
    ${cell}_${res}.cis.vecs.tsv ${cell}_${res}_expected_tracks.tsv 

# mkdir -p cooltools
# mv *.tsv cooltools/.
# mv *.txt cooltools/.

# mv ${cell}_${res}.cis.bw ${cell}_${res}_cooltools.bigwig
# bigWigToBedGraph ${cell}_${res}_cooltools.bigwig ${cell}_${res}_cooltools.bedgraph
# k27_track=`find /home/gdstantonlab/mxw010/Data/Ben/pc-chip-seq/GSL-BS-2005/H3K27AC/${cell} -name "*.bigwig" | grep -v glob | grep fc | head -1`
# sbatch --partition=himem --wrap="hicPCA -m ${cell}_${res}_corrected.h5 --outputFileName ${cell}_${res}_pca1.bedgraph ${cell}_${res}_pca2.bedgraph --format bedgraph \
# --pearsonMatrix ${cell}_${res}_pearson.h5 \
# --extraTrack $k27_track \
# --histonMarkType active"

#TADs
# cd ..
# mkdir -p TADs
# cd TADs
# resolution=10000
# if [[ $resolution -lt 1000000 ]]; then
#     res=`echo "$resolution/1000" | bc`
#     res=`echo ${res}kb`
# else 
#     res=`echo "$resolution/1000000" | bc`
#     res=`echo ${res}mb`
# fi 

# hicConvertFormat --matrices ../cools/${cell}.mcool::resolutions/$resolution  --outFileName ${cell}_${res}.h5 --inputFormat cool \
#     --outputFormat h5 --load_raw_values --resolutions $resolution
# hicCorrectMatrix diagnostic_plot --matrix ${cell}_${res}.h5 -o ${cell}_correction_plot.png 2> ${cell}_${res}_output.txt
# min=`grep "mad threshold" ${cell}_${res}_output.txt | awk '{print $NF}'`
# echo "min is $min"
# hicCorrectMatrix correct --correctionMethod ICE --matrix ${cell}_${res}.h5 -o ${cell}_${res}_corrected.h5  --filterThreshold $min 4.5
# rm ${cell}_${res}_output.txt ${cell}_correction_plot.png
# #saddle 
# hicConvertFormat --matrices ${cell}_${res}_corrected.h5  --outFileName ${cell}_${res}_hicexplorer_corrected.cool \
#                   --inputFormat h5 --outputFormat cool
# hicFindTADs --matrix ${cell}_${res}_hicexplorer_corrected.cool --correctForMultipleTesting fdr --outPrefix ${cell}_${res} \
# --numberOfProcessors 8

# cd ..	
