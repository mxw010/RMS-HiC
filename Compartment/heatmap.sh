#!/bin/bash

#SBATCH --cpus-per-task=10
ref=$1
resolution=$2
cell=$3


#cwd:
#/home/gdstantonlab/mxw010/Data/SPICE-C/RH30/Compartments

ml GCCcore/9.3.0
#to generate bigwigs for compartments
ml Python/3.8.2
ml Miniconda3/4.7.10
eval "$(conda shell.bash hook)"
conda activate hic
ml ucsc/default
ml BEDTools/2.29.2

if [[ $resolution -lt 1000000 ]]; then
    res=`echo "$resolution/1000" | bc`
    res=`echo ${res}kb`
else
    res=`echo "$resolution/1000000" | bc`
    res=`echo ${res}mb`
fi

ref_dir="/home/gdstantonlab/mxw010/Data/SPICE-C/tmp/ENCFF643CGH"

echo "assigning compartment with $ref at resolutioon $res"

if [[ $ref == "gc" ]]; then
    ref_file=`echo "${ref_dir}/gc_${resolution}.txt"`
else
    ref_file=`echo "${cell}_K27ac_${res}.bedgraph"`
fi
dir="${ref}_${res}"
rm -rf $dir
mkdir -p $dir
cd $dir
#K27ac


k27_bigwig=`find /home/gdstantonlab/mxw010/Data/Ben/pc-chip-seq/GSL-BS-2005/H3K27AC/${cell}/align_primary/chip/ -name "*.nodup.fc.signal.bigwig" | grep -v glob | head -1`
multiBigwigSummary bins -b $k27_bigwig  \
--binSize $resolution --numberOfProcessors 10 --outRawCounts ${cell}_K27ac_${res}.bedgraph -out test.npz
tail -n+2 ${cell}_K27ac_${res}.bedgraph  | grep -v "random" | grep -v "chrUn" | sort -k1,1 -k2,2n > junk
mv junk ${cell}_K27ac_${res}.bedgraph


#k9me3
#new k9me3 Aug/2022
multiBigwigSummary bins -b /gpfs0/home2/gdstantonlab/mxw010/Data/Ben/H3K9me3-combined/chip/e2adea0b-6160-4fe3-97de-c09419e261a4/call-macs2_signal_track_pooled/execution/glob-8876d8ced974dc46a0c7a4fac20a3a95/rep.pooled_x_ctl.pooled.fc.signal.bigwig  \
--binSize $resolution --numberOfProcessors 10 --outRawCounts ${cell}_K9me3_${res}.bedgraph -out test.npz
tail -n+2 ${cell}_K9me3_${res}.bedgraph   | grep -v "random" | grep -v "chrUn" | sort -k1,1 -k2,2n > junk
mv junk ${cell}_K9me3_${res}.bedgraph

#MYOD
myod_bigwig=`find /home/gdstantonlab/mxw010/Data/Ben/MYOD/GSL-BS-2095/MYOD1/${cell}/align_primary/chip/ -name "*.nodup.fc.signal.bigwig" | grep -v glob | head -1`
multiBigwigSummary bins -b $myod_bigwig  \
--binSize $resolution --numberOfProcessors 10 --outRawCounts ${cell}_MYOD_${res}.bedgraph -out test.npz
tail -n+2 ${cell}_MYOD_${res}.bedgraph   | grep -v "random" | grep -v "chrUn" | sort -k1,1 -k2,2n > junk
mv junk ${cell}_MYOD_${res}.bedgraph

#PAX3-FOXO1
multiBigwigSummary bins -b /home/gdstantonlab/mxw010/Data/Ben/downsampled/data/bigwig/${cell}_pooled.bigwig \
--binSize $resolution --numberOfProcessors 10 --outRawCounts ${cell}_P3F_${res}.bedgraph -out test.npz
tail -n+2 ${cell}_P3F_${res}.bedgraph   | grep -v "random" | grep -v "chrUn" | sort -k1,1 -k2,2n > junk
mv junk ${cell}_P3F_${res}.bedgraph


rm -rf *compA_${res}.bed
rm -rf *compB_${res}.bed

#call compartment
cooltools call-compartments --reference-track $ref_file --bigwig -o ${cell}_${res} ../../cools/${cell}.mcool::/resolutions/${resolution}
tail -n+2 ${cell}_${res}.cis.vecs.tsv | awk 'BEGIN {FS="\t"; OFS="\t"} $6>0 {print $1, $2, $3, $6}' > ${cell}_compA_${res}.bed
tail -n+2 ${cell}_${res}.cis.vecs.tsv | awk 'BEGIN {FS="\t"; OFS="\t"} $6 != "" && $6<0 {print $1, $2, $3, $6}' > ${cell}_compB_${res}.bed


rm -rf  H3K27ac_compartment_by_chr_${res}.bed  H3K27ac_compartment_underpeak_${res}.bed H3K9me3_compartment_by_chr_${res}.bed H3K9me3_compartment_underpeak_${res}.bed
rm -rf MYOD_compartment_by_chr_${res}.bed MYOD_compartment_underpeak_${res}.bed P3F_compartment_by_chr_${res}.bed P3F_compartment_underpeak_${res}.bed
rm -rf gcContent_compartment_by_chr_${res}.bed
#getting average signal per compartment
#k27ac
for i in {1..22}; do
    chr=`echo "chr$i"`
    compA=`bedtools intersect -u -a ${cell}_K27ac_${res}.bedgraph -b ${cell}_compA_${res}.bed  | grep -w $chr | awk '{n++; sum+=$4} END {print NR ? sum/NR : 0}'`
    compB=`bedtools intersect -u -a ${cell}_K27ac_${res}.bedgraph -b ${cell}_compB_${res}.bed  | grep -w $chr | awk '{n++; sum+=$4} END {print NR ? sum/NR : 0}'`
    lengthA=`grep -w $chr ${cell}_compA_${res}.bed | awk '$4 > 0 {print}' | wc -l`
    lengthB=`grep -w $chr ${cell}_compB_${res}.bed | awk '$4 < 0 {print}' | wc -l`
    echo "$chr $compA $compB $lengthA $lengthB" >> H3K27ac_compartment_by_chr_${res}.bed
done


#k27, region under peak
for i in {1..22}; do
    chr=`echo "chr$i"`
    compA=`bedtools intersect -a ../H3K27ac.bed -b ${cell}_compA_${res}.bed | grep -w $chr | awk '{f=$3-$2; k+=$7/f} END {print k}'`
    compB=`bedtools intersect -a ../H3K27ac.bed -b ${cell}_compB_${res}.bed | grep -w $chr | awk '{f=$3-$2; k+=$7/f} END {print k}'`
    lengthA=`grep -w $chr ${cell}_compA_${res}.bed | awk '$4 > 0 {print}' | wc -l`
    lengthB=`grep -w $chr ${cell}_compB_${res}.bed | awk '$4 < 0 {print}' | wc -l`
    echo "$chr $compA $compB $lengthA $lengthB" >> H3K27ac_compartment_underpeak_${res}.bed
done


#k9me3
for i in {1..22}; do
    chr=`echo "chr$i"`
    compA=`bedtools intersect -u -a ${cell}_K9me3_${res}.bedgraph  -b ${cell}_compA_${res}.bed  | grep -w $chr | awk '{n++; sum+=$4} END {print NR ? sum/NR : 0}'`
    compB=`bedtools intersect -u -a ${cell}_K9me3_${res}.bedgraph  -b ${cell}_compB_${res}.bed  | grep -w $chr | awk '{n++; sum+=$4} END {print NR ? sum/NR : 0}'`
    lengthA=`grep -w $chr ${cell}_compA_${res}.bed | awk '$4 > 0 {print}' | wc -l`
    lengthB=`grep -w $chr ${cell}_compB_${res}.bed | awk '$4 < 0 {print}' | wc -l`
    echo "$chr $compA $compB $lengthA $lengthB" >> H3K9me3_compartment_by_chr_${res}.bed
done

#k9me3, region under peak
for i in {1..22}; do
    chr=`echo "chr$i"`
    compA=`bedtools intersect -a ../H3K9me3.bed -b ${cell}_compA_${res}.bed | grep -w $chr | awk '{f=$3-$2; k+=$7/f} END {print k}'`
    compB=`bedtools intersect -a ../H3K9me3.bed -b ${cell}_compB_${res}.bed | grep -w $chr | awk '{f=$3-$2; k+=$7/f} END {print k}'`
    lengthA=`grep -w $chr ${cell}_compA_${res}.bed | awk '$4 > 0 {print}' | wc -l`
    lengthB=`grep -w $chr ${cell}_compB_${res}.bed | awk '$4 < 0 {print}' | wc -l`
    echo "$chr $compA $compB $lengthA $lengthB" >> H3K9me3_compartment_underpeak_${res}.bed
done

#MYOD
for i in {1..22}; do
    chr=`echo "chr$i"`
    compA=`bedtools intersect -u -a ${cell}_P3F_${res}.bedgraph  -b ${cell}_compA_${res}.bed  | grep -w $chr | awk '{n++; sum+=$4} END {print NR ? sum/NR : 0}'`
    compB=`bedtools intersect -u -a ${cell}_P3F_${res}.bedgraph  -b ${cell}_compB_${res}.bed  | grep -w $chr | awk '{n++; sum+=$4} END {print NR ? sum/NR : 0}'`
    lengthA=`grep -w $chr ${cell}_compA_${res}.bed | awk '$4 > 0 {print}' | wc -l`
    lengthB=`grep -w $chr ${cell}_compB_${res}.bed | awk '$4 < 0 {print}' | wc -l`
    echo "$chr $compA $compB $lengthA $lengthB" >> MYOD_compartment_by_chr_${res}.bed
done

#MYOD1, region under peak
for i in {1..22}; do
    chr=`echo "chr$i"`
    compA=`bedtools intersect -a ../MYOD1.bed -b ${cell}_compA_${res}.bed | grep -w $chr | awk '{f=$3-$2; k+=$7/f} END {print k}'`
    compB=`bedtools intersect -a ../MYOD1.bed -b ${cell}_compB_${res}.bed | grep -w $chr | awk '{f=$3-$2; k+=$7/f} END {print k}'`
    lengthA=`grep -w $chr ${cell}_compA_${res}.bed | awk '$4 > 0 {print}' | wc -l`
    lengthB=`grep -w $chr ${cell}_compB_${res}.bed | awk '$4 < 0 {print}' | wc -l`
    echo "$chr $compA $compB $lengthA $lengthB" >> MYOD_compartment_underpeak_${res}.bed
done

#P3F
for i in {1..22}; do
    chr=`echo "chr$i"`
    compA=`bedtools intersect -u -a ${cell}_MYOD_${res}.bedgraph  -b ${cell}_compA_${res}.bed  | grep -w $chr | awk '{n++; sum+=$4} END {print NR ? sum/NR : 0}'`
    compB=`bedtools intersect -u -a ${cell}_MYOD_${res}.bedgraph  -b ${cell}_compB_${res}.bed  | grep -w $chr | awk '{n++; sum+=$4} END {print NR ? sum/NR : 0}'`
    lengthA=`grep -w $chr ${cell}_compA_${res}.bed | awk '$4 > 0 {print}' | wc -l`
    lengthB=`grep -w $chr ${cell}_compB_${res}.bed | awk '$4 < 0 {print}' | wc -l`
    echo "$chr $compA $compB $lengthA $lengthB" >> P3F_compartment_by_chr_${res}.bed
done
#P3F, region under peak
for i in {1..22}; do
    chr=`echo "chr$i"`
    compA=`bedtools intersect -a ../P3F.bed -b ${cell}_compA_${res}.bed | grep -w $chr | awk '{f=$3-$2; k+=$7/f} END {print k}'`
    compB=`bedtools intersect -a ../P3F.bed -b ${cell}_compB_${res}.bed | grep -w $chr | awk '{f=$3-$2; k+=$7/f} END {print k}'`
    lengthA=`grep -w $chr ${cell}_compA_${res}.bed | awk '$4 > 0 {print}' | wc -l`
    lengthB=`grep -w $chr ${cell}_compB_${res}.bed | awk '$4 < 0 {print}' | wc -l`
    echo "$chr $compA $compB $lengthA $lengthB" >> P3F_compartment_underpeak_${res}.bed
done


#gc content
for i in {1..22}; do
    chr=`echo "chr$i"`
    compA=`grep -w $chr ${cell}_compA_${res}.bed | bedtools intersect -u -a ${ref_dir}/gc_${resolution}.txt -b - | awk '{n++; sum+=$4} END {print NR ? sum/NR : 0}' `
    compB=`grep -w $chr ${cell}_compB_${res}.bed | bedtools intersect -u -a ${ref_dir}/gc_${resolution}.txt -b - | awk '{n++; sum+=$4} END {print NR ? sum/NR : 0}' `
    lengthA=`grep -w $chr ${cell}_compA_${res}.bed | awk '$4 > 0 {print}' | wc -l`
    lengthB=`grep -w $chr ${cell}_compB_${res}.bed | awk '$4 < 0 {print}' | wc -l`
    echo "$chr $compA $compB $lengthA $lengthB" >> gcContent_compartment_by_chr_${res}.bed
done

#plotting
computeMatrix scale-regions  \
-m $resolution -bs 1000 \
-S /home/gdstantonlab/mxw010/Data/Ben/P3F_hetero/H3K9me3_result/file_transfer/modified/modified_rep1_FC.bigwig \
/home/gdstantonlab/mxw010/Data/Ben/pc-chip-seq/GSL-BS-2005/H3K27AC/${cell}/align_primary/chip/8b868a24-1eac-425e-9524-88c372566b95/call-macs2_signal_track/shard-0/execution/SEG0134_S2_L003_R1_001_val_1.srt.nodup_x_SEG0133_S1_L003_R1_001_val_1.srt.nodup.fc.signal.bigwig \
/home/gdstantonlab/mxw010/Data/Ben/MYOD/GSL-BS-2095/MYOD1/${cell}/align_primary/chip/e353e904-834c-4c17-80e1-955a4b00deb9/call-macs2_signal_track/shard-0/execution/SEG0181_S5_L001_R1_001_val_1.srt.nodup_x_SEG0177_S1_L001_R1_001_val_1.srt.nodup.fc.signal.bigwig \
/home/gdstantonlab/mxw010/Data/Ben/downsampled/data/bigwig/${cell}_pooled.bigwig \
-R ${cell}_compA_${res}.bed ${cell}_compB_${res}.bed \
-o tf.matrix.gz \
-p 10 \
-bl /home/gdstantonlab/mxw010/GenRef/hg38/hg38.blacklist.bed

plotHeatmap -m tf.matrix.gz \
--samplesLabel "H3K9me3" "H3K27ac" "MYOD" "PAX3-FOXO1" \
--regionsLabel "Comp A" "Comp B" \
--plotTitle "TFs in Compartment Analysis (resolution ${res})" \
-o tfs.pdf \
--startLabel 0 --endLabel $res

plotProfile -m tf.matrix.gz \
--samplesLabel "H3K9me3" "H3K27ac" "MYOD" "PAX3-FOXO1" \
--regionsLabel "Comp A" "Comp B" \
--plotTitle "TFs in Compartment Analysis (resolution ${res})" \
--startLabel 0 --endLabel $res \
-o tf_${res}_profile.pdf

computeMatrix reference-point  \
-a 3000 -b 3000  -bs 50 \
--referencePoint center \
-S /home/gdstantonlab/mxw010/Data/Ben/downsampled/data/bigwig/${cell}_pooled.bigwig \
/home/gdstantonlab/mxw010/Data/Ben/MYOD/GSL-BS-2095/MYOD1/${cell}/align_primary/chip/e353e904-834c-4c17-80e1-955a4b00deb9/call-macs2_signal_track/shard-0/execution/SEG0181_S5_L001_R1_001_val_1.srt.nodup_x_SEG0177_S1_L001_R1_001_val_1.srt.nodup.fc.signal.bigwig \
/home/gdstantonlab/mxw010/Data/Ben/pc-chip-seq/GSL-BS-2005/H3K27AC/${cell}/align_primary/chip/8b868a24-1eac-425e-9524-88c372566b95/call-macs2_signal_track/shard-0/execution/SEG0134_S2_L003_R1_001_val_1.srt.nodup_x_SEG0133_S1_L003_R1_001_val_1.srt.nodup.fc.signal.bigwig \
/home/gdstantonlab/mxw010/Data/Ben/P3F_hetero/H3K9me3_result/file_transfer/modified/modified_rep1_FC.bigwig \
/home/gdstantonlab/mxw010/Data/Ben/P3F_hetero/H3K9me3_result/file_transfer/standard/standard_rep1_FC.bigwig  \
-R P3F_compA.bed P3F_compB.bed \
-o compartment_intensity.matrix.gz \
-p 10 \
-bl /home/gdstantonlab/mxw010/GenRef/hg38/hg38.blacklist.bed

plotHeatmap -m compartment_intensity.matrix.gz \
--samplesLabel "P3F" "MYOD" "K27ac" "K9me3 2X" "K9me3 standard" \
--regionsLabel "P3F Comp A" "P3F Comp B" \
--plotTitle "called with $type (resolution ${res})" \
-o ../${type}_${res}_intensity.pdf \
--yMin 0 --yMax 2 \
--refPointLabel "" \
--startLabel 0 --endLabel $res


cd ..
