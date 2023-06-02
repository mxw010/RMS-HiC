#!/bin/bash

ml NeoLoopFinder/0.2.4.post2
#joblib-0.17.0
#sklearn downgrade as well
#calculate-cnv -H ../RH4.mcool::/resolutions/50000 --cachefolder cache --output RH4_cnv_norm.txt
#segment-cnv --cnv-file RH4_cnv_norm.txt --binsize 50000 --output RH4_segment.txt --nproc 10
#plot-cnv --cnv-profile RH4_cnv_norm.txt --cnv-segment RH4_segment.txt --output-figure-name cnv.pdf
#correct-cnv --cnv-file RH4_cnv_norm.txt -H ../RH4.mcool::/resolutions/50000 --force


awk '{print $2, $6, $5$9, ($5=="+")?$3:$4, ($9=="+")?$8:$7, "note" }' \
    ~/Data/SPICE-C/hic_breakfinder/RH4/RH4.breaks.txt  > RH4_breaks.txt

python prepare-SV-breakpoints.py ~/Data/SPICE-C/hic_breakfinder/RH4/RH4.breaks.txt RH4_breaks-10kb.txt 10kb
assemble-complexSVs -O SV-10kb -H ../RH4.mcool::/resolutions/50000 -B RH4_breaks-10kb.txt
neoloop-caller -O RH4_neoloops-10kb.txt -H ../RH4.mcool::/resolutions/50000 --assembly SV-10kb.assemblies.txt --no-clustering --prob 0.95 -R 500000
neotad-caller -O RH4_neoTADs-10kb.txt -H ../RH4.mcool::/resolutions/50000 --assembly SV-10kb.assemblies.txt -R 50000 --window-size 50000
#searchSVbyGene --loop-file RH4_neoloops-10kb.txt  --ensembl-release  97 -G PAX3 FOXO1


python prepare-SV-breakpoints.py ~/Data/SPICE-C/hic_breakfinder/RH4_1kb/RH4.breaks.txt RH4_breaks_1kb.txt 1kb
assemble-complexSVs -O SV-1kb -H ../RH4.mcool::/resolutions/50000 -B RH4_breaks_1kb.txt
neoloop-caller -O RH4_neoloops-1kb.txt -H ../RH4.mcool::/resolutions/50000 --assembly SV-1kb.assemblies.txt --no-clustering --prob 0.95
neotad-caller -O RH4_neoTADs-1kb.txt -H ../RH4.mcool::/resolutions/50000 --assembly SV-1kb.assemblies.txt -R 50000 --window-size 50000
# searchSVbyGene --loop-file RH4_neoloops-1kb.txt  --ensembl-release  97 -G PAX3
# searchSVbyGene --loop-file RH4_neoloops-1kb.txt  --ensembl-release  97 -G FOXO1
