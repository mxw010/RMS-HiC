#!/gpfs0/export/apps/opt/NeoLoopFinder/0.2.4.post2/envs/neoloop/bin/python3.7
#
# import io
from pyensembl import EnsemblRelease
from neoloop.visualize.core import *
import cooler
import re
import subprocess
from neoloop.tadtool.core import TADcaller
import joblib
import os
import pandas as pd
import numpy as np
import glob
# clr = cooler.Cooler(
#     '/home/gdstantonlab/mxw010/Data/SPICE-C/RH30/cools/RH30.mcool::/resolutions/10000')
# clr = cooler.Cooler(
#     '/home/gdstantonlab/mxw010/Data/SPICE-C/RH4/cools/RH4.mcool::/resolutions/10000')

type = 'RH4'

res='10kb'

gene_list = [line.rstrip() for line in open(
     '/home/gdstantonlab/mxw010/GenRef/hg38/Annotation/Genes/gene_names')]
clr = cooler.Cooler('/home/gdstantonlab/mxw010/Data/SPICE-C/' +
                    type + '/cools/' + type + '.mcool' + '::/resolutions/10000')

from os.path import exists

df = pd.DataFrame()
#get all potential interesting regions:
loop_file=type + '_neoloops-' + res + '.txt'
tad_file = type + '_neoTADs-' + res + '.txt'
if exists(loop_file):
    temp = pd.read_csv(loop_file, sep="\t", header=None)
    df = pd.concat([df,temp], axis=1)

if exists(tad_file):
    temp = pd.read_csv(tad_file, sep="\t",header=None)
    df = pd.concat([df,temp], axis=0)

#merge the reigon ID:
candidate = []
for itx, row in df.iterrows():
    x = np.array(row[6].split(',')).reshape(int(len(row[6].split(','))/3),3)
    for index, ind in np.ndenumerate(x[:,2]):
        if ind == '1' and (not any(x[index,0].tolist()[0] == y for y in candidate)):
            candidate.append(x[index,0].tolist()[0])

candidate.sort()

def plot_sv(assembly, cell, res, correct='weight'):
    ID=assembly.split("   ")[0]
    outname = cell + "_" +  ID + "_" + res +  '.pdf'
    if cell == 'RH4' or cell == 'RH30':
        vis = Triangle(clr, assembly, n_rows=8, figsize=(7, 4.2), track_partition=[5, 1, 1,0.8,1, 0.5, 1, 0.5], correct=correct)
        #subplot1
        vis.matrix_plot(vmin=0)
        vis.plot_chromosome_bounds(linewidth=2.5)
        #vis.plot_loops('RH4_neoloops-1kb.txt',
        #                face_color='none', marker_size=40, cluster=True)
        #subplot2
        vis.plot_neoTAD(ws=1000000)
        #subplot3, standard DI plot is buggy
        hmm_folder = os.path.join(os.path.split(neoloop.__file__)[0], 'data')
        hmm = joblib.load(os.path.join(hmm_folder, 'HMM-model.pkl'))
        tad_ax = vis.fig.add_subplot(vis.grid[vis.track_count])
        vis.track_count += 1
        tad_ax.set_xlim(vis.hx.min(), vis.hx.max())
        x = range(int(vis.hx.min()), int(vis.hx.max()))
        work = TADcaller(vis.matrix, vis.res, hmm, window_size=1000000)
        work.callDomains()
        y = work.DIs
        d = np.zeros(len(y))
        tad_ax.fill_between(x, d, y, where=y >= d, interpolate=True, color='blue')
        tad_ax.fill_between(x, d, y, where=y <= d, interpolate=True, color='red')
        tad_ax.set_axis_off()
        #subplot4
        vis.plot_signal('P3F', '/home/gdstantonlab/mxw010/Data/Ben/downsampled/data/bigwig/' + cell + '_pval.bigwig',
                        label_size=6, data_range_size=9, max_value=3, color='#6A3D9A')
        #subplot5 and 6
        motif_bed = '/home/gdstantonlab/mxw010/Data/Ben/CTCF/GSL-BS-2483/result/motif/' + '*/*' + cell + '_CTCF_clean.bed'
        motif_path = glob.glob(motif_bed, recursive=True)[0]
        vis.plot_motif(motif_path, subset="+")
        vis.plot_motif(motif_path, subset="-")
        #vis.plot_motif('/home/gdstantonlab/mxw010/Data/Ben/CTCF/GSL-BS-2483/result/motif/' + cell + '/' + cell + '_CTCF_clean.bed', subset='-')
        #vis.plot_signal('CTCF', bigwigPath, label_size=6, data_range_size=9, max_value=3, color='#6A3D9A')
        #subplot7
        vis.plot_genes(release=79, filter_=gene_list, fontsize=3, style='other')
        #subplot8
        vis.plot_chromosome_bar(name_size=11, coord_size=3)
        vis.outfig(outname)
    else:
        vis = Triangle(clr, assembly, n_rows=7, figsize=(
            7, 4.2), track_partition=[5, 1, 1,0.8,0.5, 1, 0.5], correct=correct)
        #subplot1
        vis.matrix_plot(vmin=0)
        vis.plot_chromosome_bounds(linewidth=2.5)
        #vis.plot_loops('RH4_neoloops-1kb.txt',
        #                face_color='none', marker_size=40, cluster=True)
        #subplot2
        vis.plot_neoTAD(ws=1000000)
        #subplot3, standard DI plot is buggy
        hmm_folder = os.path.join(os.path.split(neoloop.__file__)[0], 'data')
        hmm = joblib.load(os.path.join(hmm_folder, 'HMM-model.pkl'))
        tad_ax = vis.fig.add_subplot(vis.grid[vis.track_count])
        vis.track_count += 1
        tad_ax.set_xlim(vis.hx.min(), vis.hx.max())
        x = range(int(vis.hx.min()), int(vis.hx.max()))
        work = TADcaller(vis.matrix, vis.res, hmm, window_size=1000000)
        work.callDomains()
        y = work.DIs
        d = np.zeros(len(y))
        tad_ax.fill_between(x, d, y, where=y >= d, interpolate=True, color='blue')
        tad_ax.fill_between(x, d, y, where=y <= d, interpolate=True, color='red')
        tad_ax.set_axis_off()
        #subplot5
        #bigwigName = '/home/gdstantonlab/mxw010/Data/Ben/CTCF/GSL-BS-2483/result/file_transfer/' + '*/*' + type + '_pooled_pval.bigwig'
        motif_bed = '/home/gdstantonlab/mxw010/Data/Ben/CTCF/GSL-BS-2483/result/motif/' + '*/*' + cell + '_CTCF_clean.bed'
        motif_path = glob.glob(motif_bed, recursive=True)[0]
        vis.plot_motif(motif_path, subset="+")
        vis.plot_motif(motif_path, subset="-")
        #vis.plot_motif('/home/gdstantonlab/mxw010/Data/Ben/CTCF/GSL-BS-2483/result/motif/RH4/CTCF_motif1_clean.bed', subset='-')
        #subplot6
        vis.plot_genes(release=79, filter_=gene_list, fontsize=3, style='other')
        #subplot7
        vis.plot_chromosome_bar(name_size=11, coord_size=3)
        vis.outfig(outname)
    plt.close(fig=vis.fig)


#P3F locus is C36 at 1kb
# add ctcf
with open('SV-'+ res + '.assemblies.txt') as file1:
    for line in file1:
        #indicator = line.decode().strip().replace('\t', '  ')
        assembly = line.strip().replace('\t', '   ')
        ID=assembly.split("   ")[0]
        if any(ID == item for item in candidate):
            print('plotting ' + ID + "...")
            plot_sv(assembly)
        else:
            print('skipping ' + ID + "...")


#RH30, C44, 1kb
assembly = 'A0      translocation,2,222203510,-,9,86635521,-        translocation,9,86702560,+,13,40620000,+        2,222340000     13,39431592'
type = 'RH4'
clr = cooler.Cooler('/home/gdstantonlab/mxw010/Data/SPICE-C/' +
                    type + '/cools/' + type + '.mcool' + '::/resolutions/10000')
ID=assembly.split("   ")[0]
plot_sv(assembly, 'RH4', res)

#RD, C44, 1kb
type = 'RD'
clr = cooler.Cooler('/home/gdstantonlab/mxw010/Data/SPICE-C/' +
                    type + '/cools/' + type + '.mcool' + '::/resolutions/10000')
assembly = 'C44   translocation,13,40614000,+,2,222203000,-   13,39000000   2,223460000'
ID=assembly.split("   ")[0]
plot_sv(assembly, 'RD', res)

type = 'CTR'
clr = cooler.Cooler('/home/gdstantonlab/mxw010/Data/SPICE-C/' +
                    type + '/cools/' + type + '.mcool' + '::/resolutions/10000')
assembly = 'C44   translocation,13,40614000,+,2,222203000,-   13,39000000   2,223460000'
ID=assembly.split("   ")[0]
plot_sv(assembly, 'CTR', res)



#PAX3-FOXO1, C36, 1kb
res='1kb'
assembly='C36   translocation,13,40607000,+,2,222203000,-   13,39440000   2,222440000'
ID=assembly.split("   ")[0]
outname = 'RH4.' + ID + "_" +  res + '.pdf'
vis = Triangle(clr, assembly, n_rows=8, figsize=(
    7, 4.2), track_partition=[5, 1, 1,0.8,1, 1, 1, 0.5], correct=None)
#subplot1
vis.matrix_plot(vmin=0)
vis.plot_chromosome_bounds(linewidth=2.5)
vis.plot_loops('RH4_neoloops-1kb.txt',
                face_color='none', marker_size=40, cluster=True)
#subplot2
vis.plot_neoTAD(ws=1000000)
#subplot3, standard DI plot is buggy
hmm_folder = os.path.join(os.path.split(neoloop.__file__)[0], 'data')
hmm = joblib.load(os.path.join(hmm_folder, 'HMM-model.pkl'))
tad_ax = vis.fig.add_subplot(vis.grid[vis.track_count])
vis.track_count += 1
tad_ax.set_xlim(vis.hx.min(), vis.hx.max())
x = range(int(vis.hx.min()), int(vis.hx.max()))
work = TADcaller(vis.matrix, vis.res, hmm, window_size=50000)
work.callDomains()
y = work.DIs
d = np.zeros(len(y))
tad_ax.fill_between(x, d, y, where=y >= d, interpolate=True, color='blue')
tad_ax.fill_between(x, d, y, where=y <= d, interpolate=True, color='red')
tad_ax.set_axis_off()
vis.plot_arcs(lw=1.5, cutoff='top', gene_filter=['FOXO1','PAX3'])
#subplot4
vis.plot_signal('P3F', '/home/gdstantonlab/mxw010/Data/Ben/downsampled/data/bigwig/RH4_pval.bigwig',
                label_size=6, data_range_size=9, max_value=3, color='#6A3D9A')
#subplot5
vis.plot_signal('CTCF', '/home/gdstantonlab/mxw010/Data/Ben/CTCF/GSL-BS-2483/result/file_transfer/RH4/RH4_pooled_pval.bigwig',
                label_size=6, data_range_size=9, max_value=3, color='#6A3D9A')
#subplot6
vis.plot_genes(release=79, filter_=gene_list, fontsize=3)
#subplot7
vis.plot_chromosome_bar(name_size=11, coord_size=3)
vis.outfig(outname)
plt.close(fig=vis.fig)

#C38,10kb
res='10kb'
assembly='C38   translocation,13,40610000,+,2,222300000,+   13,39430000   2,222190000'
ID=assembly.split("   ")[0]
outname = 'RH4.' + ID + "_" +  res+ '.pdf'
vis = Triangle(clr, assembly, n_rows=8, figsize=(
    7, 4.2), track_partition=[5, 1, 1,0.8,1, 1, 1, 0.5], correct='weight')
#subplot1
vis.matrix_plot(vmin=0)
vis.plot_chromosome_bounds(linewidth=2.5)
vis.plot_loops('RH4_neoloops-10kb.txt',
                face_color='none', marker_size=40, cluster=True)
#subplot2
vis.plot_neoTAD(ws=1000000)
#subplot3, standard DI plot is buggy
hmm_folder = os.path.join(os.path.split(neoloop.__file__)[0], 'data')
hmm = joblib.load(os.path.join(hmm_folder, 'HMM-model.pkl'))
tad_ax = vis.fig.add_subplot(vis.grid[vis.track_count])
vis.track_count += 1
tad_ax.set_xlim(vis.hx.min(), vis.hx.max())
x = range(int(vis.hx.min()), int(vis.hx.max()))
work = TADcaller(vis.matrix, vis.res, hmm, window_size=50000)
work.callDomains()
y = work.DIs
d = np.zeros(len(y))
tad_ax.fill_between(x, d, y, where=y >= d, interpolate=True, color='blue')
tad_ax.fill_between(x, d, y, where=y <= d, interpolate=True, color='red')
tad_ax.set_axis_off()
vis.plot_arcs(lw=1.5, cutoff='top', gene_filter=['FOXO1','PAX3'])
#subplot4
vis.plot_signal('P3F', '/home/gdstantonlab/mxw010/Data/Ben/downsampled/data/bigwig/RH4_pval.bigwig',
                label_size=6, data_range_size=9, max_value=3, color='#6A3D9A')
#subplot5
vis.plot_signal('CTCF', '/home/gdstantonlab/mxw010/Data/Ben/CTCF/GSL-BS-2483/result/file_transfer/RH4/RH4_pooled_pval.bigwig',
                label_size=6, data_range_size=9, max_value=3, color='#6A3D9A')
#subplot6
vis.plot_genes(release=79, filter_=gene_list, fontsize=3)
#subplot7
vis.plot_chromosome_bar(name_size=11, coord_size=3)
vis.outfig(outname)
plt.close(fig=vis.fig)

# release = 97
# species = 'human'
# ref = EnsemblRelease(release, species=species)
# ID = 'C38'
# assembly = 'C38  translocation,13,40610000,+,2,222900000,+  13,39430000  2,222190000'
# contig1 = assembly.split()[2].split(",")[0]
# start1 = int(assembly.split()[1].split(",")[2])
# end1 = int(assembly.split()[2].split(",")[1])
# contig2 = assembly.split()[-1].split(",")[0]
# start2 = int(assembly.split()[1].split(",")[5])
# end2 = int(assembly.split()[-1].split(",")[1])

# if start1 > end1:
#     temp = start1
#     start1 = end1
#     end1 = temp

# if start2 > end2:
#     temp = start2
#     start2 = end2
#     end2 = temp

# gene_list = []
# for item in ref.genes_at_locus(contig1, start1, end1) + ref.genes_at_locus(contig2, start2, end2):
#     if item.is_protein_coding:
#         gene_list.append(item.gene_name)

# hmm_folder = os.path.join(os.path.split(neoloop.__file__)[0], 'data')
# hmm = joblib.load(os.path.join(hmm_folder, 'HMM-model.pkl'))


# #gene_list = ['PROSER1', 'UFM1', 'FREM2', 'LHFPL6', 'EPHA4', 'PAX3']
# vis = Triangle(clr, assembly, n_rows=8, figsize=(14, 12), track_partition=[
#                5, 0.8, 0.3, 0.8, 0.3, 0.3, 0.5, 0.5], correct='sweight')
# vis.matrix_plot(vmin=0)
# vis.plot_chromosome_bounds(linewidth=2.5)
# vis.plot_loops('RH4_neoloops-10kb.txt', face_color='none',
#                marker_size=40, cluster=True)
# # plot TADs, sub fig 2

# # P3F signal, sub fig 4
# vis.plot_signal('P3F', '/home/gdstantonlab/mxw010/Data/Ben/downsampled/data/bigwig/RH4_pval.bigwig',
#                 label_size=6, data_range_size=9, max_value=3, color='#6A3D9A')
# # ctcf, sub fig 5
# vis.plot_motif(
#     '/home/gdstantonlab/mxw010/Data/Ben/CTCF/GSL-BS-2483/result/motif/RH4/CTCF_motif_orientation.bed', subset='+')
# vis.plot_motif(
#     '/home/gdstantonlab/mxw010/Data/Ben/CTCF/GSL-BS-2483/result/motif/RH4/CTCF_motif_orientation.bed', subset='-')
# # gene annot, subfig 6
# vis.plot_genes(fontsize=10, filter_=gene_list)
# #vis.plot_genes(filter_=['PRAME','BCRP4', 'RAB36', 'BCR', 'ABL1', 'NUP214'],label_aligns={'PRAME':'right','RAB36':'right'}, fontsize=9)
# # chr, subfigure 7
# vis.plot_chromosome_bar(name_size=10, coord_size=4.8)
# vis.outfig('RH4.+ ID + '.pdf', dpi=300)

# # RH30

# ID = 'C37'
# assembly = 'C37  translocation,13,40620000,+,2,222200000,-   13,39000000   2,223460000'


# release = 97
# species = 'human'
# ref = EnsemblRelease(release, species=species)
# ID = 'C36'
# contig1 = assembly.split()[2].split(",")[0]
# start1 = int(assembly.split()[1].split(",")[2])
# end1 = int(assembly.split()[2].split(",")[1])
# contig2 = assembly.split()[-1].split(",")[0]
# start2 = int(assembly.split()[1].split(",")[5])
# end2 = int(assembly.split()[-1].split(",")[1])

# if start1 > end1:
#     temp = start1
#     start1 = end1
#     end1 = temp

# if start2 > end2:
#     temp = start2
#     start2 = end2
#     end2 = temp

# gene_list = []
# for item in ref.genes_at_locus(contig1, start1, end1) + ref.genes_at_locus(contig2, start2, end2):
#     if item.is_protein_coding:
#         gene_list.append(item.gene_name)



# gene_list = ['UFM1', 'FREM2', 'LHFPL6', 'EPHA4', 'PAX3', 'FOXO1']

# vis = Triangle(clr, assembly, n_rows=7, figsize=(14, 12), track_partition=[
#                5, 0.8, 0.3, 0.8, 0.8, 0.5, 0.5], correct='weight')
# vis.matrix_plot(vmin=0)
# vis.plot_chromosome_bounds(linewidth=2.5)
# vis.plot_loops('RH30_neoloops-10kb.txt', face_color='none',
#                marker_size=40, cluster=True)
# # plot TADs, sub fig 2


# # neoTADs, sub fig 3
# vis.plot_neoTAD(ws=1000000)
# # P3F signal, sub fig 4
# vis.plot_signal('P3F', '/home/gdstantonlab/mxw010/Data/Ben/downsampled/data/bigwig/RH30_pval.bigwig',
#                 label_size=6, data_range_size=9, max_value=3, color='#6A3D9A')
# # ctcf, sub fig 5
# vis.plot_signal('CTCF', '/home/gdstantonlab/mxw010/Data/Ben/CTCF/GSL-BS-2483/CTCF/RH30/align_primary/chip/92ea8e57-ca43-40e4-92a5-5ff47483d886/call-macs2_signal_track_pooled/execution/rep.pooled_x_ctl.pooled.pval.signal.bigwig',
#                 label_size=6, data_range_size=9, max_value=3, color='#6A3D9A')
# # gene annot, subfig 6
# vis.plot_genes(fontsize=10, filter_=gene_list)
# #vis.plot_genes(filter_=['PRAME','BCRP4', 'RAB36', 'BCR', 'ABL1', 'NUP214'],label_aligns={'PRAME':'right','RAB36':'right'}, fontsize=9)
# # chr, subfigure 7
# vis.plot_chromosome_bar(name_size=10, coord_size=4.8)
# vis.outfig('RH30.C37.pdf', dpi=300)