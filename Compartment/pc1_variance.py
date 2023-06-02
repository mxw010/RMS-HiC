import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cooler
import cooltools.lib.plotting
import cooltools
import bbi
import seaborn as sns
import glob
import os


cool_file ="/home/gdstantonlab/mxw010/Data/SPICE-C/RH4/cools/RH4.mcool::resolutions/100000"

clr = cooler.Cooler(cool_file)
# import bioframe
# bins = clr.bins()[:]
# hg38_genome = bioframe.load_fasta('/home/gdstantonlab/mxw010/GenRef/hg38/Sequence/WholeGenomeFasta/genome.fa');
# ## note the next command may require installing pysam
# gc_cov = bioframe.frac_gc(bins, hg38_genome)
# gc_cov.to_csv('hg38_gc_cov_100kb.tsv',index=False,sep='\t')

gc_cov = pd.read_csv('/home/gdstantonlab/mxw010/Data/SPICE-C/RH4/cools/hg38_gc_cov_100kb.tsv',sep="\t")




#maximum length to calcualte:
#sum(~np.isnan(clr.bins().fetch('chr21')['weight'].values))

#at least 100
# length = []
# for chr in clr.chromnames:
#     length.append(sum(~np.isnan(clr.bins().fetch(chr)['weight'].values)))


colnames = [ 'eigen' + str(i) for i in range(1,101)]
cumvar_genome = pd.DataFrame(columns = colnames)

view_df = pd.DataFrame({'chrom': clr.chromnames,
                    'start': 0,
                    'end': clr.chromsizes.values,
                    'name': clr.chromnames}
                    )  
k=0
for i in range(k,len(clr.chromnames)):
    chr=clr.chromnames[i]
    print("i=", i, "chr=", chr)
    # bins = clr.bins().fetch(chr)
    # bins.index = bins.index-clr.offset(chr)
    # pixels = clr.matrix(as_pixels=True).fetch(chr,chr)
    # pixels.loc[:,'bin1_id'] = pixels.loc[:,'bin1_id'] - clr.offset(chr)
    # pixels.loc[:,'bin2_id'] = pixels.loc[:,'bin2_id'] - clr.offset(chr)
    # cooler.create_cooler("test.cool", bins=bins, pixels=pixels)
    #clr=cooler.Cooler("test.cool")
    chr_df = view_df.iloc[i:i+1,:].copy()
    chr_df.index=[0]
    length = sum(~np.isnan(clr.bins().fetch(chr)['weight'].values))
    cis_eigs = cooltools.eigs_cis(
                            clr,
                            gc_cov,
                            view_df=chr_df,
                            n_eigs=length,
                            phasing_track_col='GC',
                            )
    x = cis_eigs[0].iloc[:,4:].to_numpy()
    y = x[~np.isnan(x)]
    cumvariance = np.cumsum(y**2)/np.sum(y**2)
    cumvar_genome.loc[i] = cumvariance[0:100]

cumvar_genome.index = clr.chromnames[k:]
cumvar_genome.to_csv('variance_explaned_top_100.csv')

cumvar_genome.loc[:,'chrom'] = cumvar_genome.index

# cum_index = cumvar_genome.columns.values
# cum_index[0] = 'chrom'
# cumvar_genome.columns = cum_index

df =pd.melt(cumvar_genome, id_vars='chrom')
df.loc[:,'value'] = df['value'] *100
plt.style.use('seaborn-whitegrid')
#barplot 
fig, ax = plt.subplots(1, figsize=(20, 10))
ax = sns.barplot(data=df.query("variable == 'eigen1'"),x='chrom',y='value', color='grey')
ax.set_xlabel("Chromosomes", fontsize=22)
ax.set_ylabel("Variance Explained (%)", fontsize=22)
ax.set_title("Percentage of variance explained by PC1 (RH30)", size=26)
ax.xaxis.set_tick_params(labelsize=18, rotation=90)
ax.yaxis.set_tick_params(labelsize=18)
fig.savefig("Rh30_1.pdf")

#sum(~np.isnan(clr.bins().fetch('chr21')['weight'].values)),
#cooltools algorithm
# eigvecs, eigvals = numutils.get_eig(OE, n_eigs, mask_zero_rows=True)
# numutils.get_eig:
# _eigvals, _eigvecs = scipy.sparse.linalg.eigsh(mat, _n)
# order = np.argsort(-np.abs(_eigvals))
# eigvals[:_n] = _eigvals[order]
# eigvecs[:_n,:] = _eigvecs.T[order]
# eigvecs /= np.sqrt(np.nansum(eigvecs ** 2, axis=1))[:, None]
# eigvecs *= np.sqrt(np.abs(eigvals))[:, None]


# cov =  cis_eigs[0].iloc[:,4:24].astype(float)**2
# cov.index = cis_eigs[0].chrom

# cumvar = cov.copy()

# for chrom, row in cov.iterrows():
#     cumvar.loc[chrom] = [ sum(row[0:idx+1])/sum(row) for idx,x in enumerate(row) ]
    
# cumvar.loc[:,'chrom'] = cov.index
# df =pd.melt(cumvar, id_vars='chrom')
# df.loc[:,'value'] = df['value'] *100
# plt.style.use('seaborn-whitegrid')
# fig, ax = plt.subplots(1)
# #line plot of all eigenvalues
# ax = sns.lineplot(data=df, x='variable', y='value', style='chrom')
# fig.show()




# plt.suptitle("TFs Positioning at TAD Boundaries," +
#             cell + " Cell Line At " + res, fontsize=20)
# fig.savefig("test.pdf")



    
