#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Taken from sccloud
import numpy as np
import pandas as pd
import time
from natsort import natsorted
import multiprocessing
from sklearn.cluster import KMeans
import scanpy as sc
import os
import sys


def estimate_background_probs(adt, random_state = 0):
	adt.obs['counts'] = adt.X.sum(axis = 1).A1 if adt.shape[1] > 1 else adt.X
	counts_log10 = np.log10(adt.obs['counts'].values.reshape(-1, 1))
	kmeans = KMeans(n_clusters = 2, random_state = random_state).fit(counts_log10)
	signal = 0 if kmeans.cluster_centers_[0] > kmeans.cluster_centers_[1] else 1
	adt.obs['hto_type'] = 'background'
	adt.obs.loc[kmeans.labels_ == signal, 'hto_type'] = 'signal'

	idx = np.isin(adt.obs['hto_type'], 'background')
	pvec = adt.X[idx,].sum(axis = 0).A1 if adt.shape[1] > 1 else np.array(adt.X[idx,].sum())
	pvec /= pvec.sum()

	adt.uns['background_probs'] = pvec



def estimate_probs(arr, pvec, alpha, alpha_noise, tol):
	probs = np.zeros(pvec.size + 1)
	old_probs = np.zeros(pvec.size + 1)
	z = np.zeros(pvec.size + 1)
	noise = pvec.size
	# Estimate MLE without Generalized Dirichlet prior
	probs_mle = arr / arr.sum()
	probs[noise] = (probs_mle / pvec).min() + 0.01
	probs[:-1] = np.maximum(probs_mle - probs[noise] * pvec, 0.01)
	probs = probs / probs.sum()

	# EM algorithm
	i = 0
	eps = 1.0
	while eps > tol:
		i += 1
		old_probs[:] = probs[:]
		# E step
		z[:-1] = alpha - 1.0
		z[noise] = alpha_noise - 1.0
		for j in range(pvec.size):
			if arr[j] > 0:
				p = probs[j] / (probs[noise] * pvec[j] + probs[j])
				z[j] += arr[j] * p
				z[noise] += arr[j] * (1.0 - p)
		# M step
		idx = z > 0.0
		probs[idx] = z[idx] / z[idx].sum()
		probs[~idx] = 0.0
		eps = np.linalg.norm(probs - old_probs, ord = 1)
		# print ("i = {}, eps = {:.2g}.".format(i, eps))

	return probs



def get_droplet_info(probs, sample_names):
	ids = np.nonzero(probs >= 0.1)[0]
	ids = ids[np.argsort(probs[ids])[::-1]]
	return ('singlet' if ids.size == 1 else 'doublet',
			','.join([sample_names[i] for i in ids]))

def calc_demux(data, adt, nsample, min_signal, probs = 'raw_probs'):
	demux_type = np.full(data.shape[0], 'unknown', dtype = 'object')
	assignments = np.full(data.shape[0], '', dtype = 'object')

	signals = adt.obs['counts'].reindex(data.obs_names, fill_value = 0.0).values * (1.0 - data.obsm[probs][:,nsample])
	idx = signals >= min_signal

	tmp = data.obsm[probs][idx,]
	norm_probs = tmp[:,0:nsample] / (1.0 - tmp[:,nsample])[:,None]

	values1 = []
	values2 = []
	for i in range(norm_probs.shape[0]):
		droplet_type, droplet_id = get_droplet_info(norm_probs[i,], adt.var_names)
		values1.append(droplet_type)
		values2.append(droplet_id)

	demux_type[idx] = values1
	data.obs['demux_type'] = pd.Categorical(demux_type, categories = ['singlet', 'doublet', 'unknown'])
	assignments[idx] = values2
	data.obs['assignment'] = pd.Categorical(assignments, categories = natsorted(np.unique(assignments)))



def demultiplex(data, adt, min_signal = 10.0, alpha = 0.0, alpha_noise = 1.0, tol = 1e-6, n_threads = 1):
	start = time.time()
	
	nsample = adt.shape[1]
	data.uns['background_probs'] = adt.uns['background_probs']

	idx_df = data.obs_names.isin(adt.obs_names)
	adt.obs['rna_type'] = 'background'
	adt.obs.loc[data.obs_names[idx_df], 'rna_type'] = 'signal'

	if nsample == 1:
		print("Warning: detect only one barcode, no need to demultiplex!")
		data.obsm['raw_probs'] = np.zeros((data.shape[0], nsample + 1))
		data.obsm['raw_probs'][:, 0] = 1.0
		data.obsm['raw_probs'][:, 1] = 0.0
		data.obs['demux_type'] = 'singlet'
		data.obs['assignment'] = adt.var_names[0]
	else:
		if nsample == 2:
			print("Warning: detect only two barcodes, demultiplexing accuracy might be affected!")

		ncalc = idx_df.sum()
		if ncalc < data.shape[0]:
			nzero = data.shape[0] - ncalc
			print("Warning: {} cells do not have ADTs, percentage = {:.2f}%.".format(nzero, nzero * 100.0 / data.shape[0]))
		adt_small = adt[data.obs_names[idx_df],].X.toarray()

		data.obsm['raw_probs'] = np.zeros((data.shape[0], nsample + 1))
		data.obsm['raw_probs'][:, nsample] = 1.0

		iter_array = [(adt_small[i,], adt.uns['background_probs'], alpha, alpha_noise, tol) for i in range(ncalc)]
		with multiprocessing.Pool(n_threads) as pool:
			data.obsm['raw_probs'][idx_df, :] = pool.starmap(estimate_probs, iter_array)

		calc_demux(data, adt, nsample, min_signal)
	
	end = time.time()
	print("demuxEM.demultiplex is finished. Time spent = {:.2f}s.".format(end - start))


# In[2]:


multiseqdata=sc.read_10x_mtx('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_kallisto/E90-2019_Multi-seq_kOut/multiseq_outs',gex_only=False)


# In[3]:


adata=sc.read_10x_mtx('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_kallisto/E90-2019_Multi-seq_kOut/all/')


# In[4]:


adata=adata[adata.obs.index.isin(multiseqdata.obs.index),:]
multiseqdata=multiseqdata[multiseqdata.obs.index.isin(adata.obs.index),:]


# In[5]:


estimate_background_probs(multiseqdata)
demultiplex(adata,multiseqdata,min_signal=5.0,alpha=0.1,alpha_noise=1.0)
print(adata.obs['assignment'].value_counts())
adata.obs.loc[:,['demux_type','assignment']]


# In[6]:


min_genes=700


# In[8]:


adata.raw=adata
adata.var_names_make_unique()
ribo_genes=[name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL') ]
adata.obs['percent_ribo'] = np.sum(
adata[:, ribo_genes].X, axis=1) / np.sum(adata.X, axis=1)
mito_genes = [name for name in adata.var_names if name in ['ND1','ND2','ND4L','ND4','ND5','ND6','ATP6','ATP8','CYTB','COX1','COX2','COX3'] or name.startswith('chrM-') or name.startswith('MT-')]
adata.obs['percent_mito'] = np.sum(
adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
print(adata)


#adata=adata[adata.obs['percent_ribo']<.5,:]
#adata=adata[adata.obs['percent_ribo']<.2,:]
adata._inplace_subset_obs(adata.obs['percent_ribo']<.4)
adata._inplace_subset_obs(adata.obs['percent_mito']<.15)
adata._inplace_subset_var(~adata.var.index.isin(mito_genes))
print(adata)
sc.pp.filter_genes(adata,min_cells=10)
#sc.pp.filter_cells(adata,min_counts=1000)
sc.pp.filter_cells(adata,min_genes=min_genes)
sc.pl.highest_expr_genes(adata, n_top=20, )
print(adata)
sc.pp.normalize_total(adata,exclude_highly_expressed=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=6000,subset=False)
sc.pp.scale(adata,max_value=10)
sc.pp.pca(adata,n_comps=100)
sc.pp.neighbors(adata)
#bbknn.bbknn(adata,batch_key='batch_name',n_pcs=100,neighbors_within_batch=3)


# In[9]:


sc.tl.leiden(adata,resolution=4)
sc.tl.umap(adata,spread=2)
sc.pl.umap(adata,color=['leiden'])


# In[17]:


adata.obs.loc[adata.obs['demux_type']=='singlet','region']=adata.obs['assignment'][adata.obs['demux_type']=='singlet']


# In[19]:


sc.pl.umap(adata,color=['region'])


# In[ ]:


for alpha in [0.0,.5,1,5]:
    for alpha_noise in [.5,1,4,10]:
        print(alpha,alpha_noise)
        demultiplex(adata,multiseqdata,min_signal=15.0,alpha=alpha,alpha_noise=alpha_noise)
        try:
            print(adata.obs['assignment'].value_counts())
            print(adata.obs['assignment'].value_counts()['E90-2019_CGE'])
            print(adata.obs['assignment'].value_counts()['E90-2019_Septum'])
        except:
            pass


# In[ ]:


adata.obs['demux_type'].value_counts()


# In[50]:


import os
multiseq=pd.read_csv(os.path.join('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_kallisto/E80-2019_Multi-seq_kOut','multiseq_calls.csv'),index_col=0)


# In[53]:


multiseq=multiseq.loc[multiseq['demux_type']=='singlet',:]


# In[54]:


[sc_utils.region_format_macaque(x) for x in multiseq['assignment']]


# In[3]:


df=pd.read_csv('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_kallisto/E90-2019_Multi-seq_kOut/MULTIseq_counts.txt',sep='\t')


# In[9]:


df['E90-2019_Parietal']
seaborn.distplot(np.log(df['E90-2019_Septum']+1))


# In[10]:


seaborn.distplot(np.log(df['E90-2019_Parietal']+1))


# In[11]:


seaborn.distplot(np.log(df['E90-2019_PFC']+1))


# In[12]:


seaborn.scatterplot(np.log(df['E90-2019_Parietal']+1),np.log(df['E90-2019_Septum']+1))


# In[ ]:




