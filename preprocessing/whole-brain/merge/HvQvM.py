import os
import scanpy
import anndata
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from scipy import stats
import re
import sklearn
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
from collections import Counter
import random
import seaborn
import sys
import shutil
import scvelo as scv
import tqdm
import bbknn
#Load my pipeline functions
import importlib
import importlib.util
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../utils/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)
sc.settings.figdir='/wynton/group/ye/mtschmitz/WBfigures/HvQvM_dev/'
scv.settings.figdir='/wynton/group/ye/mtschmitz/WBfigures/HvQvM_dev/'
sc.settings.file_format_figs='pdf'
sc.settings.autosave=True
sc.settings.autoshow=False

mdata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT_chang_h5ad/KDCbVelocityMouseWbPresupervision.h5ad')
mdata.obs['species']='mouse'
mdata.var.index=[x.upper() for x in mdata.var.index]
mdata=sc_utils.sum_duplicate_var(mdata)
#mdata.X=mdata.raw.X[:,mdata.raw.var.index.isin(mdata.var.index)]#.todense()
#mdata.X=scipy.sparse.csr_matrix(mdata.X)

orthos=pd.read_csv('/wynton/home/ye/mschmitz1/utils/HOM_AllOrganism.rpt',sep='\t')
orthos=orthos.loc[orthos['NCBI Taxon ID'].isin([10090,9606]),:]
classcounts=orthos['DB Class Key'].value_counts()
one2one=classcounts.index[list(classcounts==2)]
orthos=orthos.loc[orthos['DB Class Key'].isin(one2one),:]

htab=orthos.loc[orthos['NCBI Taxon ID']==9606,:]
mtab=orthos.loc[orthos['NCBI Taxon ID']==10090,:]
genemapping=dict(zip([x.upper() for x in mtab['Symbol']],htab['Symbol']))
mdata.var.index=[x.upper() for x in mdata.var.index]
mdata=mdata[:,mdata.var.index.isin(genemapping.keys())]
mdata.var.index=[genemapping[x] for x in mdata.var.index]
#mdata=mdata[np.random.choice(mdata.obs.index,300000,replace=False),:]

qdata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT_chang_h5ad/KDCbVelocityMacaqueWbPresupervise.h5ad')
qdata.obs['species']='macaque'
#qdata.X=qdata.raw.X[:,qdata.raw.var.index.isin(qdata.var.index)]#.todense()
#qdata.X=scipy.sparse.csr_matrix(qdata.X)
#qdata=qdata[np.random.choice(qdata.obs.index,300000,replace=False),:]
qdata=sc_utils.sum_duplicate_var(qdata)

hdata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT_chang_h5ad/KDCbVelocityPanHumanPresupervise.h5ad')
hdata.obs['species']='human'
#hdata.X=hdata.raw.X[:,hdata.raw.var.index.isin(hdata.var.index)]#.todense()
#hdata.X=scipy.sparse.csr_matrix(hdata.X)
#hdata=hdata[np.random.choice(hdata.obs.index,300000,replace=False),:]
hdata=sc_utils.sum_duplicate_var(hdata)

qdata.var_names_make_unique()
mdata.var_names_make_unique()
hdata.var_names_make_unique()
qdata.obs_names_make_unique()
mdata.obs_names_make_unique()
hdata.obs_names_make_unique()

adata=anndata.AnnData.concatenate(mdata,qdata,hdata,batch_key='species_batch')
adata.raw=adata
print(adata)
print(adata.X)
print(adata.obs)
print(adata.var)

ribo_genes=[name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL') ]
adata.obs['percent_ribo'] = np.sum(
adata[:, ribo_genes].X, axis=1) / np.sum(adata.X, axis=1)
mito_genes = [name for name in adata.var_names if name in ['ND1','ND2','ND4L','ND4','ND5','ND6','ATP6','ATP8','CYTB','COX1','COX2','COX3'] or name.startswith('CHRM-') or name.startswith('MT-')]
adata.obs['percent_mito'] = np.sum(
adata[:, adata.var.index.isin(mito_genes)].X, axis=1) / np.sum(adata.X, axis=1)
print(adata)
print((adata.obs['percent_ribo'].to_numpy()<0.40)[0])
adata=adata[adata.obs.index[(adata.obs['percent_ribo'].to_numpy()<0.4)[0]],:]
adata=adata[adata.obs.index[(adata.obs['percent_mito'].to_numpy()<0.15)[0]],:]
adata=adata[:,~adata.var.index.isin(mito_genes)]

sc.pp.filter_genes(adata,min_cells=10)
sc.pp.normalize_total(adata,exclude_highly_expressed=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=2000,batch_key='species')
#markers=pd.read_csv(os.path.expanduser('~/markers.txt'),sep='\t')
##adata=adata[:,adata.var.index.isin(markers['hgnc_symbol'])]
#adata.var['highly_variable']=list(adata.var.index.isin(markers['hgnc_symbol']))
adata.var['highly_variable']=(adata.var['highly_variable']&(adata.var['highly_variable_nbatches']>1))
sc.pp.scale(adata,max_value=10)
sc.pp.pca(adata,n_comps=100)
bbknn.bbknn(adata,batch_key='species',n_pcs=100,neighbors_within_batch=9)
sc.tl.leiden(adata,resolution=1.5)#.8 works for cross species when not marker subsetted
sc.tl.umap(adata)
adata.write('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT_chang_h5ad/HvQvMall.h5ad')
sc.pl.umap(adata,color=['leiden'],save='leiden')
sc.pl.umap(adata,color=['species'],save='species')
sc.pl.umap(adata,color=['batch_name'],save='batch_name')
sc.pl.umap(adata,color=['region'],save='region')
sc.pl.umap(adata,color=['timepoint'],save='timepoint')
viewgenes= ['MKI67','SOX2','AQP4','EDNRB','IL33','PDGFRA','HOPX','ERBB4','FOXG1','CALB1','CALB2','GAD1','GAD2','GADD45G','LHX1','LHX5','LHX6','LHX8','SST','NKX2-1','MAF','SP8','SP9','PROX1','VIP','CCK','NPY','LAMP5','NR2F1','NR2F2','ETV1','SCGN','FOXP1','FOXP2','FOXP4','TH','DDC','SLC18A2','PAX6','MEIS2','ISL1','PENK','CRABP1','ZIC1','ZIC4','EBF1','DLX5','TFAP2A','RBFOX3','PPP1R17','EOMES','DCX','TUBB3','BCL11B','TLE4','FEZF2','SATB2','TBR1','SLC17A6','SLC17A7','RORB','AIF1','RELN','PECAM1','HBZ']
easygenes=[x for x in viewgenes if x in adata.var.index]
sc.pl.umap(adata,color=viewgenes,use_raw=False,save='FaveGenes')



