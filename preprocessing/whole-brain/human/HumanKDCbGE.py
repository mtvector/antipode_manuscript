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
import bbknn
#Load my pipeline functions
import importlib
import importlib.util
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../utils/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)
figdir='/wynton/group/ye/mtschmitz/figures/PanHumanGE/'
sc.settings.figdir=figdir
scv.settings.figdir=figdir
sc.settings.file_format_figs='pdf'
sc.settings.autosave=True
sc.settings.autoshow=False


def most_frequent(List): 
    return max(set(List), key = List.count)

def variance(X,axis=1):
    return( (X.power(2)).mean(axis=axis)-(np.power(X.mean(axis=axis),2)) )

filepath='/wynton/group/ye/mtschmitz/humanmotor/cathg38_gencodev33_kallisto'
files=os.listdir(filepath)
files=[f for f in files if 'kOut' in f]
fileList=files
#regionlist=['mix9','mix10','ge','mix8','amy','str','septum','otor','dfc','cing','pfc','parietal','v1','temp','somato','re-optic','poa','mix7','mix5','insula','hippo','hip','hippocampus']
#fileList=[f for f in files if any([x in f.lower() for x in regionlist]) ]
#tplist=['E40','E50','E65','E80','E90','E100','PEC_Yale','Mac2','Mix']
#fileList=[f for f in fileList if any([x in f for x in tplist]) ]
print(fileList,flush=True)
newfile='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityPanHumanPresupervise.h5ad'
min_genes=800

adata=sc.read(re.sub('\.h5ad','Processed.h5ad',newfile))

dlxmean=[]
for c in adata.obs.leiden.unique():
    dlxmean.append(np.nansum(adata[adata.obs.leiden==c,['DLX1','DLX2','DLX5','DLX6','GAD1','GAD2']].X,axis=1).mean())
adata=adata[adata.obs.leiden.isin(adata.obs.leiden.unique()[np.array(dlxmean)>np.mean(dlxmean)]),:]
print(adata,flush=True)
print(adata.X.shape,flush=True)
print(adata.raw.X.shape,flush=True)
print(adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)].shape,flush=True)
adata.var_names_make_unique()
adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)]
if True:#not os.path.exists(re.sub('\.h5ad','Processed.h5ad',newfile)):
    ribo_genes=[name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL') ]
    adata.obs['percent_ribo'] = np.sum(
    adata[:, ribo_genes].X, axis=1) / np.sum(adata.X, axis=1)
    print(adata)
    #adata=adata[adata.obs['percent_ribo']<.5,:]
    #adata=adata[adata.obs['percent_ribo']<.2,:]
    adata._inplace_subset_obs(adata.obs['percent_ribo']<.4)
    print(adata)
    sc.pp.filter_genes(adata,min_cells=10)
    #sc.pp.filter_cells(adata,min_counts=1000)
    sc.pp.filter_cells(adata,min_genes=min_genes)
    keep_cat=adata.obs['batch_name'].value_counts().index[adata.obs['batch_name'].value_counts()>2]
    adata=adata[adata.obs.batch_name.isin(keep_cat),:]
    print(adata)
    sc.pp.normalize_total(adata,exclude_highly_expressed=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata,n_top_genes=5000,subset=False)
    sc.pp.scale(adata,max_value=10)
    sc.pp.pca(adata,n_comps=100)
    #sc.pp.neighbors(adata)
    bbknn.bbknn(adata,batch_key='batch_name',n_pcs=100,neighbors_within_batch=2)
    sc.tl.leiden(adata,resolution=10)
    sc.tl.umap(adata,spread=2)
    adata.write(re.sub('\.h5ad','Processed.h5ad',newfile))
    #adata=sc.AnnData(adata.raw.X,var=adata.raw.var,obsm=adata.obsm,uns=adata.uns,obs=adata.obs)
    #adata.raw=adata
    #sc.pp.normalize_per_cell(adata)
    #sc.pp.log1p(adata)
    #sc.pp.scale(adata,max_value=10)
    sc.pl.umap(adata,color=['leiden'],legend_loc='on data',save='leiden')
    sc.pl.umap(adata,color=['batch_name'],save='batch_name')
    sc.pl.umap(adata,color=['region'],save='region')
    sc.pl.umap(adata,color=['timepoint'],save='timepoint')
    easygenes= ['MKI67','SOX2','AQP4','EDNRB','IL33','PDGFRA','HOPX','ERBB4','CALB1','CALB2','GAD1','GAD2','GADD45G','LHX6','SST','NKX2-1','MAF','SP8','SP9','PROX1','VIP','CCK','NPY','LAMP5','HTR3A','NR2F1','NR2F2','TOX3','ETV1','SCGN','FOXP1','FOXP2','FOXP4','TH','PAX6','MEIS2','ISL1','PENK','CRABP1','ZIC1','ZIC2','EBF1','DLX5','GSX2','RBFOX3','PPP1R17','EOMES','DCX','TUBB3','BCL11B','TLE4','FEZF2','SATB2','TBR1','AIF1','RELN','PECAM1','HBZ','ROBO1','ROBO2','ROBO3','ROBO4','FOXG1','EMX1','DLX1','DLX2','DLX5','POMC']
    easygenes=[x for x in easygenes if x in adata.var.index]
    sc.pl.umap(adata,color=easygenes,use_raw=False,save='FaveGenes')
    INgenes=['MKI67','SOX2','HOPX','DCX','FOXG1','GSX2','ERBB4','CALB1','CALB2','RSPO3','RELN','DLX2','DLX5','GAD1','GAD2','GADD45G','LHX6','LHX7','LHX8','SST','NKX2-1','MAF','CHODL','SP8','NR2F2','PROX1','VIP','CCK','NPY','PAX6','ETV1','SCGN','ROBO2','CASZ1','TSHZ1','OPRM1','FOXP1','FOXP2','FOXP4','VAX1','TH','DDC','SLC18A2','MEIS2','LAMP5','TLE4','ISL1','DRD1','PENK','ADORA2A','TAC1','TAC3','CRABP1','ANGPT2','RBP4','STXBP6','ZIC1','ZIC4','EBF1','SATB2']
    INgenes=[x for x in INgenes if x in adata.var.index]
    sc.pl.umap(adata,color=INgenes,use_raw=False,save='INGenes')

    hierarchy_key='leiden'
    rgs=sc.tl.rank_genes_groups(adata,groupby=hierarchy_key,method='logreg',use_raw=False,copy=True).uns['rank_genes_groups']#,penalty='elasticnet',solver='saga')#or penalty='l1'
    result=rgs
    groups = result['names'].dtype.names
    df=pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores']})
    df.to_csv(os.path.join(sc.settings.figdir,"LogReg"+hierarchy_key+"Norm.csv"))
    topgenes=df.iloc[0:8,['_n' in x for x in df.columns]].T.values
    cols=df.columns[['_n' in x for x in df.columns]]
    cols=[re.sub('_n','',x) for x in cols]
    topdict=dict(zip(cols,topgenes))
    sc.tl.dendrogram(adata,groupby=hierarchy_key)
    var_dict=dict(zip(adata.uns["dendrogram_['"+hierarchy_key+"']"]['categories_ordered'],[topdict[x] for x in adata.uns["dendrogram_['"+hierarchy_key+"']"]['categories_ordered']]))
    sc.pl.matrixplot(adata,groupby=hierarchy_key,var_names=var_dict,cmap='RdBu_r',use_raw=False,dendrogram=True,save='leiden_marker')



 


