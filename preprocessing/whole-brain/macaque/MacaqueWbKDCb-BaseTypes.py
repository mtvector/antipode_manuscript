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
sc.settings.file_format_figs='pdf'
sc.settings.autosave=True
sc.settings.autoshow=False

min_genes=800

filepath='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/'
files=os.listdir(filepath)
files=[f for f in files if 'MacWb' in f]

for f in reversed(files):
    print(f, flush=True)
    adata=sc.read(os.path.join(filepath,f))
    base_type=re.split('_',f)[0]
    figdir='/wynton/group/ye/mtschmitz/WBfigures/macaqueWb_'+base_type
    sc.settings.figdir=figdir
    scv.settings.figdir=figdir
    adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)]
    print(adata)
    print(adata.X, flush=True)
    adata.obs.drop('n_genes',axis=1,inplace=True)
    adata.var.drop('n_cells',axis=1,inplace=True)
    sc.pp.filter_cells(adata,min_genes=800)
    sc.pp.filter_genes(adata,min_cells=10)
    sufficient_cells=adata.obs['batch_name'].value_counts().index[adata.obs['batch_name'].value_counts()>4]
    adata=adata[adata.obs['batch_name'].isin(sufficient_cells),:]
    adata.obs['batch_name'].cat.remove_unused_categories(inplace=True)
    sc.pp.normalize_total(adata,exclude_highly_expressed=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata,n_top_genes=6000,batch_key='batch_name',subset=False)
    sc.pp.scale(adata,max_value=10)
    sc.pp.pca(adata,n_comps=50)
    bbknn.bbknn(adata,batch_key='batch_name',n_pcs=50,neighbors_within_batch=3)
    sc.tl.leiden(adata,resolution=1.5)
    sc.tl.umap(adata,spread=2)
    sc.pl.umap(adata,color=['leiden'],save='leiden',legend_loc='on data')
    sc.pl.umap(adata,color=['batch_name'],save='batch_name')
    sc.pl.umap(adata,color=['region'],save='region')
    sc.pl.umap(adata,color=['timepoint'],save='timepoint')
    easygenes= ['MKI67','SOX2','AQP4','EDNRB','IL33','PDGFRA','HOPX','ERBB4','CALB1','CALB2','GAD1','GAD2','GADD45G','LHX6','SST','NKX2-1','MAF','SP8','SP9','PROX1','VIP','CCK','NPY','LAMP5','HTR3A','NR2F1','NR2F2','TOX3','ETV1','SCGN','FOXP1','FOXP2','FOXP4','FOXG1','TH','PAX6','MEIS2','SHTN1','ISL1','PENK','CRABP1','ZIC1','ZIC2','EBF1','DLX5','GSX2','RBFOX3','PPP1R17','EOMES','DCX','TUBB3','BCL11B','TLE4','FEZF2','SATB2','TBR1','AIF1','RELN','PECAM1','HBZ']
    easygenes=[x for x in easygenes if x in adata.var.index]
    sc.pl.umap(adata,color=easygenes,use_raw=False,save='FaveGenes')
    hierarchy_key='leiden'
    rgs=sc.tl.rank_genes_groups(adata,groupby=hierarchy_key,method='logreg',use_raw=False,copy=True).uns['rank_genes_groups']#,penalty='elasticnet',solver='saga')#or penalty='l1'
    result=rgs
    groups = result['names'].dtype.names
    df=pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores']})
    df.to_csv(os.path.join(sc.settings.figdir,"LogReg"+hierarchy_key+"Norm.csv"))
    topgenes=df.iloc[0:4,['_n' in x for x in df.columns]].T.values
    cols=df.columns[['_n' in x for x in df.columns]]
    cols=[re.sub('_n','',x) for x in cols]
    topdict=dict(zip(cols,topgenes))
    sc.tl.dendrogram(adata,groupby=hierarchy_key)
    var_dict=dict(zip(adata.uns["dendrogram_['"+hierarchy_key+"']"]['categories_ordered'],[topdict[x] for x in adata.uns["dendrogram_['"+hierarchy_key+"']"]['categories_ordered']]))
    sc.pl.matrixplot(adata,groupby=hierarchy_key,var_names=var_dict,save='top_degenes',cmap='RdBu_r',use_raw=False,dendrogram=True)



