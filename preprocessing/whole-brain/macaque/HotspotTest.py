#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
import tqdm
#Load my pipeline functions
import importlib
import importlib.util
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.expanduser("../../utils/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)
sc.settings.figdir='/wynton/group/ye/mtschmitz/figures/macWbSupervise/'
scv.settings.figdir='/wynton/group/ye/mtschmitz/figures/macWbSupervise/'
sc.settings.file_format_figs='pdf'
sc.settings.autosave=False
sc.settings.autoshow=True


# In[154]:


#newfile='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityMacaqueWbPresuperviseProcessed.h5ad'
newfile='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityPanHumanPresuperviseProcessed.h5ad'
newfile='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityMouseWbPresupervisionhb_subset.h5ad'

adata=sc.read(newfile)


# In[155]:


adata


# In[156]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pl.umap(adata,color='leiden', legend_loc='on data')
sc.pl.umap(adata,color=['AIF1','FOXG1'],use_raw=False)


# In[157]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pl.umap(adata,color='region')


# In[158]:


import hotspot


# In[159]:


adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)]


# In[128]:


sc.pp.filter_genes(adata,min_cells=500,inplace=True)
sc.pp.filter_genes_dispersion(adata,n_top_genes=10000,flavor='seurat_v3')


# In[129]:


adata


# In[181]:


adata=adata[0:2000,:]


# In[194]:


hopot=hotspot.Hotspot(counts=pd.DataFrame(adata.X.todense().T),latent=pd.DataFrame(adata.obsm['X_pca']),umi_counts=pd.Series(adata.X.sum(1).A1))


# In[162]:


nz=adata.obsp['distances'].nonzero()
ncols=pd.DataFrame(nz[0]).value_counts().max()
df=pd.DataFrame(index=set(nz[0]),columns=list(range(ncols)))
count=0
this_x=None
for x,y in tqdm.tqdm(zip(nz[0],nz[1])):
    if x==this_x:
        count+=1
    else:
        count=0
    df.loc[x,count]=y
    this_x=x


# In[163]:


hopot.neighbors=df


# In[164]:


nz=adata.obsp['distances'].nonzero()
ncols=pd.DataFrame(nz[0]).value_counts().max()
df=pd.DataFrame(index=set(nz[0]),columns=list(range(ncols)))
count=0
this_x=None
for x,y in tqdm.tqdm(zip(nz[0],nz[1])):
    if x==this_x:
        count+=1
    else:
        count=0
    df.loc[x,count]=adata.obsp['distances'][x,y]
    this_x=x


# In[165]:


hopot.weights=df


# In[183]:


hopot.create_knn_graph()


# In[193]:


hopot.neighbors.astype(int).dtypes


# In[192]:


hopot.weights.astype(float).dtypes


# In[185]:


hopot.graph


# In[191]:


hopot.tree


# In[49]:


#hopot.neighbors=adata.obsp['connectivities'].nonzero()
#hopot.weights=adata.obsp['connectivities']


# In[150]:


hopot.neighbors.dtypes


# In[151]:


hopot.weights.dtypes


# In[176]:


hopot.compute_hotspot(jobs=1)


# In[65]:


adata.obsp['connectivities'].nonzero()


# In[69]:


np.where((adata.obsp['connectivities']!=0).todense())


# In[76]:


hopot.neighbors


# In[113]:


hopot.weights


# In[114]:


ncols


# In[115]:


adata.obsp['distances']


# In[116]:


ncols=pd.DataFrame(nz[0]).value_counts().max()


# In[119]:


ncols


# In[118]:


df


# In[ ]:




