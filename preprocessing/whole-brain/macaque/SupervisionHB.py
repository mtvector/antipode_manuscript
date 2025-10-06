#!/usr/bin/env python
# coding: utf-8

# In[2]:


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
sc.settings.figdir='/wynton/group/ye/mtschmitz/WBfigures/macWbSuperviseHB/'
scv.settings.figdir='/wynton/group/ye/mtschmitz/WBfigures/macWbSuperviseHB/'
sc.settings.file_format_figs='pdf'
sc.settings.autosave=False
sc.settings.autoshow=True


# In[3]:


newfile='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityMacaqueWbPresupervisehb_subset.h5ad'
newfile='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityPanHumanPresupervisehb_subset.h5ad'
newfile='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityMouseWbPresupervisionhb_subset.h5ad'

adata=sc.read(newfile)


# In[17]:


adata


# In[9]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pl.umap(adata,color='leiden', legend_loc='on data')
sc.pl.umap(adata,color=['region'],use_raw=False)


# In[10]:


#sc.tl.leiden(adata)
sc.tl.umap(adata)


# In[22]:


sc.pl.umap(adata,color=['MAP2','SLC17A6','SLC17A7','GAD2','AQP4','VIM'],use_raw=False)


# In[11]:


sc.pl.umap(adata,color='leiden', legend_loc='on data')
sc.pl.umap(adata,color=['region'],use_raw=False)


# In[12]:


newfile='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityMacaqueWbPresupervisehb_subset.h5ad'

qdata=sc.read(newfile)


# In[13]:


qdata


# In[14]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pl.umap(qdata,color='leiden', legend_loc='on data')
sc.pl.umap(qdata,color=['region'],use_raw=False)


# In[15]:


#sc.tl.leiden(qdata)
sc.tl.umap(qdata)


# In[16]:


sc.pl.umap(qdata,color='leiden', legend_loc='on data')
sc.pl.umap(qdata,color=['region'],use_raw=False)


# In[21]:


sc.pl.umap(qdata,color=['MAP2','SLC17A6','SLC17A7','GAD2','AQP4','VIM'],use_raw=False)


# In[25]:


qdata=qdata[qdata.obs.leiden.isin(qdata.obs.leiden.value_counts().index[qdata.obs.leiden.value_counts()>10]),:]
hierarchy_key='leiden'
rgs=sc.tl.rank_genes_groups(qdata,groupby=hierarchy_key,method='logreg',use_raw=False,copy=True).uns['rank_genes_groups']#,penalty='elasticnet',solver='saga')#or penalty='l1'
result=rgs
groups = result['names'].dtype.names
df=pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores']})
df.to_csv(os.path.join(sc.settings.figdir,"LogReg"+hierarchy_key+"Norm.csv"))
topgenes=df.iloc[0:4,['_n' in x for x in df.columns]].T.values
cols=df.columns[['_n' in x for x in df.columns]]
cols=[re.sub('_n','',x) for x in cols]
topdict=dict(zip(cols,topgenes))
sc.tl.dendrogram(qdata,groupby=hierarchy_key)
var_dict=dict(zip(qdata.uns["dendrogram_['"+hierarchy_key+"']"]['categories_ordered'],[topdict[x] for x in qdata.uns["dendrogram_['"+hierarchy_key+"']"]['categories_ordered']]))
sc.pl.matrixplot(qdata,groupby=hierarchy_key,var_names=var_dict,save='top_degenes',cmap='RdBu_r',use_raw=False,dendrogram=True)


# In[24]:


qdata.obs.leiden.value_counts()

