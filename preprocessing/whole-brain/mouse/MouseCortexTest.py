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
sc.settings.figdir='/wynton/group/ye/mtschmitz/figures/mouseTesting/'
scv.settings.figdir='/wynton/group/ye/mtschmitz/figures/mouseTesting/'
sc.settings.file_format_figs='pdf'
sc.settings.autosave=False
sc.settings.autoshow=True


# In[2]:


adata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityMouseWbAdultCortexProcessed.h5ad')


# In[3]:


adata=adata[~adata.obs.leiden.isin(['32','11','43','45','48','44','50','37']),:]


# In[4]:


adata=adata[:,adata.var['highly_variable']]


# In[5]:


adata.obs.drop('percent_mito',axis=1,inplace=True)
adata.obs.drop('percent_ribo',axis=1,inplace=True)
adata.obs.drop('latent_RT_efficiency',axis=1,inplace=True)
adata.obs.drop('latent_cell_probability',axis=1,inplace=True)


# In[6]:


adata.var=adata.var.loc[:,adata.var.columns.str.contains('id-|feature_type|name-')]


# In[5]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pl.umap(adata,color='region')


# In[8]:


sc.pl.umap(adata,color='dataset_name')


# In[9]:


sc.pl.umap(adata,color='timepoint')


# In[10]:


sc.pl.umap(adata,color='supervised_name',legend_loc='on data',legend_fontsize=5)


# In[12]:


sc.pl.umap(adata,color='leiden',legend_loc='on data')


# In[ ]:





# In[13]:


adata.obs['Adult']=adata.obs['timepoint']>31


# In[14]:


sc.pl.umap(adata,color=['LHX6','NR2F2','PROX1','GAD1','GAD2','SLC17A6','SLC17A7','MEIS2','SATB2','BCL11B','CUX1','CUX2','ETV1','RASGRF2','CALB1','CALB2','POU3F2','RORB','FEZF2','LRATD2','SLA2','FOXP2','NXPH4','CCN2','TLE4'],use_raw=False)


# In[7]:


allencells=pd.read_csv('/wynton/group/ye/mtschmitz/mousefastqpool/BICCN_metadata.csv')
allennamedf=allencells['sample_name'].str.split('-|_',expand=True)
adatasplitnamedf=pd.DataFrame(adata.obs.index)[0].str.split('_',expand=True)
allenids=pd.DataFrame(allennamedf[0]+'_'+allennamedf[1]+'_'+allennamedf[2]+'_'+allennamedf[4])
allenids.index=allenids[0]
allenids['realind']=allencells.index
adataids=pd.DataFrame(adatasplitnamedf[0]+'_'+adatasplitnamedf[2]+'_'+adatasplitnamedf[3]+'_'+adatasplitnamedf[5])
adataids.index=adataids[0]
adataids['realind']=adata.obs.index
adata.obs['allen_cluster_label']='nan'
adata.obs['allen_class_label']='nan'
adata.obs.loc[adataids.loc[set(adataids[0])&set(allenids[0]),'realind'],'allen_cluster_label']=allencells.loc[allenids.loc[set(adataids[0])&set(allenids[0]),'realind'],:]['cluster_label'].tolist()
adata.obs.loc[adataids.loc[set(adataids[0])&set(allenids[0]),'realind'],'allen_class_label']=allencells.loc[allenids.loc[set(adataids[0])&set(allenids[0]),'realind'],:]['class_label'].tolist()
adata.obs['simplified_allen']=adata.obs['allen_cluster_label'].str.split('_',expand=True)[1]
adata.obs['simplified_allen']=[re.sub('\sCTX','',str(x)) for x in adata.obs['simplified_allen']]


# In[10]:


adata.obs['simplified_allen'].value_counts()


# In[16]:


sc.pl.dotplot(adata[adata.obs.simplified_allen.str.contains('L[0-9]|CR'),:],groupby='simplified_allen',standard_scale='var',var_names=['SLC17A6','SLC17A7','MEIS2','SATB2','BCL11B','CUX1','CUX2','ETV1','RASGRF2','CALB1','CALB2','POU3F2','RORB','FEZF2','LRATD2','SLA2','FOXP2','FOXP4','NXPH4','CCN2','TLE4'])


# In[17]:


sc.pl.umap(adata,color=['simplified_allen'],legend_loc='on data',legend_fontsize=8)


# In[25]:


adata.obs['dev_excit']=False
adata.obs['dev_excit']=adata.obs.leiden.isin(['2','5','4','18','14','17','16','38','16','12','25'])


# In[20]:


sc.pl.umap(adata,color='dev_excit')


# In[26]:


adata.obs['condit']='Garbage'
adata.obs.loc[adata.obs['dev_excit'],'condit']='DevExcitatory'
adata.obs.loc[adata.obs.simplified_allen.str.contains('L[0-9]|CR'),'condit']='AdultExcitatory'


# In[27]:


sc.pl.umap(adata,color='condit')


# In[29]:


import diffxpy
import diffxpy.api as de


# In[73]:


adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)].todense()


# In[ ]:


sf=adata.X.mean(1) 
sf = sf/sf.mean()
adata.obs['size_factors']=sf

test = de.test.wald(
    data=adata[adata.obs['condit']!='Garbage',:],
    formula_loc="~ 1 + condit + size_factors",
    factor_loc_totest="condit",#as_numeric=['latent_time'],
    size_factors='size_factors',
    as_numeric=["size_factors"]
)
resultsdiffx=test.summary()
resultsdiffx['neg_log_q']=-test.log10_qval_clean()
resultsdiffx['signif']=(resultsdiffx['neg_log_q']>2) & (np.absolute(resultsdiffx['log2fc'])>1.2)
resultsdiffx=resultsdiffx.loc[resultsdiffx['log2fc'].argsort(),:]


# In[42]:


resultsdiffx['mean']


# In[57]:


with pd.option_context('display.max_rows', 30, 'display.max_columns', 10):
    print(resultsdiffx.loc[resultsdiffx['mean']>.3,:])
    
list(resultsdiffx.loc[resultsdiffx['mean']>.3,'gene'])


# In[ ]:


sc.pp.normalize_total(adata,exclude_highly_expressed=True)
sc.pp.log1p(adata)
#sc.pp.scale(adata,max_value=10)


# In[67]:


sc.pl.umap(adata,color=list(resultsdiffx.loc[resultsdiffx['mean']>.3,'gene'])[:20],use_raw=False)


# In[68]:


sc.pl.umap(adata,color=list(resultsdiffx.loc[resultsdiffx['mean']>.3,'gene'])[-20:],use_raw=False)


# In[65]:


resultsdiffx.to_csv('/wynton/home/ye/mschmitz1/MatureVsDev.csv')


# In[6]:


sc.pl.umap(adata,color=['MKI67','GDAP1L1','FOS','FOSB','JUN','JUNB','NPAS4','MEF2C','SYN1','DLG4','GRIN2A','GRIN2B'],use_raw=False)


# In[4]:


adata


# In[ ]:




