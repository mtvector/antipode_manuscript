#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import scanpy as sc
import numpy as np
import scipy
import sklearn
import matplotlib.pyplot as plt
import matplotlib
import sys
import loompy
import scipy.optimize
import os
import seaborn
import bbknn
import re
import importlib
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.expanduser("../../utils/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)


# In[2]:


adata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityPanHumanPresuperviseProcessed.h5ad')


# In[3]:


adata


# In[4]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pl.umap(adata,color='region')
sc.pl.umap(adata,color='general_region')


# In[5]:


adata=adata[np.random.choice(adata.obs.index,size=100000,replace=False),:]


# In[6]:


adata.write('/wynton/scratch/mtschmitz/Human100k.h5ad')


# In[2]:


genomepath='/wynton/group/ye/mtschmitz/refdata2/hg38/'
t2g=pd.read_csv(os.path.join(genomepath,'cDNA_introns_t2g.txt'),sep='\t',header=None)
t2g.loc[t2g[2].isna(),2]=t2g.loc[t2g[2].isna(),1]
t2gdict=dict(zip(t2g[1],t2g[2]))


# In[17]:


get_ipython().run_cell_magic('time', '', "filepath=os.path.expanduser(os.path.join('/wynton/scratch/mtschmitz/fastqpool/humanfastqpool'))\nfileList=os.listdir(filepath)\nfileList=[x for x in fileList if 'kout' in x.lower()]\nprint(fileList)\nfileList=[x for x in fileList if ('rhomb' in x.lower()) or ('mesence' in x.lower())]\n#fileList+=['GW20CGE_kOut']\nfileList\nprint(fileList)")


# In[39]:


adatas=[]
filename='aem_cellbended_150_750_175e_V0.2'
for f in fileList[1:]:
    print(f,flush=True)
    #try:
    if os.path.exists(os.path.join(filepath,f,filename,'aem_cellbended_filtered.h5')):
        sadata = sc_utils.readCellbenderH5(os.path.join(filepath,f,filename,'aem_cellbended_filtered.h5'))
        sadata =sadata[sadata.obs['latent_cell_probability']>.99,:]
        sc.pp.filter_cells(sadata,min_genes=500)
    else:
        #continue
        print(os.path.join(filepath,f,'all_em'),flush=True)
        sadata=sc_utils.loadPlainKallisto(os.path.join(filepath,f,'all_em'),min_genes=500)
    print(1,flush=True)
    sadata.obs.index=[re.sub("-1","",x) for x in sadata.obs.index]
    sadata.uns['name']=f
    sadata.obs['batch_name']=str(sadata.uns['name'])
    print(2,flush=True)
    print(2.01,flush=True)
    sadata.obs['timepoint']=sc_utils.tp_format_human(sadata.uns['name'])
    print(2.1,flush=True)
    regionstring=sc_utils.region_format_human(sadata.uns['name'])
    print(2.2,flush=True)
    regionstring=regionstring.lower()
    print(3,flush=True)
    sadata.obs['region']=regionstring
    adatas.append(sadata)
    #except:
    #    continue

#hdatas=hdatas+[multi]
hdata=sc.concat(adatas)


# In[20]:


adatas


# In[21]:


sc.pp.filter_genes(hdata,min_cells=30)
sc.pp.filter_cells(hdata,min_genes=400)


# In[22]:


hdata


# In[25]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pl.violin(hdata,groupby='batch_name',keys='n_genes')


# In[24]:


sc.pl.violin(hdata,groupby='batch_name',keys='n_genes')


# In[26]:


sc.pp.normalize_per_cell(hdata)
sc.pp.log1p(hdata)
sc.pp.highly_variable_genes(hdata,n_top_genes=8000)
sc.pp.scale(hdata,max_value=10)


# In[27]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pp.pca(hdata,n_comps=100)
#sc.pp.neighbors(bendedadata)
bbknn.bbknn(hdata,batch_key='batch_name')
sc.tl.leiden(hdata)
sc.tl.umap(hdata)
sc.pl.umap(hdata,color=['leiden','batch_name'])


# In[31]:


sc.pl.umap(hdata,color=['batch_name'],legend_fontsize=4)


# In[41]:


sc.pl.umap(hdata,color=['n_genes'],legend_fontsize=4)


# In[29]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pl.umap(hdata,color=['region'])


# In[30]:


sc.pl.umap(hdata,color=['batch_name'])


# In[42]:


easygenes= ['AIF1','MKI67','PCNA','AQP4','PDGFRA','DCX','RBFOX3','PAX6','TH','OTX2','TUBB3','TFAP2A','HOXA2','HOXA3','MNX1','PHOX2B','EN1','EN2','ISL1','GATA3','FOXA2','FOXJ1','TTR','GAD1','GAD2','SLC17A6','SLC17A7','SLC17A8','MET','TPH1','TPH2']
easygenes=[x for x in easygenes if x in hdata.var.index]
sc.pl.umap(hdata,color=easygenes,use_raw=False)


# In[40]:


sc.tl.rank_genes_groups(hdata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(hdata, n_genes=25, sharey=False)


# In[ ]:




