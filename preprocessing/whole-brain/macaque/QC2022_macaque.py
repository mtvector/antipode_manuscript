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


get_ipython().run_cell_magic('time', '', "filepath=os.path.expanduser(os.path.join('/wynton/scratch/mtschmitz/fastqpool'))\nfileList=os.listdir(filepath)\nfileList=[x for x in fileList if 'kout' in x.lower()]\nprint(fileList)\nfileList=[x for x in fileList if ('encephal' in x.lower()) or ('encephal' in x.lower())]\n#fileList+=['GW20CGE_kOut']\nfileList\nfileList=fileList+['E65-2019A_Hypothalamus_kOut', 'E80-2019_Parietal_kOut']\nprint(fileList)")


# In[3]:


adatas=[]
filename='aem_cellbended_150_750_175e_V0.2'
for f in fileList:
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


# In[4]:


adatas


# In[5]:


sc.pp.filter_genes(hdata,min_cells=2)
sc.pp.filter_cells(hdata,min_genes=2)
sc.pp.filter_genes(hdata,min_cells=2)


# In[6]:


hdata.obs['batch_name'].value_counts()


# In[7]:


hdata.obs['n_counts']=hdata.X.sum(1)


# In[8]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pl.violin(hdata,groupby='batch_name',keys='n_counts',rotation=90)


# In[9]:


D1,D2,Mes2,Met1,Met2


# In[10]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pl.violin(hdata,groupby='batch_name',keys='n_genes',rotation=90)


# In[15]:


get_ipython().run_line_magic('matplotlib', 'inline')
#Old
sc.pl.violin(hdata,groupby='batch_name',keys='n_genes',rotation=90)


# In[11]:


sc.pp.normalize_per_cell(hdata)
sc.pp.log1p(hdata)
sc.pp.highly_variable_genes(hdata,n_top_genes=8000)
sc.pp.scale(hdata,max_value=10)


# In[12]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pp.pca(hdata,n_comps=50)
#sc.pp.neighbors(bendedadata)
bbknn.bbknn(hdata,batch_key='batch_name')
sc.tl.leiden(hdata)
sc.tl.umap(hdata)
sc.pl.umap(hdata,color=['leiden','batch_name'])


# In[19]:


sc.pl.umap(hdata[hdata.obs['batch_name'].str.contains('Myelencephalon_1'),:],color=['batch_name'],legend_fontsize=4)
sc.pl.umap(hdata[hdata.obs['batch_name'].str.contains('Myelencephalon_2'),:],color=['batch_name'],legend_fontsize=4)


# In[22]:


sc.pl.umap(hdata[hdata.obs['batch_name'].str.contains('Mesencephalon_1'),:],color=['batch_name'],legend_fontsize=4)
sc.pl.umap(hdata[hdata.obs['batch_name'].str.contains('Mesencephalon_2'),:],color=['batch_name'],legend_fontsize=4)


# In[23]:


sc.pl.umap(hdata[hdata.obs['batch_name'].str.contains('Diencephalon_1'),:],color=['batch_name'],legend_fontsize=4)
sc.pl.umap(hdata[hdata.obs['batch_name'].str.contains('Diencephalon_2'),:],color=['batch_name'],legend_fontsize=4)


# In[25]:


sc.pl.umap(hdata[hdata.obs['batch_name'].str.contains('Metencephalon_1'),:],color=['batch_name'],legend_fontsize=4)
sc.pl.umap(hdata[hdata.obs['batch_name'].str.contains('Metencephalon_2'),:],color=['batch_name'],legend_fontsize=4)


# In[21]:


sc.pl.umap(hdata,color=['batch_name'],legend_fontsize=4)


# In[20]:


sc.pl.umap(hdata,color=['n_genes'],legend_fontsize=4)


# In[15]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pl.umap(hdata,color=['region'])


# In[16]:


sc.pl.umap(hdata,color=['batch_name'])


# In[17]:


easygenes= ['AIF1','MKI67','PCNA','AQP4','PDGFRA','DCX','RBFOX3','PAX6','TH','OTX2','TUBB3','TFAP2A','HOXA2','HOXA3','MNX1','PHOX2B','EN1','EN2','ISL1','GATA3','FOXA2','FOXJ1','TTR','GAD1','GAD2','SLC17A6','SLC17A7','SLC17A8','MET','TPH1','TPH2']
easygenes=[x for x in easygenes if x in hdata.var.index]
sc.pl.umap(hdata,color=easygenes,use_raw=False)


# In[18]:


easygenes= ['NANOG','POU5F1','SOX2','KLF4','MKI67','MYCN','MYC','MYCL','BLM']
easygenes=[x for x in easygenes if x in hdata.var.index]
sc.pl.umap(hdata,color=easygenes,use_raw=False)


# In[23]:


hdata.obs['batch_name'].value_counts()


# In[13]:


destroy_cells=hdata[hdata.obs['leiden']=='1',:].obs.index


# In[14]:


hdata[hdata.obs['batch_name'].str.contains('Metencephalon_1'),:]


# In[15]:


#hdata=sc.concat(adatas)


# In[16]:


adatas[6]


# In[29]:


barcodes=pd.read_csv('/wynton/scratch/mtschmitz/fastqpool/E80-2022_Metencephalon_1_kOut/all_em/barcodes.tsv.gz',header=None)


# In[22]:


barcodes.loc[barcodes[0].isin(destroy_cells),0]=barcodes.loc[barcodes[0].isin(destroy_cells),0]+"_HUMAN"


# In[19]:


barcodes.to_csv('/wynton/scratch/mtschmitz/fastqpool/E80-2022_Metencephalon_1_kOut/all_em/human_barcodes.tsv',index=False,header=None)


# In[20]:


pd.DataFrame(destroy_cells).to_csv('/wynton/scratch/mtschmitz/fastqpool/E80-2022_Metencephalon_1_kOut/all_em/human_cells.csv',index=False,header=None)


# In[28]:


barcodes.loc[barcodes[0].isin(destroy_cells),:]


# In[18]:


sc.tl.rank_genes_groups(hdata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(hdata, n_genes=25, sharey=False)


# In[ ]:




