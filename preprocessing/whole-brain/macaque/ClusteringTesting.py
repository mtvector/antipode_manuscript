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


# In[2]:


newfile='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityDYNAMICALMouseWbGe.h5ad'
adata=sc.read(newfile)


# In[3]:


adata.obs.dtypes


# In[ ]:





# In[37]:


adata.var.dtypes


# In[ ]:


from random import random
adata=adata[np.random.choice(adata.obs.index,15000,replace=False),:]


# In[ ]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pl.umap(adata,color=['leiden','latent_time'])


# In[77]:


sc.pl.umap(adata,color=['leiden'],legend_loc='on data')


# In[18]:


sc.tl.leiden(adata,resolution=5)
sc.pl.umap(adata,color='leiden')


# In[21]:


sc.pp.regress_out(adata,keys='linear_velocity_pseudotime')


# In[22]:


adata=adata[:,np.isfinite(adata.X.sum(0))]


# In[23]:


sc.pp.pca(adata)


# In[4]:


adata


# In[43]:


sc.pl.umap(adata,color='leiden',legend_loc='on data')


# In[69]:


sufficient_cells=adata.obs['batch_name'].value_counts().index[adata.obs['batch_name'].value_counts()>10]
adata=adata[adata.obs['batch_name'].isin(sufficient_cells),:]

import bbknn
bbknn.bbknn(adata,batch_key='batch_name',use_annoy=False,metric='euclidean',neighbors_within_batch=9,pynndescent_n_neighbors=9)
sc.tl.leiden(adata,resolution=5)
sc.pl.umap(adata,color='leiden')


# In[70]:


sc.pl.umap(adata,color='supervised_name')


# In[71]:


sc.tl.umap(adata)


# In[ ]:


sc.pl.umap(adata,color='leiden')


# In[73]:


sc.pl.umap(adata,color='supervised_name')


# In[ ]:





# In[ ]:


sc.tl.leiden(adata,resolution=5)
sc.pl.umap(adata,color='leiden')


# In[7]:


sc.tl.dendrogram(adata,groupby='leiden',cor_method='pearson')
sc.pl.dendrogram(adata,groupby='leiden')


# In[7]:


adata=adata[:,adata.var.highly_variable]


# In[9]:


df=pd.DataFrame(adata[:,adata.var['highly_variable']].X)
df['leiden']=list(adata.obs.leiden)
df1=df.groupby(['leiden']).mean()
corr =df1.T.corr() 


# In[18]:


get_ipython().run_line_magic('matplotlib', 'inline')
seaborn.clustermap(df1.T.corr())


# In[11]:


adata


# In[19]:


sc.pl.umap(adata,color='leiden',legend_loc='on data')


# In[20]:


import statsmodels
import statsmodels.formula.api as smf
def rank_normalize(x):
    vptr=x.rank()
    vptr=(vptr-np.min(vptr))
    vptr=vptr/np.max(vptr)
    return(vptr)

    
for c in corr.columns:
    for p in corr[c].sort_values(ascending=False).index[1:4]:
        print(c,p)
        thisadata=adata[adata.obs.leiden.isin([c,p]),:].copy()
        thisadata.obs.loc[thisadata.obs['leiden']==c,'linear_velocity_pseudotime']=rank_normalize(thisadata.obs.loc[thisadata.obs['leiden']==c,'linear_velocity_pseudotime'])
        thisadata.obs.loc[thisadata.obs['leiden']==p,'linear_velocity_pseudotime']=rank_normalize(thisadata.obs.loc[thisadata.obs['leiden']==p,'linear_velocity_pseudotime'])
        parameters=[]
        for i in tqdm.tqdm(range(thisadata.shape[1])):
            df=pd.DataFrame({'y':thisadata.X[:,i],'linear_velocity_pseudotime':thisadata.obs['linear_velocity_pseudotime'],'leiden':thisadata.obs['leiden']})
            model = smf.ols('y ~ linear_velocity_pseudotime * leiden', data=df)
            model = model.fit()
            parameters.append(model.params)
        paramdf=pd.DataFrame(parameters)
        print(model.params)
        print(np.median(np.abs(paramdf.iloc[:,0])))
        print(np.median(np.abs(paramdf.iloc[:,1])))
        print(np.median(np.abs(paramdf.iloc[:,2])))
        print(np.median(np.abs(paramdf.iloc[:,3])))


# In[13]:


model.params


# In[55]:


df
# Initialise and fit linear regression model using `statsmodels`
model = smf.ols('y ~ latent_time + leiden', data=df)
model = model.fit()


# In[78]:


paramdf


# In[73]:


paramdf=pd.DataFrame(parameters)
seaborn.distplot(paramdf.iloc[:,0])
seaborn.distplot(paramdf.iloc[:,1])
seaborn.distplot(paramdf.iloc[:,2])


# In[74]:


np.mean(np.abs(paramdf.iloc[:,0]))


# In[76]:


np.mean(np.abs(paramdf.iloc[:,1]))


# In[75]:


np.mean(np.abs(paramdf.iloc[:,2]))


# In[33]:


import sklearn
print(sklearn.metrics.adjusted_rand_score(adata.obs['supervised_name'], adata.obs['leiden']))


# In[81]:


sc.tl.leiden(adata,resolution=3)
print(sklearn.metrics.adjusted_rand_score(adata.obs['supervised_name'], adata.obs['leiden']))
sc.tl.leiden(adata,resolution=5)
print(sklearn.metrics.adjusted_rand_score(adata.obs['supervised_name'], adata.obs['leiden']))
sc.tl.leiden(adata,resolution=7)
print(sklearn.metrics.adjusted_rand_score(adata.obs['supervised_name'], adata.obs['leiden']))
sc.tl.leiden(adata,resolution=10)
print(sklearn.metrics.adjusted_rand_score(adata.obs['supervised_name'], adata.obs['leiden']))


# In[86]:


adata.uns['linear_velocity_graph']


# In[87]:


scv.tl.paga(
    adata,vkey='linear_velocity',
    groups="leiden",
    root_key="initial_states_probs",
    end_key="terminal_states_probs",
    use_time_prior="latent_time",
)


# In[89]:


scv.pl.paga(adata,dashed_edges='connectivities_tree',basis='umap')


# In[88]:


adata.uns['paga']['connectivities_tree']


# In[94]:


most_common_sn=pd.DataFrame(adata.obs.groupby('leiden')['supervised_name'].value_counts().unstack().idxmax(1))
color_dict=dict(zip(adata.obs.supervised_name.cat.categories,adata.uns['supervised_name_colors']))
row_colors = [color_dict[x] for x in most_common_sn[0]]
mat=pd.DataFrame(adata.uns['paga']['connectivities'].todense(),index=adata.obs.leiden.cat.categories)
seaborn.clustermap(mat,method='average',row_colors=row_colors)


# In[ ]:


sc.tl.dendrogram(adata,groupby='leiden')
most_common_sn=pd.DataFrame(adata.obs.groupby('leiden')['agg_supervised_name'].value_counts().unstack().idxmax(1))
color_dict=dict(zip(adata.obs.agg_supervised_name.cat.categories,adata.uns['agg_supervised_name_colors']))
row_colors = [color_dict[x] for x in most_common_sn[0]]
mat=pd.DataFrame(adata.uns['dendrogram_leiden']['correlation_matrix'],index=adata.obs.leiden.cat.categories)
seaborn.clustermap(mat,method='average',row_colors=row_colors)

