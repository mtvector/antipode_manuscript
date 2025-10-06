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
sc.settings.figdir='/wynton/group/ye/mtschmitz/figures/macWbSupervise/'
scv.settings.figdir='/wynton/group/ye/mtschmitz/figures/macWbSupervise/'
sc.settings.file_format_figs='pdf'
sc.settings.autosave=False
sc.settings.autoshow=True
import cellrank as cr


# In[3]:


newfile='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityDYNAMICALMouseWbGe.h5ad'
adata=sc.read(newfile)


# In[62]:


from random import random
adata=adata[np.random.choice(adata.obs.index,15000,replace=False),:]


# In[63]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pl.umap(adata,color='leiden')


# In[66]:


transition_matrix=scv.utils.get_transition_matrix(adata,vkey='linear_velocity',self_transitions=True)
transition_matrix=transition_matrix.astype(np.float64)
#transition_matrix=(transition_matrix.T/transition_matrix.sum(1)).T

transition_matrix[np.arange(0,transition_matrix.shape[0]),transition_matrix.argmax(1).A1]+=(1.0-transition_matrix.sum(1).A1)


vk = cr.tl.kernels.PrecomputedKernel(transition_matrix,adata=adata)
#vk = cr.tl.kernels.VelocityKernel(adata,vkey='linear_velocity')
ck = cr.tl.kernels.ConnectivityKernel(adata)
ck.compute_transition_matrix()
vk.compute_transition_matrix()

ckvk=vk+ck
ckvk.compute_transition_matrix()


# In[70]:


g = cr.tl.estimators.GPCCA(ckvk)
#g.compute_eigendecomposition()
g.compute_schur(n_components=50,method='krylov')
print('Schur complete',flush=True)
g.plot_spectrum(real_only=True,save='spectrum')


# In[71]:


g.compute_macrostates(n_states=50, n_cells=25,cluster_key="supervised_name")
print('Macrostates computed', flush=True)
#adata.write('/wynton/home/ye/mschmitz1/GPCCAmacrostatesmouseAdult2.h5ad')
adata.uns['macrostates']=list(g.macrostates.cat.categories)
g.compute_terminal_states()
cr.tl.initial_states(adata, cluster_key='supervised_name')


# In[74]:


adata.obs['macrostates']=g.macrostates


# In[75]:


sc.pl.umap(adata,color='macrostates')


# In[13]:


get_ipython().run_line_magic('matplotlib', 'inline')
scanpy.set_figure_params(figsize=(40,40))
scv.pl.velocity_graph(adata,vkey='linear_velocity',alpha=.1)


# In[24]:


sc.tl.leiden(adata,resolution=5,directed=False)
sc.pl.umap(adata,color='leiden',size=20)


# In[27]:


sc.pl.umap(adata,edges=True,size=20,alpha=.0005)


# In[20]:


adata.uns['linear_velocity_graph']


# In[47]:


scv.tl.paga(adata,groups='leiden',vkey='linear_velocity',minimum_spanning_tree=False,use_time_prior='linear_velocity_pseudotime')


# In[61]:


seaborn.clustermap(adata.uns['paga']['transitions_confidence'].todense(),xticklabels=True,yticklabels=True,row_linkage=scipy.cluster.hierarchy.linkage(adata.uns['paga']['transitions_confidence'].todense()),col_linkage=scipy.cluster.hierarchy.linkage(adata.uns['paga']['transitions_confidence'].todense().T))


# In[59]:


scv.pl.paga(adata,basis='umap',dashed_edges='transitions_confidence')


# In[51]:


sc.pl.umap(adata,color='leiden',legend_loc='on data',size=50,legend_fontsize=45)


# In[54]:


scipy.cluster.hierarchy.dendrogram(scipy.cluster.hierarchy.linkage(adata.uns['paga']['transitions_confidence'].todense()))


# In[23]:


sc.tl.leiden(adata,resolution=5,adjacency=adata.uns['linear_velocity_graph'],key_added='velo_leiden')
sc.pl.umap(adata,color='velo_leiden',size=20)


# In[ ]:





# In[31]:


df=pd.DataFrame(adata[:,adata.var['highly_variable']].X)
df['leiden']=list(adata.obs.leiden)
df1=df.groupby(['leiden']).mean()
corr =df1.T.corr() 


# In[30]:


seaborn.clustermap(df1.T.corr())


# In[63]:


adata=adata[:,adata.var.highly_variable]


# In[79]:


import statsmodels
import statsmodels.formula.api as smf

for c in corr.columns:
    for p in corr[c].sort_values(ascending=False).index[1:4]:
        print(c,p)
        thisadata=adata[adata.obs.leiden.isin([c,p]),:].copy()
        parameters=[]
        for i in tqdm.tqdm(range(thisadata.shape[1])):
            df=pd.DataFrame({'y':thisadata.X[:,i],'latent_time':thisadata.obs['latent_time'],'leiden':thisadata.obs['leiden']})
            model = smf.ols('y ~ latent_time + leiden', data=df)
            model = model.fit()
            parameters.append(model.params)
        paramdf=pd.DataFrame(parameters)
        print(np.mean(np.abs(paramdf.iloc[:,0])))
        print(np.mean(np.abs(paramdf.iloc[:,1])))
        print(np.mean(np.abs(paramdf.iloc[:,2])))


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

