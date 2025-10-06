#!/usr/bin/env python
# coding: utf-8

# # Integrated data analysis with neuronal subclustering v3#
# This version includes neuonal subsluctering, batch distribution, but no differential QC cut-offs

# In[38]:


import scanpy as sc
import pandas as pd
import anndata
import os
import re
import numpy as np
import scipy
import seaborn
import bbknn
import matplotlib
import matplotlib.pyplot as plt
import importlib
import importlib.util
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.expanduser("../../utils/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)
figdir='/wynton/group/ye/mtschmitz/figures/macaqueWbHypoPresupervise/'
sc.settings.figdir=figdir
sc.settings.file_format_figs='pdf'
sc.settings.autosave=True
sc.settings.autoshow=True


# In[2]:


pwd


# In[3]:


#set a path to your working directory
datapath='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_kallisto'


# In[4]:


#some helper functions that Matthew has written for processing irregularities in the macaque dataset
def tp_format_macaque(name):
    if isinstance(name, str):
        name=re.sub('PEC_YALE','E110',name)
        name=re.sub('PEC_Yale','E110',name)
        name=re.sub('Mac2','E65',name)
        searched=re.search("E[0-9]+",name)
        if searched is None:
            return('nan')
        tp=re.sub("^E","",searched.group(0))
        tp=float(tp)
        return(tp)
    else:
        return('nan')

def region_format_macaque(name):
    if isinstance(name, str):
        name=name.upper()
        name=re.sub('PEC_YALE_SINGLECELLRNASEQ_RMB[0-9]+','',name)
        name=re.sub('MAC2','',name)
        name=re.sub("E80MGE","E80_GE",name)
        name=re.sub('E[0-9]+',"",name)
        name=re.sub("-2019A_","",name)
        name=re.sub("-2019B_","",name)
        name=re.sub("-2019A_AND_E65-2019B_","",name)
        name=re.sub("-2019_","",name)
        name=re.sub("^_","",name)
        name=re.sub("CGE-LGE","CGE_AND_LGE",name)
        name=re.sub("LGE_AND_CGE","CGE_AND_LGE",name)
        name=re.sub("_1","",name)
        name=re.sub("_2","",name)
        name=re.sub("_OUT","",name)
        name=re.sub("_KOUT","",name)
        name=re.sub(".LOOM","",name)
        return(name)
    else:
        region='nan'
    return(region)

def macaque_correct_regions(regionlist):
    regionkey={'thal':'thalamus',
               'hippo':'hippocampus',
               'somato':'somatosensory',
               'hypo':'hypothalamus',
               'midbrain/pons':'midbrain',
               'motor-1':'motor',
               'motor-2':'motor',
               'temp':'temporal',
               'hip':'hippocampus',
               'cbc':'cerebellum',
               'md':'thalamus'}
    newl=[]
    for l in regionlist:
        l=l.lower()
        if l in regionkey.keys():
            l=regionkey[l]
        newl.append(l)
    return(newl)


# In[11]:


fileNameList=os.listdir(datapath)
#filter the directories in your directory
fileList=[f for f in fileNameList if os.path.exists(os.path.join(datapath,f,'outs','filtered_feature_bc_matrix.h5'))]
fileList=pd.DataFrame(fileList)[0]
fileList=fileList[fileList.str.contains('ypo|Mix7|thal|Mix2|Mix6')]


# In[12]:


print(fileNameList)


# In[13]:


print(fileList)


# In[64]:


adatas=[]

for f in fileList:
    print(f)
    adata = sc.read_10x_h5(os.path.join(datapath,f,'outs','filtered_feature_bc_matrix.h5'))
    #Get rid of -1 in cell barcodes
    adata.obs.index=[re.sub("-[0-9]+","",x) for x in adata.obs.index]
    adata.uns['name']=f
    print (adata.uns)
    adata.obs['batch_name']=re.sub('_Out','',str(adata.uns['name']))
    print (adata.obs)
    adata.obs['clean_cellname']=[re.sub('-[0-9]+','',x) for x in  adata.obs.index]
    adata.obs['full_cellname']=adata.obs['clean_cellname'].astype(str)+'_'+adata.obs['batch_name'].astype(str)
    adata.obs.index=adata.obs['full_cellname']
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    regionstring=sc_utils.region_format_macaque(sc_utils.macaque_process_irregular_names(adata.uns['name']))
    regionstring=regionstring.lower()
    adata.obs['region']=regionstring
    adatas.append(adata)
  


# In[65]:


print(adatas)


# In[66]:


print(adata.obs.index)


# In[67]:


adata=anndata.AnnData.concatenate(*adatas) #This line needs to be run twice to work - I don't know why!


# In[68]:


#saving all the datasets together will allow you to pick up where you left off
#h5ad stores whole anndata data structure
#adata.write('/Users/sara/Documents/Jupyter notebooks/macaque2018integration/Macaque_Hypothalamus_integrated_v2.h5ad')        


# In[69]:


adata.raw=adata


# In[70]:


get_ipython().run_line_magic('matplotlib', 'inline')
ribo_genes=[name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL') ]
adata.obs['percent_ribo'] = np.sum(
adata[:, ribo_genes].X, axis=1) / np.sum(adata.X, axis=1)
mito_genes = [name for name in adata.var_names if name in ['ND1','ND2','ND4L','ND4','ND5','ND6','ATP6','ATP8','CYTB','COX1','COX2','COX3'] or 'MT-' in name]
adata.obs['percent_mito'] = np.sum(
adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
print(adata)
sc.pl.violin(adata,groupby='batch_name',keys=['percent_ribo','percent_mito'],rotation=45,save='ribomito')


# In[71]:


adata.obs.groupby(['batch_name'])['percent_ribo'].median()


# In[72]:


mitochondria_cutoffs = {
    'E40_Hypothalamus': 0.15,
    'E50_Hypothalamus': 0.15,
    'E65A_Hypothalamus': 0.15,
    'E65B_Hypothalamus': 0.15,
    'E90_Hypothalamus': 0.3,
    'E100_Hypothalamus+POA': 0.15,
}


# In[73]:


print(mitochondria_cutoffs['E40_Hypothalamus'])


# In[74]:


'''cells=[]
for b in adata.obs['batch_name'].unique():
    cells+=list(adata.obs.index[(adata.obs['percent_mito']<mitochondria_cutoffs[b]) & adata.obs['batch_name']==b])
print(cells)
adata=adata._inplace_subset_obs(cells)'''


# In[75]:


#set these cutoffs depending on the dataset
adata._inplace_subset_obs(adata.obs['percent_ribo']<.4)
adata._inplace_subset_obs(adata.obs['percent_mito']<.2) #I used either 0.15 or 0.3 for seperate batch analsyis
#Get rid of mito genes
#adata._inplace_subset_var(~adata.var.index.isin(mito_genes))
print(adata)


# In[76]:


#Plot the distribution of number of genes per cell
seaborn.distplot((adata.X>0).sum(0).A1)
plt.title('Genes Per Cell')


# In[77]:


#basic filtering to get n_genes added to obs
sc.pp.filter_cells(adata, min_genes=800)
sc.pp.filter_genes(adata, min_cells=5)


# In[78]:


# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1)


# In[79]:


#n_genes after basic filtering
sc.pl.violin(adata,groupby='batch_name',keys=['n_counts','n_genes'],rotation=45,save='n_genes') #Not working if not doing the basic filtering to get obs: n_genes 


# In[80]:


print(adata)


# In[81]:


sc.pp.filter_genes(adata,min_cells=10) #I used 3 for seperate batch analsyis
#sc.pp.filter_cells(adata,min_counts=1000)
sc.pp.filter_cells(adata,min_genes=800) #I used 600 or 1000 for seperate batch analsyis
#sc.pp.filter_cells(adata,max_genes=4000) #I introduced this line
#sc.pl.highest_expr_genes(adata, n_top=20, )
print(adata)
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=5000,batch_key='batch_name',subset=False)
sc.pp.scale(adata,max_value=10)
sc.pp.pca(adata,n_comps=50) #I used 50 for seperate batch analsyis
#sc.pp.neighbors(adata)
bbknn.bbknn(adata,batch_key='batch_name',n_pcs=50,neighbors_within_batch=3)
sc.tl.leiden(adata,resolution=1.5)
sc.tl.umap(adata,spread=2)


# In[82]:


sc.pl.umap(adata,color=['leiden'],save='leiden')
sc.pl.umap(adata,color=['leiden'],save='leiden',legend_loc='on data')
sc.pl.umap(adata,color=['batch_name'],save='batch_name')
#sc.pl.umap(adata,color=['timepoint'],save='timepoint')


# In[83]:


sc.pl.umap(adata,color=['region'],save='region')


# # Supervised annotation of clusters
# It's often a large waste of time to cluster and recluster. The section below provides tools for analyzing the clusters and how to save and reload your annotations

# In[62]:


#Plot some marker genes so you get a sense of what is where
easygenes= ['MKI67','SOX2','DCX','VIM','FABP7','NFIA','SOX9','GFAP','AQP4','S100B','AGT','GLUL','MFGE8','SLC6A11','PTCH1','BCAS1','OLIG2','MBP','GAP43','SYT1','SLC17A6','SLC32A1','AIF1','C1QB','COL4A1','FN1','PDGFRB','FOXJ1','CLDN5','NKX2-1','FOXG1','BARHL1','PITX2','NHLH1','DLX1','LHX1','LHX6','MEIS2','ISL1','GAD2','ONECUT2','ZIC1','OTP','FEZF2','POMC','SST','AVP','TH','ADCYAP1']
easygenes=[x for x in easygenes if x in adata.var.index]
sc.pl.umap(adata,color=easygenes,use_raw=False,save='FaveGenes')


# In[92]:


sc.pl.umap(adata,color=pd.DataFrame(adata.uns['rank_genes_groups']['names'])['9'][0:10],use_raw=False)


# In[93]:


sc.pl.umap(adata,color=pd.DataFrame(adata.uns['rank_genes_groups']['names'])['0'][0:10],use_raw=False)


# In[91]:


sc.pl.umap(adata,color=pd.DataFrame(adata.uns['rank_genes_groups']['names'])['1'][0:10],use_raw=False)


# In[94]:


for i in pd.DataFrame(adata.uns['rank_genes_groups']['names']).columns:
    print('CLUSTER:')
    print(i)
    sc.pl.umap(adata,color=pd.DataFrame(adata.uns['rank_genes_groups']['names'])[i][0:10],use_raw=False)


# In[ ]:


#Change the values of this dictionary to assign the name you want to each of the leiden clusters
labeldict={'0':'Astrocytes (MFGE8+/GLUL+/SLC6A11+)','1':'Neuronal cluster 1','2':'OPCs','3':'Neuronal cluster 2',
           '4':'Ependymal cells','5':'Neuronal cluster 3','6':'Cycling cells (VIM+ or OLIG2+)','7':'Radial glia',
 '8':'Astrocytes (AGT+/GFAP+/S100B+)','9':'Early progenitor cells','10':'Vascular cells','11':'Neuronal cluster 4',
           '12':'Oligodendrocytes','13':'Microglia','14':'Endothelial cells', '15':'Neuronal cluster 5', '16':'Glial progenitors'}


# In[ ]:


#Change the values of this dictionary to assign the name you want to each of the leiden clusters
labeldict={'0':'PTPRZ1/EDNRB','1':'AST_MGFE8/GLUL','2':'','3':'',
           '4':'AST_MLF1/VIM','5':'','6':'','7':'NEU',
 '8':'NEU','9':'AST_ANKRD65/S100B/EFCAB10','10':'NEU','11':'',
           '12':'','13':'MG_AIF1/CCL3','14':'NEU', '15':'', '16':'','17':'', '18':'EPENDY_FOXJ1'}


# In[84]:


#basic logistic regression to get list of differential expression
#sc.tl.filter_rank_genes_groups(adata,groupby='leiden',min_in_group_fraction=.2,max_out_group_fraction=.1)
sc.tl.rank_genes_groups(adata,groupby='leiden',min_in_group_fraction=.2,max_out_group_fraction=.1)

sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)


# In[85]:


# Show the 50 top ranked genes per cluster in a dataframe.
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)


# In[86]:


result=adata.uns['rank_genes_groups']
groups=result['names'].dtype.names
df=pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores']})
with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print(df.iloc[0:50,:])


# In[87]:


print(adata)


# In[ ]:


#supervising cluster naming, which cluster is which?
sc.tl.paga(adata,groups='leiden')
sc.pl.paga_compare(adata,legend_fontsize=6,arrowsize=10,edge_width_scale=.4,threshold=np.quantile(adata.uns['paga']['connectivities'].data,.9))


# In[ ]:


#Assign a supervised name for each leiden cluster
#Will error out if you missed a cluster
adata.obs['supervised_name']=[labeldict[x] for x in adata.obs['leiden']]
sc.pl.umap(adata,color=['supervised_name'],save='supervised_name')


# In[ ]:


#What genes separates the neuronal cluster (cluster 1, 3, 5, 11, 15)?
sc.tl.rank_genes_groups(adata, 'leiden', groups=['1'], reference='3', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['1'], n_genes=20)


# In[ ]:


sc.tl.rank_genes_groups(adata, 'leiden', groups=['3'], reference='1', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['3'], n_genes=20)


# In[ ]:


sc.tl.rank_genes_groups(adata, 'leiden', groups=['5'], reference='11', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['5'], n_genes=20)


# In[ ]:


sc.tl.rank_genes_groups(adata, 'leiden', groups=['15'], reference='1', method='wilcoxon')
sc.pl.rank_genes_groups(adata, groups=['15'], n_genes=20)


# In[ ]:


sc.pl.violin(adata, 'n_genes', groupby='leiden')


# In[ ]:


# After annotating the cell types, visualize the marker genes
ax = sc.pl.dotplot(adata, easygenes, color_map='Blues', groupby='supervised_name', use_raw=False)


# In[ ]:


# Violin plot
ax = sc.pl.stacked_violin(adata, easygenes, groupby='supervised_name', use_raw=False, rotation=90)


# In[ ]:


#Save your supervised names with the ATGC.._batch-name cell names 
#Note, if you change the cell names or change QC, you may get nan instead
adata.obs.loc[:,['full_cellname','supervised_name']].to_csv('/Users/sara/Documents/Jupyter notebooks/macaque2018integration/MacaqueHypoSupervisedCellNames.txt',index=None)


# # Visualizing distributions across batches#

# In[ ]:


sc.tl.embedding_density(adata, groupby='batch')
sc.pl.embedding_density(adata, groupby='batch')


# In[ ]:


sc.tl.embedding_density(adata, groupby='batch_name')
sc.pl.embedding_density(adata, groupby='batch_name')


# In[ ]:


for batch_name in ['E40_Hypothalamus', 'E50_Hypothalamus', 'E65A_Hypothalamus', 'E65B_Hypothalamus', 'E90_Hypothalamus', 'E100_Hypothalamus+POA']:
    sc.pl.umap(adata, color='batch_name', groups=[batch_name])


# In[ ]:


f = plt.figure()
df_plot = adata.obs.groupby(['supervised_name', 'batch_name']).size().reset_index().pivot(columns='batch_name', index='supervised_name', values=0).apply(lambda g: g / g.sum(),1)
ax = df_plot.plot(kind='bar', legend=False,stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('proportion')
ax.get_figure().savefig(os.path.join(sc.settings.figdir,'supervisedVbatchBar.pdf'), bbox_inches="tight")


# In[ ]:


f = plt.figure()
df_plot = adata.obs.groupby(['supervised_name', 'batch_name']).size().reset_index().pivot(columns='batch_name', index='supervised_name', values=0)
ax = df_plot.plot(kind='bar', legend=False,stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('number of cells')
ax.get_figure().savefig(os.path.join(sc.settings.figdir,'supervisedVbatchBar.pdf'), bbox_inches="tight")


# In[ ]:


print(adata)


# In[ ]:


f = plt.figure()
df_plot = adata.obs.groupby(['batch_name']).size().reset_index().pivot(columns='batch_name', values=0)
ax = df_plot.plot(kind='bar', stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('number of cells')


# # Subclutering neuronal clusters#

# In[ ]:


adata_subset = adata[adata.obs['leiden'].isin(['1', '3', '5', '11', '15'])]


# In[ ]:


print(adata_subset)


# In[ ]:


sc.pp.pca(adata_subset,n_comps=100)


# In[ ]:


bbknn.bbknn(adata_subset,batch_key='batch_name',n_pcs=100,neighbors_within_batch=3)


# In[ ]:


sc.tl.leiden(adata_subset,resolution=1.5)


# In[ ]:


sc.tl.umap(adata_subset,spread=2)


# In[ ]:


sc.pl.umap(adata_subset,color=['leiden'])
sc.pl.umap(adata_subset,color=['batch_name'])


# In[ ]:


subsetgenes= ['MAP2','GAP43','SYT1','SOX2','DCX','SLC17A6','SLC32A1','FOXG1','NKX2-1','BARHL1','PITX2','LMX1A','LMX1B','ARX','SIX3','NHLH1','DLX1','DLX2','DLX5','DLX6','LHX1','LHX6','LHX9','MEIS2','PAX6','SIM1','SP9','ISL1','GAD2','ONECUT1','ONECUT2','ZIC1','OTP','FEZF2','NPY','POMC','SST','AVP','OTX','KISS1','TH','ADCYAP1' ]
subsetgenes=[x for x in subsetgenes if x in adata_subset.var.index]


# In[ ]:


sc.pl.umap(adata_subset,color=subsetgenes, use_raw=False)


# In[ ]:


ax = sc.pl.dotplot(adata_subset, subsetgenes, color_map='Blues', groupby='leiden', use_raw=False)


# In[ ]:


ax = sc.pl.stacked_violin(adata_subset, subsetgenes, groupby='leiden', use_raw=False, rotation=90)


# In[ ]:


sc.pl.violin(adata_subset, 'n_genes', groupby='leiden')


# In[ ]:


sc.tl.rank_genes_groups(adata_subset,groupby='leiden',min_in_group_fraction=.2,max_out_group_fraction=.1)
sc.pl.rank_genes_groups(adata_subset, n_genes=25, sharey=False)


# In[ ]:


# Show the 50 top ranked genes per cluster in a dataframe.
pd.DataFrame(adata_subset.uns['rank_genes_groups']['names']).head(50)


# In[ ]:


print(adata_subset)


# In[ ]:


sc.pl.embedding_density(adata_subset, groupby='batch')
sc.pl.embedding_density(adata_subset, groupby='batch_name')


# In[ ]:


for batch_name in ['E40_Hypothalamus', 'E50_Hypothalamus', 'E65A_Hypothalamus', 'E65B_Hypothalamus', 'E90_Hypothalamus', 'E100_Hypothalamus+POA']:
    sc.pl.umap(adata_subset, color='batch_name', groups=[batch_name])


# In[ ]:


f = plt.figure()
df_plot = adata_subset.obs.groupby(['leiden', 'batch_name']).size().reset_index().pivot(columns='batch_name', index='leiden', values=0).apply(lambda g: g / g.sum(),1)
ax = df_plot.plot(kind='bar', legend=False,stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('proportion')
ax.get_figure().savefig(os.path.join(sc.settings.figdir,'supervisedVbatchBar2.pdf'), bbox_inches="tight")


# In[ ]:


f = plt.figure()
df_plot = adata_subset.obs.groupby(['leiden', 'batch_name']).size().reset_index().pivot(columns='batch_name', index='leiden', values=0)
ax = df_plot.plot(kind='bar', legend=False,stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('number of cells')
ax.get_figure().savefig(os.path.join(sc.settings.figdir,'supervisedVbatchBar2.pdf'), bbox_inches="tight")


# In[ ]:


f = plt.figure()
df_plot = adata_subset.obs.groupby(['batch_name']).size().reset_index().pivot(columns='batch_name', values=0)
ax = df_plot.plot(kind='bar', stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('number of cells')


# In[ ]:




