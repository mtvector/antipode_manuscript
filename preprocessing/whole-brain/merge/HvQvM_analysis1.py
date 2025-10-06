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
import tqdm
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
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.expanduser("../../utils/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)
spec.loader.exec_module(sc_utils)
sc.settings.figdir='/wynton/group/ye/mtschmitz/figures/hvqvmAnalysis1/'
scv.settings.figdir='/wynton/group/ye/mtschmitz/figures/hvqvmAnalysis1/'
sc.settings.file_format_figs='pdf'
sc.settings.autosave=False
sc.settings.autoshow=True
import matplotlib.font_manager
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Nimbus Sans','Arial']})
matplotlib.rc('text', usetex=False)
scanpy.set_figure_params(scanpy=True,dpi_save=300)
pd.set_option('display.max_rows', 500)


# In[2]:


mdata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT_chang_h5ad/KDCbVelocityMouseWbPresupervision.h5ad')
mdata.obs['species']='mouse'
mdata.var.index=[x.upper() for x in mdata.var.index]
print(mdata.var.index)
mdata=sc_utils.sum_duplicate_var(mdata)
#mdata.X=mdata.raw.X[:,mdata.raw.var.index.isin(mdata.var.index)]#.todense()
#mdata.X=scipy.sparse.csr_matrix(mdata.X)

orthos=pd.read_csv('/wynton/home/ye/mschmitz1/utils/HOM_AllOrganism.rpt',sep='\t')
orthos=orthos.loc[orthos['NCBI Taxon ID'].isin([10090,9606]),:]
classcounts=orthos['DB Class Key'].value_counts()
one2one=classcounts.index[list(classcounts==2)]
orthos=orthos.loc[orthos['DB Class Key'].isin(one2one),:]

htab=orthos.loc[orthos['NCBI Taxon ID']==9606,:]
mtab=orthos.loc[orthos['NCBI Taxon ID']==10090,:]
genemapping=dict(zip([x.upper() for x in mtab['Symbol']],htab['Symbol']))
mdata=mdata[:,mdata.var.index.isin(genemapping.keys())]
mdata.var.index=[genemapping[x] for x in mdata.var.index]
#mdata=mdata[np.random.choice(mdata.obs.index,300000,replace=False),:]

qdata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT_chang_h5ad/KDCbVelocityMacaqueWbPresupervise.h5ad')
qdata.obs['species']='macaque'
#qdata.X=qdata.raw.X[:,qdata.raw.var.index.isin(qdata.var.index)]#.todense()
#qdata.X=scipy.sparse.csr_matrix(qdata.X)
#qdata=qdata[np.random.choice(qdata.obs.index,300000,replace=False),:]
qdata=sc_utils.sum_duplicate_var(qdata)

hdata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT_chang_h5ad/KDCbVelocityPanHumanPresupervise.h5ad')
hdata.obs['species']='human'
#hdata.X=hdata.raw.X[:,hdata.raw.var.index.isin(hdata.var.index)]#.todense()
#hdata.X=scipy.sparse.csr_matrix(hdata.X)
#hdata=hdata[np.random.choice(hdata.obs.index,300000,replace=False),:]
hdata=sc_utils.sum_duplicate_var(hdata)

qdata.var_names_make_unique()
mdata.var_names_make_unique()
hdata.var_names_make_unique()
qdata.obs_names_make_unique()
mdata.obs_names_make_unique()
hdata.obs_names_make_unique()

adata=anndata.AnnData.concatenate(mdata,qdata,hdata,batch_key='species_batch')


# In[13]:


mdata.var.index=[x.upper() for x in mdata.var.index]


# In[14]:


mdata.var.index


# In[8]:


[genemapping[x] for x in mdata.var.index]


# In[22]:


adata=anndata.AnnData.concatenate(mdata,qdata,hdata,batch_key='species_batch')


# In[3]:


mdata


# In[ ]:





# In[2]:


adata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/HvQvMsmall.h5ad')


# In[15]:


adata.obs['species']


# In[180]:


gw20fix={'GW20V1':'GW20parietal',
 'GW20somato':'GW20PFC',
 'GW20cingulate-corpus':'GW20hypoventalmedial',
 'GW20claustrum':'GW20MGE',
 'GW20preoptic':'GW20caudalthalamus',
 'GW20motor':'GW20CGE',
 'GW20LGE':'GW20dorsalthalamus',
 'GW20putanum':'GW20ventralthalamus'}


# In[181]:


gw20fixnew={}
for x in gw20fix.keys():
    gw20fixnew[x+'_kOut']=gw20fix[x]+'_kOut'
    gw20fixnew[gw20fix[x]+'_kOut']=x+'_kOut'


# In[182]:


adata.obs['batch_name']=adata.obs['batch_name'].replace(gw20fixnew)


# In[183]:


list(adata.obs.loc[adata.obs.batch_name.str.contains('GW20'),'batch_name'].unique())


# In[184]:


list(adata.obs['batch_name'].unique())


# In[215]:


reg_dict={'cerebel|cbc|ce|vermis|cb':'Cb','vz|cortex|ctx|wt[0-9]+|Neurons_Sample_':'Ctx','ob|GSM3449':'OB','head|all|nan':'head','subpall|sp':'FB','thal|pulv|lgn':'Thal',
 'stria|stiatum|putanum|putamen|caud|accumb|nac|basal|SAMN107673':'Str','hypo':'Hypo','gp':'GP','clau':'Ctx-Clau','amy':'Amy',
 'ge':'GE','lge':'LGE','mge':'MGE','somato|som':'Ctx-Somato','cge':'CGE','motor|mop|m1|dfc':'Ctx-MOp',
 'pfc|prefront':'Ctx-PFC','cing':'Ctx-Cing','v1|occ':'Ctx-Occ','central|par':'Ctx-Par','temp':'Ctx-Temp',
 'hippo|ca[0-4]|dent|hc':'Ctx-Hip','forebrain|prosen|diencepha|telenceph':'FB','midbrain|md|mb':'MB','hind':'HB','medull':'Medulla','choroid|cp':'Choroid',
 'sept':'Sept','pre-opt|poa|pos|preoptic':'POA','pons':'Pons','drg':'DRG','org|cultured|rna':'cultured','insula':'Ctx-Ins'
}


# In[216]:


adata.obs['batch_name']=list(adata.obs['batch_name'].astype(str))
adata.obs['region']=list(adata.obs['region'].astype(str))
adata.obs['msregion']=list(adata.obs['msregion'].astype(str))

msd=adata.obs['batch_name'].str.contains('multi',case=False)
adata.obs.loc[msd,'msregion']=adata.obs.loc[msd,'region']
adata.obs['msregion']=adata.obs['msregion'].astype(str)
for k in reg_dict.keys():
    print(k)
    adata.obs.loc[adata.obs['msregion'].str.contains(k,case=False),'region']=reg_dict[k]


# In[217]:



for k in reg_dict.keys():
    print(k)
    print(adata.obs.loc[adata.obs['batch_name'].str.contains(k,case=False),'batch_name'])
    adata.obs.loc[adata.obs['batch_name'].str.contains(k,case=False),'region']=reg_dict[k]
    


# In[218]:


sorted(list(adata.obs['region'].unique()))


# In[219]:


adata[adata.obs['region']=='cortex',:].obs


# In[279]:


region_groupings={'ctx': ['Ctx-PFC','Ctx-Cing','Ctx-Clau','Ctx-Ins','Ctx-Hip','Ctx-MOp', 'Ctx-Somato','Ctx','Ctx-Temp','Ctx-Par','Ctx-Occ'],
 'bn':[ 'POA','GP','Sept','Str','Amy'],
 'ob':['OB'],
 'cp':['Choroid'],
 'mb':['MB'],
 'hb':['Pons','Cb','HB','DRG'],
 'ge':['MGE','GE','LGE','CGE'],
 'de':['Hypo', 'FB','Thal'],
 'h':['head'],
 'xvctx':['cultured']}


# In[280]:


palettes={
'ctx':seaborn.color_palette("YlOrBr",n_colors=len(region_groupings['ctx'])+2).as_hex()[2:],
'bn':seaborn.color_palette("Blues",n_colors=len(region_groupings['bn'])+2).as_hex()[2:],
'ge':seaborn.color_palette("Reds",n_colors=len(region_groupings['ge'])+2).as_hex()[2:],
'hb':seaborn.color_palette("Greens",n_colors=len(region_groupings['hb'])+2).as_hex()[2:],
'de':seaborn.color_palette("Purples",n_colors=len(region_groupings['de'])+2).as_hex()[2:],
'ob':['magenta'],
'mb':['cyan'],#turquoise
'cp':['slategray'],
'xvctx':['lightgray'],
'h':['tan']
}


# In[281]:


adata.obs['region']=adata.obs['region'].astype('category')
regions=[]
for k in palettes.keys():
    regions+=region_groupings[k]
    print(region_groupings[k])
region_colors=[]
for k in palettes.keys():
    region_colors+=palettes[k]
    print(palettes[k])
paldict=dict(zip(regions,region_colors))

region_colors=[paldict[x] for x in regions if x in adata.obs['region'].cat.categories]
adata.obs['region'].cat.categories=[x for x in regions if x in adata.obs['region'].cat.categories]


# In[282]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pl.umap(adata,color='region',palette=region_colors)


# In[283]:


sc.pl.umap(adata,color='species')


# In[288]:


adata.obs['batch_name']


# In[289]:


adata.obs['batch_name']


# In[77]:


def corr2_coeff(A, B):
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]
    # Sum of squares across rows
    ssA = (A_mA**2).sum(1)
    ssB = (B_mB**2).sum(1)
    # Finally get corr coeff
    return np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:, None],ssB[None]))

def most_frequent(List): 
    return max(set(List), key = List.count)

def variance(X,axis=1):
    return( (X.power(2)).mean(axis=axis)-(np.power(X.mean(axis=axis),2)) )

orthos=pd.read_csv('/wynton/home/ye/mschmitz1/utils/HOM_AllOrganism.rpt',sep='\t')
orthos=orthos.loc[orthos['NCBI Taxon ID'].isin([10090,9606]),:]
classcounts=orthos['DB Class Key'].value_counts()
one2one=classcounts.index[list(classcounts==2)]
orthos=orthos.loc[orthos['DB Class Key'].isin(one2one),:]

htab=orthos.loc[orthos['NCBI Taxon ID']==9606,:]
mtab=orthos.loc[orthos['NCBI Taxon ID']==10090,:]
genemapping=dict(zip([x.upper() for x in mtab['Symbol']],htab['Symbol']))
print(len(genemapping.keys()),flush=True)
print(htab)
print(mtab,flush=True)
paldict={'CGE_NR2F2/PROX1': 'slategray',
    'G1-phase_SLC1A3/ATP1A1': '#17344c',
    'G2-M_UBE2C/ASPM': '#19122b',
    'LGE_FOXP1/ISL1': 'cyan',
    'LGE_FOXP1/PENK': 'navy',
    'LGE_FOXP2/TSHZ1': 'goldenrod',
    'LGE_MEIS2/PAX6': 'orangered',
    'LGE_MEIS2/PAX6/SCGN': 'orange',
    'LGE-OB_MEIS2/PAX6': 'red',
    'MGE_CRABP1/MAF': 'indigo',
    'MGE_CRABP1/TAC3': 'fuchsia',
    'MGE_LHX6/MAF': 'darksalmon',
    'MGE_LHX6/NPY': 'maroon',
    'RMTW_ZIC1/RELN': 'yellow',
    'S-phase_MCM4/H43C': 'lawngreen',
    'Transition': '#3c7632',
    'VMF_ZIC1/ZIC2': 'teal',
    'VMF_CRABP1/LHX8':'skyblue',
    'VMF_NR2F2/LHX6':'lightseagreen',
    'VMF_LHX1/POU6F2':'seagreen',
    'VMF_TMEM163/OTP':'lightsteelblue',
    'VMF_PEG10/DLK1':'steelblue',
    'LGE_FOXP1/ISL1/NPY1R':'mediumpurple',
    'nan':'white',
    'Cholinergic':'olivedrab',
    'Amy/Hypo_HAP1':'darkslateblue',
    'GP_GBX1':'teal',
    'NAc_HAP1':'darkslategray',
    'Excitatory':'whitesmoke',
    'Lamp5':'darkgoldenrod', 
    'Lamp5 Lhx6':'rosybrown', 
    'OB-GC_RPRM':'palegoldenrod', 
    'OB-GC_STXBP6/PENK':'olive', 
    'OB-PGC_FOXP2/CALB1':'aquamarine', 
    'OB-PGC_TH/SCGN':'orange', 
    'OB-PGC_ZIC':'saddlebrown', 
    'Pvalb':'hotpink', 
    'Sncg':'lightslategray', 
    'Sst':'salmon', 
    'Sst Chodl':'maroon',
    'Vip':'tan',
    'PAX6':'sienna',
    'Pvalb Vipr2':'lightcoral',
    'eSPN':'springgreen', 
    'dSPN':'cyan', 
    'iSPN':'navy',
    'v-dSPN':'violet',
    'Str_IN':'black',
    'Glia':'lightgoldenrodyellow'}


# In[3]:


qdata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityMacaqueGeAllcortexHippocampusProcessed.h5ad')
qdata.X=qdata.raw.X[:,qdata.raw.var.index.isin(qdata.var.index)]#.todense()
qdata.obs['supervised_name']=qdata.obs['noctx_supervised_name']

#hdata=sc.read(os.path.join('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityHumanMotor.h5ad'))
#mdata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/MouseBrain.h5ad')
#mdata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityMouseGEClean.h5ad')
mdata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityMouseWbGeProcessed.h5ad')
mdata.X=mdata.raw.X[:,mdata.raw.var.index.isin(mdata.var.index)]#.todense()
print(mdata.var)
mdata=mdata[:,mdata.var.index.isin(genemapping.keys())]


mdata.var.index=[x.upper() for x in mdata.var.index]
qdata.var.index=[x.upper() for x in qdata.var.index]
mdata.var.index=[genemapping[x] for x in mdata.var.index]
qdata.obs['species']='Macaque'
mdata.obs['species']='Mouse'
qdata.var_names_make_unique()
mdata.var_names_make_unique()
qdata.obs_names_make_unique()
mdata.obs_names_make_unique()
print(mdata.var.index)
print(qdata.var.index)
adata=sc.AnnData.concatenate(qdata,mdata,batch_key='species_batch',index_unique=None)

adata.raw=adata
adata.var_names_make_unique()
adata.obs_names_make_unique()
ribo_genes=[name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL') ]
adata.obs['percent_ribo'] = np.sum(
adata[:, ribo_genes].X, axis=1) / np.sum(adata.X, axis=1)
mito_genes = [name for name in adata.var_names if name in ['ND1','ND2','ND4L','ND4','ND5','ND6','ATP6','ATP8','CYTB','COX1','COX2','COX3'] or name.startswith('CHRM-') or name.startswith('MT-')]


# In[4]:


print(adata)
#sc.pl.violin(adata,groupby='batch_name',keys=['percent_ribo','percent_mito'],save='ribomito')

#adata=adata[adata.obs['percent_ribo']<.5,:]
#adata=adata[adata.obs['percent_ribo']<.2,:]
adata._inplace_subset_obs(adata.obs['percent_ribo']<.4)
adata._inplace_subset_obs(adata.obs['percent_mito']<.15)
adata._inplace_subset_var(~adata.var.index.isin(mito_genes))


adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)]
sc.pp.filter_genes(adata,min_cells=10)
sc.pp.normalize_total(adata,exclude_highly_expressed=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=5000,batch_key='species')
#markers=pd.read_csv(os.path.expanduser('~/markers.txt'),sep='\t')
##adata=adata[:,adata.var.index.isin(markers['hgnc_symbol'])]
#adata.var['highly_variable']=list(adata.var.index.isin(markers['hgnc_symbol']))
adata.var['highly_variable']=(adata.var['highly_variable']&(adata.var['highly_variable_nbatches']>1))
sc.pp.scale(adata,max_value=10)


# In[5]:


adata.var.index[adata.var['highly_variable']]


# In[6]:


df=pd.DataFrame(adata[:,adata.var['highly_variable']].X)
df['supervised_name']=list(adata.obs.supervised_name)


# In[7]:


df


# In[ ]:





# In[8]:


df1=df.loc[list(adata.obs.species.isin(['Mouse'])),:].groupby(['supervised_name']).mean()
corr = 1 - df1.T.corr() 
corr_condensed = scipy.cluster.hierarchy.distance.squareform(corr) # convert to condensed
df1_link = scipy.cluster.hierarchy.linkage(corr_condensed, method='average')
df2=df.loc[list(adata.obs['species'].isin(['Macaque'])),:].groupby(['supervised_name']).mean()
corr = 1 - df2.T.corr() 
corr_condensed = scipy.cluster.hierarchy.distance.squareform(corr) # convert to condensed
df2_link = scipy.cluster.hierarchy.linkage(corr_condensed, method='average')

corrmat=pd.DataFrame(corr2_coeff(np.array(df1),np.array(df2)),index=df1.index,columns=df2.index)


# In[11]:


get_ipython().run_line_magic('matplotlib', 'inline')
xind=list(corrmat.index)
print(xind)
yind=list(corrmat.columns)
print(yind)
xind=list(corrmat.index)
xind.remove('Transition')
xind.insert(0, 'Transition')
xind.remove('G2-M_UBE2C/ASPM')
xind.insert(0, 'G2-M_UBE2C/ASPM')
xind.remove('S-phase_MCM4/H43C')
xind.insert(0, 'S-phase_MCM4/H43C')
xind.remove('LGE-OB_MEIS2/PAX6')
xind.insert(8, 'LGE-OB_MEIS2/PAX6')
yind=list(corrmat.columns)
yind.remove('Transition')
yind.insert(0, 'Transition')
yind.remove('G2-M_UBE2C/ASPM')
yind.insert(0, 'G2-M_UBE2C/ASPM')
yind.remove('S-phase_MCM4/H43C')
yind.insert(0, 'S-phase_MCM4/H43C')
yind.remove('G1-phase_SLC1A3/ATP1A1')
yind.insert(0, 'G1-phase_SLC1A3/ATP1A1')
seaborn.heatmap(corrmat.loc[xind,yind],xticklabels=yind,yticklabels=xind,cmap='coolwarm')
seaborn.set(font_scale=.6)
plt.savefig(os.path.join(sc.settings.figdir,'MouseVsMacaqueCorrNoclustermap.pdf'), bbox_inches="tight")


# In[ ]:





# In[2]:


adata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbProcessedQvMMotor.h5ad')


# In[3]:


get_ipython().run_line_magic('matplotlib', 'inline')
adata


# In[4]:


import scipy.spatial as sp, scipy.cluster.hierarchy as hc
seaborn.set(font_scale=.7)
df=pd.DataFrame(adata[:,adata.var['highly_variable']].X)
df['noctx_supervised_name']=list(adata.obs.noctx_supervised_name)

def corr2_coeff(A, B):
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]
    # Sum of squares across rows
    ssA = (A_mA**2).sum(1)
    ssB = (B_mB**2).sum(1)
    # Finally get corr coeff
    return np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:, None],ssB[None]))

df1=df.loc[list(((~adata.obs.region.isin(['ob','OB']) &  ~adata.obs.noctx_supervised_name.str.contains('"')) | adata.obs.noctx_supervised_name.str.contains('OB')) & adata.obs.species.isin(['Mouse'])),:].groupby(['noctx_supervised_name']).mean()
corr = 1 - df1.T.corr() 
corr_condensed = hc.distance.squareform(corr) # convert to condensed
df1_link = hc.linkage(corr_condensed, method='average')
df2=df.loc[list(adata.obs['species'].isin(['Macaque'])),:].groupby(['noctx_supervised_name']).mean()
corr = 1 - df2.T.corr() 
corr_condensed = hc.distance.squareform(corr) # convert to condensed
df2_link = hc.linkage(corr_condensed, method='average')

corrmat=pd.DataFrame(corr2_coeff(np.array(df1),np.array(df2)),index=df1.index,columns=df2.index)
ax =seaborn.clustermap(corrmat,xticklabels=True,yticklabels=True,cmap='coolwarm',row_linkage=df1_link, col_linkage=df2_link)
#plt.ylabel('Mouse Cluster')
#plt.xlabel('Developing Macaque Cluster')
plt.savefig(os.path.join(sc.settings.figdir,'MouseVsMacaqueCorr.pdf'), bbox_inches="tight")


# In[5]:


corrmat


# In[8]:


xind=list(corrmat.index)
print(xind)
yind=list(corrmat.columns)
print(yind)


# In[9]:


corrmat=pd.DataFrame(corr2_coeff(np.array(df1),np.array(df2)),index=df1.index,columns=df2.index)
xind=list(corrmat.index)
xind.remove('Transition')
xind.insert(0, 'Transition')
xind.remove('G2-M_UBE2C/ASPM')
xind.insert(0, 'G2-M_UBE2C/ASPM')
xind.remove('S-phase_MCM4/H43C')
xind.insert(0, 'S-phase_MCM4/H43C')
yind=list(corrmat.columns)
yind.remove('Transition')
yind.insert(0, 'Transition')
yind.remove('G2-M_UBE2C/ASPM')
yind.insert(0, 'G2-M_UBE2C/ASPM')
yind.remove('S-phase_MCM4/H43C')
yind.insert(0, 'S-phase_MCM4/H43C')
yind.remove('G1-phase_SLC1A3/ATP1A1')
yind.insert(0, 'G1-phase_SLC1A3/ATP1A1')
seaborn.heatmap(corrmat.loc[xind,yind],xticklabels=True,yticklabels=True,cmap='coolwarm')
seaborn.set(font_scale=.4)
plt.savefig(os.path.join(sc.settings.figdir,'MouseVsMacaqueCorrNoclustermap.pdf'), bbox_inches="tight")


# In[14]:


ykeep=['LGE_MEIS2/PAX6','LGE_FOXP2/TSHZ1','LGE_FOXP1/ISL1','LGE_FOXP1/PENK']
xkeep=['LGE-OB_MEIS2/PAX6','LGE_MEIS2/PAX6','LGE_FOXP2/TSHZ1','LGE_FOXP1/ISL1','LGE_FOXP1/PENK']

seaborn.heatmap(corrmat.loc[xkeep,ykeep],xticklabels=True,yticklabels=True,cmap='coolwarm')
seaborn.set(font_scale=.9)
plt.savefig(os.path.join(sc.settings.figdir,'MouseVsMacaqueCorrNoclustermapLGEonly.pdf'), bbox_inches="tight")


# In[6]:


import scipy.spatial as sp, scipy.cluster.hierarchy as hc
xdata=adata[adata.obs.noctx_supervised_name.str.contains('LGE|CGE'),:][:,adata.var['highly_variable']]
df=pd.DataFrame(xdata.X)
df['noctx_supervised_name']=list(xdata.obs.noctx_supervised_name)

def corr2_coeff(A, B):
    # Rowwise mean of input arrays & subtract from input arrays themeselves
    A_mA = A - A.mean(1)[:, None]
    B_mB = B - B.mean(1)[:, None]
    # Sum of squares across rows
    ssA = (A_mA**2).sum(1)
    ssB = (B_mB**2).sum(1)
    # Finally get corr coeff
    return np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:, None],ssB[None]))

df1=df.loc[list(((~xdata.obs.region.isin(['ob','OB']) &  ~xdata.obs.noctx_supervised_name.str.contains('"')) | xdata.obs.noctx_supervised_name.str.contains('OB')) & xdata.obs.species.isin(['Mouse'])),:].groupby(['noctx_supervised_name']).mean()
corr = 1 - df1.T.corr() 
corr_condensed = hc.distance.squareform(corr) # convert to condensed
df1_link = hc.linkage(corr_condensed, method='average')
df2=df.loc[list(xdata.obs['species'].isin(['Macaque'])),:].groupby(['noctx_supervised_name']).mean()
corr = 1 - df2.T.corr() 
corr_condensed = hc.distance.squareform(corr) # convert to condensed
df2_link = hc.linkage(corr_condensed, method='average')

corrmat=pd.DataFrame(corr2_coeff(np.array(df1),np.array(df2)),index=df1.index,columns=df2.index)
seaborn.heatmap(corrmat,xticklabels=True,yticklabels=True,cmap='coolwarm')
seaborn.set(font_scale=.4)
plt.savefig(os.path.join(sc.settings.figdir,'MouseVsMacaqueCorrXGENoclustermap.pdf'), bbox_inches="tight")
plt.close()

corrmat=pd.DataFrame(corr2_coeff(np.array(df1),np.array(df2)),index=df1.index,columns=df2.index)
seaborn.clustermap(corrmat,xticklabels=True,yticklabels=True,cmap='coolwarm',row_linkage=df1_link, col_linkage=df2_link)
#plt.ylabel('Mouse Cluster')
#plt.xlabel('Developing Macaque Cluster')
plt.savefig(os.path.join(sc.settings.figdir,'MouseVsMacaqueXGECorr.pdf'), bbox_inches="tight")


# In[7]:


sc.settings.file_format_figs='png'

sc.pl.violin(adata[adata.obs.noctx_supervised_name.str.contains('LHX6/NPY'),:],keys=['NPY','CORT','CHODL','SATB2'],groupby='species',use_raw=False,save='MouseVMacaqueCHODL')


# In[8]:


sc.pl.violin(adata[adata.obs.noctx_supervised_name.str.contains('NR2F2/PROX1'),:],keys=['PROX1','VIP','SATB2'],groupby='species',use_raw=False,save='MouseVMacaqueSATB2')


# In[9]:


sc.pl.stacked_violin(adata[adata.obs.noctx_supervised_name.str.contains('MEIS2/PAX6'),:],var_names=['PAX6','MEIS2','NR2F2','ETV1','FOXP2','TSHZ1'],groupby='species',use_raw=False,save='MouseVMacaqueGC-like')


# In[10]:


matplotlib.rcParams['figure.dpi'] = 100
sc.pl.violin(adata[adata.obs.noctx_supervised_name.str.contains('FOXP2/TSHZ1'),:],keys=['FOXP2','TSHZ1','CASZ1','OPRM1'],groupby='species',use_raw=False,save='MouseVMacaqueCASZ1')
matplotlib.rcParams['figure.dpi'] = 300
sc.settings.file_format_figs='pdf'


# In[11]:


sc.pl.stacked_violin(adata[adata.obs.noctx_supervised_name.str.contains('FOXP2/TSHZ1'),:],var_names=['FOXP2','TSHZ1','CASZ1','OPRM1'],groupby='species',use_raw=False,save='MouseVMacaqueCASZ1')


# In[71]:


newdf=pd.read_csv(os.path.join('/wynton/group/ye/mtschmitz/figures/macaqueBranchDE/',"CoefsTable.csv"),index_col=0)
cols=newdf.columns[['__' in x for x in newdf.columns]]
cols=[x.split('__')[0] for x in cols]


# In[72]:


macaquecoefdfs={}
for c in cols:
    macaquecoefdfs[c]=newdf.loc[:,newdf.columns.str.contains(c)]
    macaquecoefdfs[c].columns=['coef','intercept','pval','gene','log_reg_coef','signif','qval','log10qval']
    macaquecoefdfs[c].index=macaquecoefdfs[c]['gene']


# In[73]:


newdf=pd.read_csv(os.path.join('/wynton/group/ye/mtschmitz/figures/mouseWbGeCB202002/',"CoefsTable.csv"),index_col=0)
cols=newdf.columns[['__' in x for x in newdf.columns]]
cols=[x.split('__')[0] for x in cols]


# In[74]:


mousecoefdfs={}
for c in cols:
    mousecoefdfs[c]=newdf.loc[:,newdf.columns.str.contains(c)]
    mousecoefdfs[c].columns=['coef','intercept','pval','gene','log_reg_coef','signif','qval','log10qval']
    mousecoefdfs[c].index=mousecoefdfs[c]['gene']


# In[7]:


macaquecoefdfs


# In[7]:


import scipy.spatial as sp, scipy.cluster.hierarchy as hc
seaborn.set(font_scale=.7)
def jaccard_similarity(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union

mouse_activate_geneset={}
for k in mousecoefdfs.keys():
    mouse_activate_geneset[k]=mousecoefdfs[k].index[(mousecoefdfs[k]['coef']>1)]#& (mousecoefdfs[k]['qval']<.05)

macaque_activate_geneset={}
for k in macaquecoefdfs.keys():
    macaque_activate_geneset[k]=macaquecoefdfs[k].index[(macaquecoefdfs[k]['coef']>1) ]#& (macaquecoefdfs[k]['qval']<.05)
activate_setcounts=np.zeros((len(mouse_activate_geneset.keys()),len(macaque_activate_geneset.keys())))

#activate_sums={}
for ini,i in enumerate(sorted(mouse_activate_geneset.keys())):
    for inj,j in enumerate(sorted(macaque_activate_geneset.keys())):
        #activate_setcounts[ini,inj]=len(set(activate_geneset[i]).intersection(set(activate_geneset[j])))
        activate_setcounts[ini,inj]=jaccard_similarity(mouse_activate_geneset[i],macaque_activate_geneset[j])
    #activate_sums[i]=len(activate_geneset[i])
activate_setcounts_dist=1-activate_setcounts
#np.fill_diagonal(activate_setcounts,0)

#corr_condensed = hc.distance.squareform(activate_setcounts_dist) # convert to condensed
#activate_setcounts_link = hc.linkage(corr_condensed, method='average')

mouse_inactivate_geneset={}
for k in mousecoefdfs.keys():
    mouse_inactivate_geneset[k]=mousecoefdfs[k].index[(mousecoefdfs[k]['coef']<-1)]# & (mousecoefdfs[k]['qval']<.05)

macaque_inactivate_geneset={}
for k in macaquecoefdfs.keys():
    macaque_inactivate_geneset[k]=macaquecoefdfs[k].index[(macaquecoefdfs[k]['coef']<-1)]#& (macaquecoefdfs[k]['qval']<.05)
inactivate_setcounts=np.zeros((len(mouse_inactivate_geneset.keys()),len(macaque_inactivate_geneset.keys())))

#inactivate_sums={}
for ini,i in enumerate(sorted(mouse_inactivate_geneset.keys())):
    for inj,j in enumerate(sorted(macaque_inactivate_geneset.keys())):
        #inactivate_setcounts[ini,inj]=len(set(inactivate_geneset[i]).intersection(set(inactivate_geneset[j])))
        inactivate_setcounts[ini,inj]=jaccard_similarity(mouse_inactivate_geneset[i],macaque_inactivate_geneset[j])
    #inactivate_sums[i]=len(inactivate_geneset[i])
inactivate_setcounts_dist=1-inactivate_setcounts
#np.fill_diagonal(inactivate_setcounts,0)
plt.figure()
ax=seaborn.heatmap(pd.DataFrame(activate_setcounts,index=sorted(mouse_activate_geneset.keys()),columns=sorted(macaque_activate_geneset.keys())),cmap=matplotlib.cm.Reds,xticklabels=True,yticklabels=True)#row_linkage=activate_setcounts_link,col_linkage=inactivate_setcounts_link,
ax.xaxis.tick_top() # x axis on top
ax.xaxis.set_label_position('top')
ax.set_xticklabels(ax.get_xticklabels(),rotation=90,horizontalalignment='center')
plt.title('Jaccard Index of Activating Genesets')
plt.xlabel('Macaque')
plt.ylabel('Mouse')
plt.savefig(os.path.join(sc.settings.figdir,'ActivateJaccard.pdf'), bbox_inches="tight")
plt.show()

plt.figure()
ax=seaborn.heatmap(pd.DataFrame(inactivate_setcounts,index=sorted(mouse_activate_geneset.keys()),columns=sorted(macaque_activate_geneset.keys())),cmap=matplotlib.cm.Blues,xticklabels=True,yticklabels=True)#row_linkage=activate_setcounts_link,col_linkage=inactivate_setcounts_link,
ax.xaxis.set_label_position('bottom')
ax.set_xticklabels(ax.get_xticklabels(),rotation=90,horizontalalignment='center')
plt.title('Jaccard Index of Inactivating Genesets')
plt.xlabel('Macaque')
plt.ylabel('Mouse')
plt.savefig(os.path.join(sc.settings.figdir,'InactivateJaccard.pdf'), bbox_inches="tight")
plt.show()


# In[49]:


d={}
for k in mousecoefdfs.keys():
    if (k in mousecoefdfs.keys()) and (k in macaquecoefdfs.keys()):
        inds=set(mousecoefdfs[k].index).intersection(set(macaquecoefdfs[k].index))
        qmat=macaquecoefdfs[k].loc[inds,:]
        mmat=mousecoefdfs[k].loc[inds,:]
        df=pd.DataFrame([mmat['coef'],qmat['coef'],mmat['signif'],qmat['signif']],index=['mouse','macaque','msig','qsig']).T
        df['msig']=df['msig'].astype(bool)
        df['qsig']=df['qsig'].astype(bool)
        sigo=[]
        for i in df.index:
            if df.loc[i,'msig'] & df.loc[i,'qsig']& (np.sign(df.loc[i,'mouse'])==np.sign(df.loc[i,'macaque'])):
                sigo.append('M&Q')
            elif df.loc[i,'msig'] & df.loc[i,'qsig']& (np.sign(df.loc[i,'mouse'])!=np.sign(df.loc[i,'macaque'])):
                sigo.append('Flip')
            elif df.loc[i,'msig'] & ~df.loc[i,'qsig']:
                sigo.append('M')
            elif ~df.loc[i,'msig'] & df.loc[i,'qsig']:
                sigo.append('Q')
            else:
                sigo.append('None')
        df['signif']=sigo
        d[k]=df
        seaborn.lmplot(x='mouse',y='macaque',data=df,hue='signif',fit_reg=False)
        plt.title(k)
        plt.show()


# In[52]:


import matplotlib.pyplot as plt
from matplotlib_venn import venn2

for k in mousecoefdfs.keys():
    if (k in mousecoefdfs.keys()) and (k in macaquecoefdfs.keys()):
        inds=set(mousecoefdfs[k].index).intersection(set(macaquecoefdfs[k].index))
        qmat=macaquecoefdfs[k].loc[inds,:]
        mmat=mousecoefdfs[k].loc[inds,:]
        df=pd.DataFrame([mmat['coef'],qmat['coef'],mmat['signif'],qmat['signif']],index=['mouse','macaque','msig','qsig']).T
        df['msig']=df['msig'].astype(bool)
        df['qsig']=df['qsig'].astype(bool)
        print(df.index[df['msig']])
        print(df.index[df['qsig']])
        plt.figure()
        venn2([set(df.index[df['msig']]), set(df.index[df['qsig']])],set_labels=['Mouse','Macaque'])
        plt.title(k)
        plt.show()
        


# ### Figure 3 Scatterplots

# In[76]:


newdf=pd.read_csv(os.path.join('/wynton/group/ye/mtschmitz/figures/macaqueGeAllcortexHippoCB202002/',"DiffxpyMarkers.csv"),index_col=0)
usecols=newdf.columns[['__' in x for x in newdf.columns]]
cols=[x.split('__')[0] for x in usecols]


# In[77]:


macaquecoefdfs={}
for c in set(cols):
    titles=[x.split('__')[1] for x in usecols[usecols.str.contains(c)]]
    macaquecoefdfs[c]=newdf.loc[:,newdf.columns.str.contains(c)]
    macaquecoefdfs[c].columns=titles
    macaquecoefdfs[c].index=macaquecoefdfs[c]['gene']


# In[78]:


newdf=pd.read_csv(os.path.join('/wynton/group/ye/mtschmitz/figures/mouseWbGeCB202002/',"DiffxpyMarkers.csv"),index_col=0)
usecols=newdf.columns[['__' in x for x in newdf.columns]]
cols=[x.split('__')[0] for x in usecols]


# In[79]:


mousecoefdfs={}
for c in set(cols):
    titles=[x.split('__')[1] for x in usecols[np.array([x.split('__')[0] for x in usecols])==c]]
    mousecoefdfs[c]=newdf.loc[:,np.array([x.split('__')[0] for x in newdf.columns])==c]
    mousecoefdfs[c].columns=titles
    mousecoefdfs[c].index=mousecoefdfs[c]['gene']


# In[80]:


mousecoefdfs[k]


# In[81]:


get_ipython().run_line_magic('matplotlib', 'inline')
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
seaborn.set(font_scale=.7)
def jaccard_similarity(list1, list2):
    intersection = len(list(set(list1).intersection(list2)))
    union = (len(list1) + len(list2)) - intersection
    return float(intersection) / union

mouse_activate_geneset={}
for k in mousecoefdfs.keys():
    print(k)
    mouse_activate_geneset[k]=mousecoefdfs[k].index[(mousecoefdfs[k]['log2fc']>1) & (mousecoefdfs[k]['qval']<.05)]#

macaque_activate_geneset={}
for k in macaquecoefdfs.keys():
    macaque_activate_geneset[k]=macaquecoefdfs[k].index[(macaquecoefdfs[k]['log2fc']>1)& (macaquecoefdfs[k]['qval']<.05) ]#
activate_setcounts=np.zeros((len(mouse_activate_geneset.keys()),len(macaque_activate_geneset.keys())))

#activate_sums={}
for ini,i in enumerate(sorted(mouse_activate_geneset.keys())):
    for inj,j in enumerate(sorted(macaque_activate_geneset.keys())):
        #activate_setcounts[ini,inj]=len(set(activate_geneset[i]).intersection(set(activate_geneset[j])))
        activate_setcounts[ini,inj]=jaccard_similarity(mouse_activate_geneset[i],macaque_activate_geneset[j])
    #activate_sums[i]=len(activate_geneset[i])
activate_setcounts_dist=1-activate_setcounts
#np.fill_diagonal(activate_setcounts,0)

#corr_condensed = hc.distance.squareform(activate_setcounts_dist) # convert to condensed
#activate_setcounts_link = hc.linkage(corr_condensed, method='average')

plt.figure()
ax=seaborn.heatmap(pd.DataFrame(activate_setcounts,index=sorted(mouse_activate_geneset.keys()),columns=sorted(macaque_activate_geneset.keys())),cmap=matplotlib.cm.Reds,xticklabels=True,yticklabels=True)#row_linkage=activate_setcounts_link,col_linkage=inactivate_setcounts_link,
ax.xaxis.tick_top() # x axis on top
ax.xaxis.set_label_position('top')
ax.set_xticklabels(ax.get_xticklabels(),rotation=90,horizontalalignment='center')
plt.title('Jaccard Index of DiffxpyMarkers')
plt.xlabel('Macaque')
plt.ylabel('Mouse')
plt.savefig(os.path.join(sc.settings.figdir,'MarkersJaccard.pdf'), bbox_inches="tight")
plt.show()


# In[82]:


d={}
for k in mousecoefdfs.keys():
    if (k in mousecoefdfs.keys()) and (k in macaquecoefdfs.keys()):
        sanitize_k="".join(x for x in k if x.isalnum())
        print(sanitize_k)
        inds=set(mousecoefdfs[k].index).intersection(set(macaquecoefdfs[k].index))
        qmat=macaquecoefdfs[k].loc[inds,:]
        mmat=mousecoefdfs[k].loc[inds,:]
        df=pd.DataFrame([np.clip(mmat['log2fc'],-4,5),np.clip(qmat['log2fc'],-4,5),mmat['signif'],qmat['signif']],index=['mouse','macaque','msig','qsig']).T
        df['msig']=df['msig'].astype(bool)
        df['qsig']=df['qsig'].astype(bool)
        sigo=[]
        for i in df.index:
            if df.loc[i,'msig'] & df.loc[i,'qsig'] & (np.sign(df.loc[i,'mouse'])>0) & ((np.sign(df.loc[i,'mouse'])==np.sign(df.loc[i,'macaque']))):
                sigo.append(1)
            elif df.loc[i,'msig'] & ~df.loc[i,'qsig'] & (np.sign(df.loc[i,'mouse'])>0) & (np.sign(df.loc[i,'macaque'])<0):
                sigo.append(2)
            elif ~df.loc[i,'msig'] & df.loc[i,'qsig'] & (np.sign(df.loc[i,'macaque'])>0) & (np.sign(df.loc[i,'mouse'])<0):
                sigo.append(3)
            else:
                sigo.append(4)
        df['signif']=sigo
        d[k]=df
        df.signif=df.signif.astype('category')
        df.signif.cat.reorder_categories(sorted(df.signif.cat.categories),inplace=True)
        #seaborn.distplot(np.clip(mmat['log2fc'],-4,5))
        #seaborn.distplot(np.clip(qmat['log2fc'],-4,5))
        #plt.show()
        coolgenes=df.index[(df['mouse']>1.2)&(df['macaque']>1.5)]
        print(list(coolgenes))
        coolgenes=df.index[(df['mouse']>1.2)&~(df['macaque']>1.5)]
        print(list(coolgenes))
        coolgenes=df.index[~(df['mouse']>1.2)&(df['macaque']>1.5)]
        print(list(coolgenes))


# In[83]:


mousecoefdfs['MGE_CRABP1/TAC3']=mousecoefdfs['MGE_CRABP1/MAF']


# In[84]:


get_ipython().run_line_magic('matplotlib', 'inline')
import adjustText

d={}
for k in mousecoefdfs.keys():
    if (k in mousecoefdfs.keys()) and (k in macaquecoefdfs.keys()):
        sanitize_k="".join(x for x in k if x.isalnum())
        print(sanitize_k)
        inds=set(mousecoefdfs[k].index).intersection(set(macaquecoefdfs[k].index))
        qmat=macaquecoefdfs[k].loc[inds,:]
        mmat=mousecoefdfs[k].loc[inds,:]
        df=pd.DataFrame([np.clip(mmat['log2fc'],-4,5),np.clip(qmat['log2fc'],-4,5),mmat['signif'],qmat['signif']],index=['mouse','macaque','msig','qsig']).T
        df['msig']=df['msig'].astype(bool)
        df['qsig']=df['qsig'].astype(bool)
        sigo=[]
        for i in df.index:
            if df.loc[i,'msig'] & df.loc[i,'qsig'] & (np.sign(df.loc[i,'mouse'])>0) & ((np.sign(df.loc[i,'mouse'])==np.sign(df.loc[i,'macaque']))):
                sigo.append(1)
            elif df.loc[i,'msig'] & ~df.loc[i,'qsig'] & (np.sign(df.loc[i,'mouse'])>0) & (np.sign(df.loc[i,'macaque'])<0):
                sigo.append(2)
            elif ~df.loc[i,'msig'] & df.loc[i,'qsig'] & (np.sign(df.loc[i,'macaque'])>0) & (np.sign(df.loc[i,'mouse'])<0):
                sigo.append(3)
            else:
                sigo.append(4)
        df['signif']=sigo
        print(list(df.loc[np.array(sigo)==1,:].index))
        print(list(df.loc[np.array(sigo)==2,:].index))
        print(list(df.loc[np.array(sigo)==3,:].index))
        d[k]=df
        df.signif=df.signif.astype('category')
        df.signif.cat.reorder_categories(sorted(df.signif.cat.categories),inplace=True)
        #seaborn.distplot(np.clip(mmat['log2fc'],-4,5))
        #seaborn.distplot(np.clip(qmat['log2fc'],-4,5))
        #plt.show()
        
        coolgenes=df.index[(df['mouse']>1.2)|(df['macaque']>1.5)]
        tolabel=[x for x in df.index[df['signif'].isin([1,2,3])] if any([i == x for i in coolgenes])]
        coolgenes=['OPR.*','^FOX.*','PENK','^ISL[0-9]+','^PAX.*','^DRD.*','^NKX.*','^HTR.*','^ST18$','^TRH$','^ZIC.*',
                   '^ZFHX.*','^TH$','^SALL.*','FGF.*','^CASZ.*','OTOF','ANGPT2','^BMP.*','^MEF.*','^IGFBP.*','^SPOCK.*',
                   'TGFB.*','WNT.*','^PDYN','EBF.*','IKZF.*','VAX.*','^CASTOR','^ADORA.*','SST','^GABR.*','^SHH$','^RYR.*','^SLIT.*',
                   '^TSHZ.*','^TLE.*','^SP8','^NR2.*','^ADARB.*','LHX.*','PROX.*','^CCK$','^VIP$','MAF','NPY','^STXBP.*','^NPTX.*',
                   'CHODL','SCGN','NXPH.*','^RBP.*','CRABP1','^NOG$','^RXR.*','^RAR.*','^GNG8$','SLC4A4','^ETV.*','^CHRN.*','^TACR.*']
        tolabel=[x for x in tolabel if any([bool(re.search(i,x)) for i in coolgenes])]
        
        seaborn.set(style="whitegrid")
        g=seaborn.relplot(data=df,x='mouse',y='macaque',hue=df['signif'].values,palette=['green','orange','blue','lightgray'])#size=pvalDF['percent_cells'].astype(float).values
        ax = g.axes[0,0]
        leg = g._legend
        leg.set_bbox_to_anchor([1.5, 0.7])  # coordinates of lower left of bounding box
        plt.xlabel('Mouse (Class / Rest ) Log2FC')
        plt.ylabel('Macaque (Class / Rest ) Log2FC')
        texts=[]
        for line in tolabel:
            texts.append(ax.text(df.mouse[line], df.macaque[line], df.index[df.index==line].values[0], horizontalalignment='left', size='small', color='black'))#+np.random.uniform(-.04,.04)
        adjustText.adjust_text(texts,arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
        plt.title(k)
        plt.savefig(os.path.join(sc.settings.figdir,sanitize_k+'_DE_relplot.png'), bbox_inches="tight")
        plt.show()
        #plt.clf()
        #plt.figure(figsize=(20,20))
        #seaborn.lmplot(x='mouse',y='macaque',data=df,hue='signif',fit_reg=False)
        #plt.title(k)
        #plt.show()


# In[85]:


import matplotlib.pyplot as plt
from matplotlib_venn import venn2

for k in mousecoefdfs.keys():
    if (k in mousecoefdfs.keys()) and (k in macaquecoefdfs.keys()):
        inds=set(mousecoefdfs[k].index).intersection(set(macaquecoefdfs[k].index))
        qmat=macaquecoefdfs[k].loc[inds,:]
        mmat=mousecoefdfs[k].loc[inds,:]
        df=pd.DataFrame([mmat['log2fc'],qmat['log2fc'],mmat['signif'],qmat['signif']],index=['mouse','macaque','msig','qsig']).T
        df['msig']=(df['mouse']>1.2)&df['msig'].astype(bool)
        df['qsig']=(df['macaque']>1.5)&df['qsig'].astype(bool)
        print(df.index[df['msig']])
        print(df.index[df['qsig']])
        plt.figure()
        venn2([set(df.index[df['msig']]), set(df.index[df['qsig']])],set_labels=['Mouse','Macaque'])
        plt.title(k)
        plt.show()
        


# In[86]:


mousecoefdfs


# In[ ]:





# In[ ]:





# In[22]:


supercell=pd.read_csv('/wynton/group/ye/mtschmitz/macaquedevbrain/MacaqueGEsupervisednamesHippoPredictedEnd.txt')
ind=adata.obs.index[adata.obs.index.isin(supercell['full_cellname'])]
supercell.index=supercell['full_cellname']
supercell=supercell.loc[ind,:]
adata.obs['latent_time']=np.nan
adata.obs['predicted_end']=''
adata.obs.loc[ind,['noctx_supervised_name','predicted_end','latent_time']]=supercell.loc[:,['noctx_supervised_name','predicted_end','latent_time']]
adata.obs.loc[adata.obs['predicted_end'].isna(),'noctx_supervised_name']=adata.obs['noctx_supervised_name'][adata.obs['predicted_end'].isna()]
adata.obs['supervised_name']=adata.obs['noctx_supervised_name']

supercell=pd.read_csv('/wynton/group/ye/mtschmitz/macaquedevbrain/MouseGEsupervisednamesPredictedEnd.txt')
ind=adata.obs.index[adata.obs.index.isin(supercell['full_cellname'])]
supercell.index=supercell['full_cellname']
supercell=supercell.loc[ind,:]
adata.obs.loc[ind,['noctx_supervised_name','predicted_end','latent_time']]=supercell.loc[:,['noctx_supervised_name','predicted_end','latent_time']]
adata.obs.loc[adata.obs['predicted_end'].isna(),'noctx_supervised_name']=adata.obs['noctx_supervised_name'][adata.obs['predicted_end'].isna()]
adata.obs['supervised_name']=adata.obs['noctx_supervised_name']
adata.obs['supervised_name']=adata.obs['supervised_name'].replace({'LGE-OB_MEIS2/PAX6_GC':'LGE_MEIS2/PAX6_PC'})


# In[23]:


adata.obs['supervised_name'].unique()


# In[124]:


#number of windows
n=100

allgenes=[]
celldict={}
celllist=[]
for k in d.keys():
    print(k)
    genes=list(d[k].index[(d[k]['signif']=='M&Q')]) + list(d[k].index[(d[k]['signif']=='M')]) +list(d[k].index[(d[k]['signif']=='Q')])
    allgenes=allgenes+genes
    mdata=adata[(adata.obs.species=='Mouse') & (adata.obs['supervised_name']==k),:]
    celllist=celllist+list(mdata.obs.index[mdata.obs['latent_time'].to_numpy().argsort()])
    celldict[k]=np.array_split(np.array(mdata.obs.index[mdata.obs['latent_time'].to_numpy().argsort()]),n)


# In[125]:


len(allgenes)


# In[172]:


#number of windows
n=500

allgenes=[]
celldict={}
celllist=[]

for k in d.keys():
    print(k)
    mdata=adata[(adata.obs.species=='Mouse') & (adata.obs['supervised_name']==k),:]
    l=[]
    for x in ['M&Q','M','Q']:
        genes=d[k].index[(d[k]['signif']==x)]
        mmat=mdata[:,genes].X[mdata.obs['latent_time'].to_numpy().argsort(),:]
        #values=[]
        #for x in np.array_split(list(range(mmat.shape[0])),n):
        #    if len(x)>5:
        #        values.append(mmat[x,:].mean(0))
        #mmat=np.array(values)
        mmat=mmat[:,mmat.argmax(0).argsort()]
        allgenes=allgenes+list(genes[mmat.argmax(0).argsort()])
        l.append(mmat)
    celllist=celllist+list(mdata.obs.index[mdata.obs['latent_time'].to_numpy().argsort()])
    celldict[k]=np.array_split(np.array(mdata.obs.index[mdata.obs['latent_time'].to_numpy().argsort()]),n)
    plt.figure()
    seaborn.heatmap(np.concatenate(l,axis=1),cmap='coolwarm')
    plt.title(k)
    plt.show()


# In[173]:


#number of windows
n=500

allgenes=[]
celldict={}
celllist=[]

for k in d.keys():
    print(k)
    mdata=adata[(adata.obs.species=='Macaque') & (adata.obs['supervised_name']==k),:]
    l=[]
    for x in ['M&Q','M','Q']:
        genes=d[k].index[(d[k]['signif']==x)]
        mmat=mdata[:,genes].X[mdata.obs['latent_time'].to_numpy().argsort(),:]
        #values=[]
        #for x in np.array_split(list(range(mmat.shape[0])),n):
        #    if len(x)>5:
        #        values.append(mmat[x,:].mean(0))
        #mmat=np.array(values)
        mmat=mmat[:,mmat.argmax(0).argsort()]
        allgenes=allgenes+list(genes[mmat.argmax(0).argsort()])
        l.append(mmat)
    celllist=celllist+list(mdata.obs.index[mdata.obs['latent_time'].to_numpy().argsort()])
    celldict[k]=np.array_split(np.array(mdata.obs.index[mdata.obs['latent_time'].to_numpy().argsort()]),n)
    plt.figure()
    seaborn.heatmap(np.concatenate(l,axis=1),cmap='coolwarm')
    plt.title(k)
    plt.show()


# In[129]:


adata[celllist,:].obs.index


# In[74]:


genes=d[k].index[d[k]['signif']==x]
mmat=mdata[:,genes].X[mdata.obs['latent_time'].to_numpy().argsort(),:]


# In[85]:


.shape


# In[84]:


.shape


# In[58]:


mouse_cells=((~adata.obs.region.isin(['ob','OB']) &  ~adata.obs.noctx_supervised_name.str.contains('"')) | adata.obs.noctx_supervised_name.str.contains('OB')) & adata.obs.species.isin(['Mouse'])
macaque_cells=adata.obs['species'].isin(['Macaque'])
ob_cells=adata.obs.region.isin(['ob','OB']) 


# In[71]:


gcdata=adata[mouse_cells&adata.obs.supervised_name.str.contains('MEIS2/PAX6'),:].copy()


# In[81]:


import diffxpy
import diffxpy.api as de
covarname='supervised_name'
gcdata.X=gcdata.raw.X[:,gcdata.var.index.isin(gcdata.var.index)]
sf=gcdata.X.mean(axis=1)
sf = sf/sf.mean()
gcdata.obs['sf']=sf
test = de.test.wald(
    data=gcdata,
    formula_loc='~ 1 + supervised_name',
    factor_loc_totest='supervised_name',
    noise_model="nb",
    size_factors='sf')


# In[123]:


resultsdiffx=pd.DataFrame([test.gene_ids,test.log2_fold_change(),test.qval,-test.log10_qval_clean()]).T
resultsdiffx.columns=['gene','log2fc','qval','log10qval']
resultsdiffx['log2fc']=np.clip(resultsdiffx['log2fc'],-5,5)
resultsdiffx['signif']=(resultsdiffx['qval']<.01) & (np.abs(resultsdiffx['log2fc'])>2) & (resultsdiffx['log10qval']<25)
resultsdiffx=resultsdiffx.loc[np.argsort(resultsdiffx['log2fc']),:]


# In[115]:


seaborn.scatterplot(x='log2fc',y='log10qval',data=resultsdiffx)


# In[124]:


resultsdiffx.loc[resultsdiffx.signif,:]


# In[121]:


sc.pl.umap(gcdata,color=['supervised_name'])


# In[125]:


sc.pl.umap(gcdata,color=['IQGAP1','PDZD11','FUT1','RARS2','SMOC1'])


# In[126]:


comparables=adata[~ob_cells,:].copy()


# In[139]:


clustdata.obs['species']


# In[128]:


clustdata=comparables[comparables.obs[covarname]==c,:].copy()
print(clustdata.obs['species'].value_counts())
sf=clustdata.X.mean(axis=1)
sf = sf/sf.mean()
clustdata.obs['sf']=sf
test = de.test.wald(
    data=clustdata,
    formula_loc='~ 1 + species',
    factor_loc_totest='species',
    noise_model="nb",
    size_factors='sf')


# In[ ]:


comparables.obs['protocol']=['v3' if '2019' in x else 'v2' for x in comparables.obs['batch_name']]
for c in comparables.obs[covarname].unique():
    print(c)
    clustdata=comparables[comparables.obs[covarname]==c,:].copy()
    print(clustdata.obs['species'].value_counts())
    try:
        sf=clustdata.X.mean(axis=1)
        sf = sf/sf.mean()
        clustdata.obs['sf']=sf
        test = de.test.wald(
            data=clustdata,
            formula_loc='~ 1 + species',
            factor_loc_totest='species',
            noise_model="nb",
            size_factors='sf')
        np.save('/wynton/scratch/mtschmitz/'+covarname+'CrossSpeciesDE_'+c+'.npy',test,allow_pickle=True)
    except:
        pass


# In[4]:


adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)]
adata.obs['n_counts']=np.log10(adata.X.sum(1))
seaborn.violinplot(data=adata.obs,x='species',y='n_counts')
plt.title('before')
plt.show()
sc.pp.downsample_counts(adata,counts_per_cell=600)
adata.obs['n_counts']=np.log10(adata.X.sum(1))
seaborn.violinplot(data=adata.obs,x='species',y='n_counts')
plt.title('after')
plt.show()
sc.pp.filter_genes(adata,min_cells=10)
sc.pp.normalize_total(adata,exclude_highly_expressed=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=1000,batch_key='species')
adata.var['highly_variable']=(adata.var['highly_variable']&(adata.var['highly_variable_nbatches']>1))
print(adata.X.shape,flush=True)
sc.pp.scale(adata,max_value=10)


# In[6]:


lineage_genes= ['DCX','ERBB4','GAD1','GAD2','GADD45G','LHX6','LHX8','SST','NKX2-1','MAF','SP8','SP9','PROX1','VIP','CCK','LAMP5','NR2F1','NR2F2','TOX1','TOX3','ETV1','SCGN','FOXP1','FOXP2','FOXP4','TAC1','TAC3','NPY','PCP4','PAX6','VAX1','MEIS2','ISL1','PENK','ADORA2A','DRD1','DRD2','EBF1','TLE4','ZIC1','ZIC2','ZIC4','GSX2','RELN']
s_genes=['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8']
g2m_genes=['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA']
r_genes=list(np.random.choice(adata.var.index,size=100,replace=False))
allgenes=lineage_genes+s_genes+g2m_genes+r_genes
allgenes=[x for x in allgenes if x in adata.var.index]


# In[7]:


cols={'l':'purple','s':'cyan','g':'orange','r':'grey'}


# In[8]:


genecategories=[]
for g in allgenes:
    if g in lineage_genes:
        genecategories.append('l')
    elif g in s_genes:
        genecategories.append('s')
    elif g in g2m_genes:
        genecategories.append('g')
    else:
        genecategories.append('r')
genecolors=[cols[x] for x in genecategories]
genecategories=np.array(genecategories)


# In[9]:


# Calculate marker expression for progenitors
mdata=adata[adata.obs.species.str.contains('ouse')&adata.obs.noctx_supervised_name.str.contains('ividing'),:].copy()
#sc.pp.downsample_counts(mdata,counts_per_cell=600)
mmat=pd.DataFrame(mdata[:,allgenes].X,index=mdata.obs.index,columns=allgenes)
mcorr=mmat.corr()
mcorr.fillna(0,inplace=True)
seaborn.set(font_scale=.5)
seaborn.clustermap(mcorr,yticklabels=True,xticklabels=True,cmap='coolwarm',col_colors=genecolors,row_colors=genecolors,vmin=-1, vmax=1)
plt.savefig(os.path.join(sc.settings.figdir,'MouseDividingCellCorrelationsClustermap.pdf'), bbox_inches="tight")
plt.show()


# In[10]:


qdata=adata[adata.obs.species.str.contains('caque')&adata.obs.noctx_supervised_name.str.contains('ividing'),:].copy()
#sc.pp.downsample_counts(qdata,counts_per_cell=600)
qmat=pd.DataFrame(qdata[:,allgenes].X,index=qdata.obs.index,columns=allgenes)
qcorr=qmat.corr()
qcorr.fillna(0,inplace=True)
seaborn.set(font_scale=.5)
seaborn.clustermap(qcorr,yticklabels=True,xticklabels=True,cmap='coolwarm',col_colors=genecolors,row_colors=genecolors,vmin=-1, vmax=1)
plt.savefig(os.path.join(sc.settings.figdir,'MacaqueDividingCellCorrelationsClustermap.pdf'), bbox_inches="tight")
plt.show()


# In[11]:


qdata.obs


# In[12]:


mcorr


# In[13]:


qcorr


# In[14]:


d={}
for c in set(genecategories):
    subcor=qcorr.loc[genecategories==c,genecategories==c].to_numpy()
    np.fill_diagonal(subcor,np.nan)
    d[c]=np.nanmax(np.abs(subcor),axis=1).flatten()

qdfs=[]
for c in d.keys():
    flattened_df=pd.DataFrame(d[c],columns=['value'])
    flattened_df['categor']=c
    flattened_df['species']='q'
    qdfs.append(flattened_df)
    
d={}
for c in set(genecategories):
    subcor=mcorr.loc[genecategories==c,genecategories==c].to_numpy()
    np.fill_diagonal(subcor,np.nan)
    d[c]=np.nanmax(np.abs(subcor),axis=1).flatten()

mdfs=[]
for c in d.keys():
    flattened_df=pd.DataFrame(d[c],columns=['value'])
    flattened_df['categor']=c
    flattened_df['species']='m'
    mdfs.append(flattened_df)


# In[15]:


seaborn.violinplot(x='categor',y='value',hue='species',data=pd.concat(mdfs+qdfs))
plt.title('Dividing cell correlations')
plt.ylabel('Marker highest correlation')
plt.show()
catdf=pd.concat(mdfs+qdfs)
seaborn.violinplot(x='categor',y='value',hue='species',data=catdf.loc[catdf['categor'].str.contains('l'),:])
plt.show()
print(catdf.groupby(['categor','species']).mean())


# In[16]:


seaborn.heatmap(qcorr-mcorr,yticklabels=True,xticklabels=True,cmap='coolwarm',vmin=-1, vmax=1)
plt.show()


# In[107]:


qdata.obs['n_counts']=np.log10(qdata.X.sum(1))
seaborn.violinplot(data=qdata.obs,x='species',y='n_counts')


# In[4]:


import memento


# In[5]:


divadata=adata[adata.obs.noctx_supervised_name.str.contains('ividing'),:].copy()


# In[9]:


divadata.X=divadata.raw.X[:,divadata.raw.var.index.isin(divadata.var.index)]


# In[10]:


divadata.obs['capture_rate'] = 0.07
divadata.obs.loc[divadata.obs.species.str.contains('que'),'capture_rate']=0.15
memento.setup_memento(divadata, q_column='capture_rate')
memento.create_groups(divadata, label_columns=['species'])
memento.compute_1d_moments(divadata,
    min_perc_group=.9) # percentage of groups that satisfy the condition for a gene to be considered.


# In[14]:


memento.


# In[ ]:


import itertools
lineage_genes= ['DCX','ERBB4','GAD1','GAD2','GADD45G','LHX6','LHX8','SST','NKX2-1','MAF','SP8','SP9','PROX1','VIP','CCK','LAMP5','NR2F1','NR2F2','TOX1','TOX3','ETV1','SCGN','FOXP1','FOXP2','FOXP4','TAC1','TAC3','NPY','PCP4','PAX6','VAX1','TSHZ1','CASZ1','MEIS2','ISL1','PENK','ADORA2A','DRD1','DRD2','EBF1','TLE4','ZIC1','ZIC2','ZIC4','GSX2','RELN']
lineage_genes=[x for x in lineage_genes if x in divadata.var.index]
gene_pairs = list(itertools.product(lineage_genes,lineage_genes))

memento.compute_2d_moments(divadata, gene_pairs)
memento.ht_2d_moments(
    divadata, 
    formula_like='1', 
    cov_column='', 
    num_cpus=1, 
    num_boot=50000)

result_2d = memento.get_2d_ht_result(divadata)

result_2d.sort_values('corr_pval').head(10)


# In[19]:


seaborn.distplot(result_2d['corr_coef'])


# In[104]:


zicdata=adata[adata.obs.noctx_supervised_name.str.contains('ZIC1\\+ Neuron|NR2F2\\+|LHX8'),:].copy()#|PAX6|adata.obs.region.str.contains('POA|eptum|ypo')
zicdata.obs


# In[105]:


zicdata.X=zicdata.raw.X[:,zicdata.raw.var.index.isin(adata.var.index)]
sc.pp.filter_genes(adata,min_cells=20)
sc.pp.normalize_total(zicdata,exclude_highly_expressed=True)
sc.pp.log1p(zicdata)
sc.pp.highly_variable_genes(zicdata,n_top_genes=1000,batch_key='species')
zicdata.var['highly_variable']=(zicdata.var['highly_variable']&(zicdata.var['highly_variable_nbatches']>1))
sc.pp.scale(zicdata,max_value=10)
sc.pp.pca(zicdata,n_comps=50)


# In[106]:


bbknn.bbknn(zicdata,batch_key='species',n_pcs=50,neighbors_within_batch=9)
sc.tl.leiden(zicdata,resolution=1.5)#.8 works for cross species when not marker subsetted
sc.tl.umap(zicdata)


# In[107]:


sc.pl.umap(zicdata,color=['noctx_supervised_name','species','region','leiden'])


# In[108]:


sc.pl.umap(zicdata,color=['leiden'],legend_loc='on data')


# In[116]:


sc.pl.umap(zicdata,color=['LHX6','LHX8','NR2F2','CRABP1','LAMP5','PAX6','NR4A2','SCGN','TH','TRH','NR4A2','ZIC1','ZIC4','FOXG1','SLC17A6','SLC17A7'],use_raw=False)


# In[110]:


df_plot = zicdata.obs.groupby(['leiden', 'species']).size().reset_index().pivot(columns='species', index='leiden', values=0).apply(lambda g: g / g.sum(),1)
ax = df_plot.plot(kind='bar', legend=False,stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('proportion of cells in region')
ax.grid(False)


# In[ ]:





# In[115]:


sc.pl.stacked_violin(zicdata[zicdata.obs.leiden.isin(['3','5']) ,:],groupby='species',var_names=['LHX6','LHX8','CRABP1','LAMP5','PAX6','NR2F2','NR4A2','SCGN','TH','TRH','NR4A2','ZIC1','ZIC4','SLC17A6','SLC17A7'],use_raw=False,rotation=90)


# In[96]:


sc.pl.stacked_violin(zicdata[zicdata.obs.leiden.isin(['10']),:],groupby='species',var_names=['LHX6','LHX8','NR2F2','CRABP1','LAMP5','PAX6','NR4A2','SCGN','TH','TRH','NR4A2','ZIC1','ZIC4','SLC17A6','SLC17A7'],use_raw=False,rotation=90)


# In[ ]:





# In[91]:


hierarchy_key='leiden'
rgs=sc.tl.rank_genes_groups(zicdata,groupby=hierarchy_key,method='logreg',use_raw=False,copy=True).uns['rank_genes_groups']#,penalty='elasticnet',solver='saga')#or penalty='l1'
result=rgs
groups = result['names'].dtype.names
df=pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores']})
df.to_csv(os.path.join(sc.settings.figdir,"LogReg"+hierarchy_key+"Norm.csv"))
topgenes=df.iloc[0:4,['_n' in x for x in df.columns]].T.values
cols=df.columns[['_n' in x for x in df.columns]]
cols=[re.sub('_n','',x) for x in cols]
topdict=dict(zip(cols,topgenes))
sc.tl.dendrogram(zicdata,groupby=hierarchy_key)
var_dict=dict(zip(zicdata.uns["dendrogram_['"+hierarchy_key+"']"]['categories_ordered'],[topdict[x] for x in zicdata.uns["dendrogram_['"+hierarchy_key+"']"]['categories_ordered']]))
sc.pl.matrixplot(zicdata,groupby=hierarchy_key,var_names=var_dict,save='top_degenes',cmap='RdBu_r',use_raw=False,dendrogram=True)


# In[ ]:





# In[51]:


genez=['FOXP1','FOXP2','FOXP4','PCP4','CBLN2','LGR4','TSHZ1','SCGN','MEIS2','PAX6','TH','TRH']


# In[28]:


sc.pl.stacked_violin(adata[adata.obs.noctx_supervised_name.str.contains('FOXP2|PGC-2') & adata.obs.species.str.contains('Mouse'),:],groupby='noctx_supervised_name',var_names=genez,use_raw=False,rotation=90,standard_scale='obs')


# In[27]:


sc.pl.stacked_violin(adata[adata.obs.noctx_supervised_name.str.contains('FOXP2') & adata.obs.species.str.contains('Macaque'),:],groupby='noctx_supervised_name',var_names=genez,use_raw=False,rotation=90,standard_scale='obs')


# In[28]:


foxdata=adata[adata.obs.noctx_supervised_name.str.contains('FOXP2|PGC-2|Adult'),:].copy()


# In[36]:


foxdata.X=foxdata.raw.X[:,foxdata.raw.var.index.isin(adata.var.index)]
sc.pp.filter_genes(adata,min_cells=20)
sc.pp.normalize_total(foxdata,exclude_highly_expressed=True)
sc.pp.log1p(foxdata)
sc.pp.highly_variable_genes(foxdata,n_top_genes=1000,batch_key='species')
foxdata.var['highly_variable']=(foxdata.var['highly_variable']&(foxdata.var['highly_variable_nbatches']>1))
sc.pp.scale(foxdata,max_value=10)
sc.pp.pca(foxdata,n_comps=50)


# In[37]:


bbknn.bbknn(foxdata,batch_key='species',n_pcs=50,neighbors_within_batch=9)
sc.tl.leiden(foxdata,resolution=1.5)#.8 works for cross species when not marker subsetted
sc.tl.umap(foxdata)


# In[38]:


sc.pl.umap(foxdata,color='noctx_supervised_name')


# In[39]:


sc.pl.umap(foxdata,color='species')


# In[52]:


genez=['FOXP1','FOXP2','FOXP4','DRD1','DRD2','DRD3','DRD5','PCP4','CBLN2','LGR4','TSHZ1','SCGN','MEIS2','ETV1','PAX6','TH','TRH','SP8','NR2F1','NR2F2','PROX1','LHX6','ZIC1','ZIC2','CRABP1','ID2', 'KLF7', 'SALL3', 'NR4A2']


# In[53]:


sc.pl.umap(foxdata,color=genez,use_raw=False)


# In[44]:


foxdata.obs['species_supervised_name']=foxdata.obs['species'].astype(str)+' '+foxdata.obs['noctx_supervised_name'].astype(str)


# In[55]:


sc.pl.umap(foxdata,color='species_supervised_name')


# In[47]:


sc.pl.dendrogram(foxdata,groupby='species_supervised_name')


# In[64]:


sc.pl.umap(foxdata,color='leiden',legend_loc='on data')


# In[65]:


df_plot = foxdata.obs.groupby(['leiden', 'species']).size().reset_index().pivot(columns='species', index='leiden', values=0).apply(lambda g: g / g.sum(),1)
ax = df_plot.plot(kind='bar', legend=False,stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('proportion of cells in region')
ax.grid(False)


# In[68]:


sc.pl.violin(foxdata[foxdata.obs.leiden.str.contains('22'),:],keys=genez,groupby='species')


# In[70]:


sc_utils.log_reg_diff_exp(foxdata,obs_name='species_supervised_name')


# In[26]:


alldata=sc.read('/wynton/group/ye/mtschmitz/marmosetfastqpool/QvRvM.h5ad')


# In[29]:


foxdata.X=foxdata.raw.X[:,foxdata.raw.var.index.isin(adata.var.index)]


# In[30]:


#foxdata.obs.species.replace({'Macaque':'macaque','Mouse':'mouse'})


# In[31]:


alldata=anndata.AnnData.concatenate(alldata,foxdata)


# In[32]:


alldata.raw=alldata


# In[33]:


alldata=alldata[alldata.obs.species.isin(['Macaque','marmoset']),:]


# In[39]:


alldata.X=alldata.raw.X[:,alldata.raw.var.index.isin(alldata.var.index)]
sc.pp.filter_genes(alldata,min_cells=10)
sc.pp.filter_cells(alldata,min_genes=500)

sc.pp.normalize_total(alldata,exclude_highly_expressed=True)
sc.pp.log1p(alldata)
sc.pp.highly_variable_genes(alldata,n_top_genes=3000,batch_key='species')
alldata.var['highly_variable']=(alldata.var['highly_variable']&(alldata.var['highly_variable_nbatches']>1))
print(alldata.X.shape,flush=True)
sc.pp.scale(alldata,max_value=10)
sc.pp.pca(alldata,n_comps=50)
#alldata.obsm['X_pca']=alldata.obsm['X_pca'][:,1:]
#sc.pp.neighbors(alldata)
scv.pp.remove_duplicate_cells(alldata)
bbknn.bbknn(alldata,batch_key='interspec_batch',n_pcs=50,neighbors_within_batch=6)
sc.tl.leiden(alldata,resolution=2.5)#.8 works for cross species when not marker subsetted
sc.tl.umap(alldata,spread=3,min_dist=.1)


# In[40]:


sc.pl.umap(alldata,color='species')


# In[41]:


sc.pl.umap(alldata,color='leiden',legend_loc='on data')


# In[42]:


genez=['FOXP1','FOXP2','FOXP4','PCP4','CBLN2','LGR4','VAX1','TSHZ1','SCGN','MEIS2','ETV1','PAX6','TH','TRH','SP8','NR2F1','NR2F2','PROX1','LHX6','ZIC1','ZIC2','CRABP1','ID2', 'KLF7', 'SALL3', 'NR4A2']


# In[43]:


sc.pl.umap(alldata,color=genez,use_raw=False)


# In[51]:


sc.pl.umap(alldata,color=['DRD1','DRD2','DRD3','DRD5'],use_raw=False)


# In[44]:


sc.pl.umap(alldata,color='supervised_name',legend_loc='on data')
sc.pl.umap(alldata,color='noctx_supervised_name',legend_loc='on data')


# In[45]:


scanpy.tl.score_genes(alldata, ['FOXP2','FOXP4','VAX1','TSHZ1','TH','MEIS2','SCGN'], ctrl_size=50, gene_pool=None, n_bins=25, score_name='FOXP2_score', use_raw=False)


# In[46]:


scanpy.tl.score_genes(alldata, ['FOXP1','ISL1','PENK','TLE4'], ctrl_size=50, gene_pool=None, n_bins=25, score_name='FOXP1_score', use_raw=False)


# In[47]:


scanpy.tl.score_genes(alldata, ['AIF1','CRABP1','LHX6','LHX8','PDGFRA','MOG','NR2F2','PROX1','AQP4','IL33','SLC17A6','ZIC1','ZIC2'], ctrl_size=50, gene_pool=None, n_bins=25, score_name='Other_score', use_raw=False)


# In[48]:


sc.pl.umap(alldata,color=['FOXP1_score','FOXP2_score','Other_score'])
sc.pl.violin(alldata,groupby='leiden',keys=['FOXP1_score','FOXP2_score','Other_score'])


# In[54]:


df=alldata.obs.loc[:,['leiden','FOXP1_score','FOXP2_score','Other_score']].groupby('leiden').mean()
keep_clusts=df.idxmax(1)
keep_clusts=keep_clusts.index[keep_clusts.str.contains('FOXP2')]


# In[65]:


df


# In[58]:


sc.pl.umap(alldata[alldata.obs.leiden.isin(keep_clusts)|alldata.obs.supervised_name.str.contains('SCGN|FOXP2'),:],color=genez,use_raw=False)


# In[73]:


allfoxdata=alldata[alldata.obs.leiden.isin(list(keep_clusts)+['16','31'])|alldata.obs.supervised_name.str.contains('SCGN|FOXP2'),:].copy()
allfoxdata=allfoxdata[~allfoxdata.obs.supervised_name.str.contains('CRABP'),:]


# In[86]:


allfoxdata.X=allfoxdata.raw.X[:,allfoxdata.raw.var.index.isin(allfoxdata.var.index)].todense()
sc.pp.filter_genes(allfoxdata,min_cells=10)
sc.pp.filter_cells(allfoxdata,min_genes=500)

sc.pp.normalize_total(allfoxdata,exclude_highly_expressed=True)
sc.pp.log1p(allfoxdata)
sc.pp.highly_variable_genes(allfoxdata,n_top_genes=3000,batch_key='species')
allfoxdata.var['highly_variable']=(allfoxdata.var['highly_variable']&(allfoxdata.var['highly_variable_nbatches']>3))
print(allfoxdata.X.shape,flush=True)
sc.pp.scale(allfoxdata,max_value=10)
sc.pp.pca(allfoxdata,n_comps=50)
#allfoxdata.obsm['X_pca']=allfoxdata.obsm['X_pca'][:,1:]
#sc.pp.neighbors(allfoxdata)
scv.pp.remove_duplicate_cells(allfoxdata)
bbknn.bbknn(allfoxdata,batch_key='interspec_batch',n_pcs=50,neighbors_within_batch=6)
sc.tl.leiden(allfoxdata,resolution=2.5)#.8 works for cross species when not marker subsetted
sc.tl.umap(allfoxdata,spread=3,min_dist=.1)


# In[87]:


sc.pl.umap(allfoxdata,color=genez,use_raw=False)


# In[104]:


sc.pl.umap(allfoxdata,color='noctx_supervised_name',use_raw=False)


# In[89]:


sc.pl.umap(allfoxdata,color='species',use_raw=False)


# In[93]:


sc.pl.umap(allfoxdata,color='leiden',legend_loc='on data')


# In[91]:


df_plot = allfoxdata.obs.groupby(['leiden', 'species']).size().reset_index().pivot(columns='species', index='leiden', values=0).apply(lambda g: g / g.sum(),1)
ax = df_plot.plot(kind='bar', legend=False,stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('proportion of cells in region')
ax.grid(False)


# In[101]:


sc.pl.stacked_violin(allfoxdata[allfoxdata.obs.leiden.str.contains('29'),:],groupby='species',var_names=genez,use_raw=False,rotation=90)


# In[102]:


sc.pl.stacked_violin(allfoxdata[allfoxdata.obs.leiden.isin(['20']),:],groupby='species',var_names=genez,use_raw=False,rotation=90)


# In[103]:


sc.pl.stacked_violin(allfoxdata[allfoxdata.obs.leiden.isin(['4']),:],groupby='species',var_names=genez,use_raw=False,rotation=90)


# In[105]:


sc.pl.stacked_violin(allfoxdata[allfoxdata.obs.leiden.isin(['16']),:],groupby='species',var_names=genez,use_raw=False,rotation=90)


# In[ ]:


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

