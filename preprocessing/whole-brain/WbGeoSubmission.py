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
sc.settings.figdir='/wynton/group/ye/mtschmitz/figures/geo/'
scv.settings.figdir='/wynton/group/ye/mtschmitz/figures/geo/'
sc.settings.file_format_figs='pdf'
sc.settings.autosave=False
sc.settings.autoshow=True
import matplotlib.font_manager
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Nimbus Sans','Arial']})
matplotlib.rc('text', usetex=False)
scanpy.set_figure_params(scanpy=True,dpi_save=300)


# In[26]:


adata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityMacaqueGeAllcortexHippocampusProcessed.h5ad')


# In[90]:


get_ipython().run_line_magic('matplotlib', 'inline')
adata


# In[29]:


sc.pl.umap(adata,color='region')
cortex_order=['pfc','cingulate','motor','somato','temporal','insula','hippocampus','parietal','v1']
cortex_colors=seaborn.color_palette("YlOrBr",n_colors=len(cortex_order)+2).as_hex()[2:]
ventral_tel=['ge','cge', 'cge_and_lge', 'lge','mge_and_cge_and_lge', 'mge']
vt_colors=seaborn.blend_palette(('fuchsia','dodgerblue'),n_colors=len(ventral_tel)).as_hex()
med_tel=['septum','pre-optic', 'hypothalamusandpoa','septumandnucleusaccumbens']
mt_colors=seaborn.blend_palette(('grey','black'),n_colors=len(med_tel)).as_hex()
basal_gang=['putamen_and_septum','str', 'putamenandclaustrumandbasalganglia', 'amy']
bg_colors=seaborn.color_palette("Greens",n_colors=len(basal_gang)).as_hex()

all_regions=cortex_order+ventral_tel+med_tel+basal_gang
all_regions_colors=cortex_colors+vt_colors+mt_colors+bg_colors

all_regions_dict={'pfc':'PFC',
 'cingulate':'Cingulate',
 'motor':'Motor',
 'somato':'Somato',
 'temporal':'Temporal',
 'insula':'Insula',
 'hippocampus':'Hippocampus',
 'parietal':'Parietal',
 'v1':'V1',
 'cge':'CGE',
 'cge_and_lge':'CGE+LGE',
 'lge':'LGE',
 'mge_and_cge_and_lge':'MGE+CGE+LGE',
 'mge':'MGE',
 'septum':'Septum',
 'ge':'dLGE',
 'xLGE':'dLGE',
 'pre-optic':'POA',
 'hypothalamusandpoa':'POH+POA',
 'septumandnucleusaccumbens':'Septum+NAc',
 'putamen_and_septum':'Septum+Putamen',
 'str':'Striatum',
 'putamenandclaustrumandbasalganglia':'Striatum+Claustrum',
 'amy':'Amygdala'}
all_regions=[all_regions_dict[x] for x in all_regions]
cortex_order=[all_regions_dict[x] for x in cortex_order]
region_color_dict=dict(zip(all_regions,all_regions_colors))

adata.obs['region']=adata.obs['region'].replace(all_regions_dict)
adata.obs['region']=adata.obs['region'].astype('category')
adata.obs.region.cat.reorder_categories(all_regions,inplace=True,ordered=True)
adata.uns['region_colors']=all_regions_colors

sc.pl.umap(adata,color='region')


# In[30]:


supercell=pd.read_csv('/wynton/group/ye/mtschmitz/macaquedevbrain/MacaqueGEsupervisednamesHippo3PreserveCompoundname.txt')
ind=adata.obs.index[adata.obs.index.isin(supercell['full_cellname'])]
supercell.index=supercell['full_cellname']
supercell=supercell.loc[ind,:]
adata.obs.loc[ind,'supervised_name']=supercell['supervised_name']
adata.obs.loc[ind,'hires_leiden']=supercell['leiden']
adata.obs['supervised_name']=adata.obs['supervised_name'].astype(str)
adata.obs['supervised_name']=[re.sub('Cortical ','',x) for x in adata.obs['supervised_name'] ]
adata.obs['supervised_name']=adata.obs['supervised_name'].astype('category')
adata.obs['noctx_supervised_name']=adata.obs['supervised_name']


# In[35]:


obs=adata.obs


# In[32]:


supercell=pd.read_csv('/wynton/group/ye/mtschmitz/macaquedevbrain/MacaqueGEsupervisednamesHippoPredictedEnd.txt')
ind=adata.obs.index[adata.obs.index.isin(supercell['full_cellname'])]
supercell.index=supercell['full_cellname']
supercell=supercell.loc[ind,:]
adata.obs['latent_time']=np.nan
adata.obs['predicted_end']=''
adata.obs.loc[ind,['noctx_supervised_name','predicted_end','latent_time']]=supercell.loc[:,['noctx_supervised_name','predicted_end','latent_time']]
adata.obs.loc[adata.obs['predicted_end'].isna(),'noctx_supervised_name']=adata.obs['noctx_supervised_name'][adata.obs['predicted_end'].isna()]


# In[33]:


pd.DataFrame(adata.obsm['X_umap']).to_csv(os.path.expanduser('~/MacaqueUMAP.csv'))


# In[38]:


adata.obs['intronic_umi']=adata.layers['unspliced'].sum(1)
adata.obs['exonic_umi']=adata.layers['spliced'].sum(1)
adata.obs['umi_per_cell']=adata.obs['intronic_umi']+adata.obs['exonic_umi']
adata.obs['percent_introns']=adata.obs['intronic_umi']/adata.obs['umi_per_cell']
adata.obs['percent_exons']=adata.obs['exonic_umi']/adata.obs['umi_per_cell']


# In[39]:


adata.obs.loc[:,['batch','timepoint','region','file_name','supervised_name','leiden','n_genes','latent_cell_probability','phase','latent_time','intronic_umi','exonic_umi','umi_per_cell','percent_introns','percent_exons']].to_csv(os.path.expanduser('~/MacaqueCellMetadata.csv'))


# In[28]:


adata.obs=adata.obs.loc[:,['batch','batch_name','file_name','timepoint','region','supervised_name','hires_leiden','leiden','n_genes','latent_cell_probability','phase','latent_time']]


# In[29]:


import json
filenames=[]
jsons=[]
for root, dirs, files in os.walk('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_kallisto', topdown=True):
    dirs[:] = [d for d in dirs if '_Out' not in d]
    for name in files:
        if name =='run_info.json':
            print(root)
            filenames.append(root)
            f=open(os.path.join(root, name))
            jsons.append(json.load(f))
            f.close()

jsonframe=pd.DataFrame(jsons)
jsonframe.index=[x.split('/')[-1] for x in filenames]
jsonframe=jsonframe.loc[jsonframe.index.isin(adata.obs['file_name'].unique()),:]


# In[30]:


jsonframe.drop(columns=['n_targets','n_bootstraps','start_time','call']).to_csv(os.path.expanduser('~/MacaqueSampleMetadata.csv'))


# In[31]:


adata.obs=adata.obs.rename(columns={'supervised_name':'class'})


# In[32]:


adata.var=adata.var.loc[:,['feature_type-0','id-0']]


# In[33]:


adata.uns['class_colors']=adata.uns['noctx_supervised_name_colors']


# In[34]:


adata.uns['regions']=list(adata.obs.region.cat.categories)
adata.uns['classes']=list(adata.obs['class'].cat.categories)


# In[35]:


adata.obs


# In[36]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.pl.umap(adata,color=['class'])


# In[37]:


del adata.layers['ambiguous']
del adata.layers['spliced']
del adata.layers['unspliced']


# In[38]:


adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)]


# In[ ]:


adata.write_csvs('/wynton/home/ye/mschmitz1/INcsvs/',skip_data=False)


# In[39]:


sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.scale(adata,max_value=10)


# In[40]:


adata.write('/wynton/home/ye/mschmitz1/MacaqueDevInhibitoryNeurons.h5ad')


# In[24]:


adata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityDYNAMICALMouseWbGe.h5ad')


# In[31]:


adata.obs['clippedtp']=np.clip(adata.obs['timepoint'],0,24)


# In[32]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.settings.figdir='/wynton/group/ye/mtschmitz/figures/mouseWbGeCB202002/'
scv.settings.figdir='/wynton/group/ye/mtschmitz/figures/mouseWbGeCB202002/'
sc.pl.umap(adata,color=['linear_velocity_length', 'linear_velocity_confidence','latent_time','root_cells','end_points','linear_velocity_pseudotime','clippedtp'],color_map=matplotlib.cm.RdYlBu_r,save='linear_velocity_stats')


# In[22]:


sc.pl.umap(adata,color=['latent_time'],color_map=matplotlib.cm.RdYlBu_r)


# In[8]:


adata.obs=adata.obs.rename(columns={'supervised_name':'class'})


# In[9]:


adata.obs.loc[:,['dataset_name','batch_name','timepoint','region','class','leiden','n_genes','latent_cell_probability','phase','latent_time']].to_csv(os.path.expanduser('~/MouseCellMetadata.csv'))


# In[10]:


adata.obs=adata.obs.loc[:,['batch','batch_name','timepoint','region','class','old_leiden','leiden','n_genes','latent_cell_probability','phase','latent_time']]


# In[11]:


adata.obs.drop('old_leiden',axis=1,inplace=True)


# In[12]:


adata.var=adata.var.loc[:,['feature_type-0-0','id-0-0','name-0-0']]
adata.var.index=list(adata.var['name-0-0'])
adata.raw.var.index=list(adata.raw.var['name-0-0'])


# In[13]:


adata.var


# In[14]:


oldvar=adata.var


# In[15]:


adata.layers


# In[16]:


del adata.layers['spliced']
del adata.layers['unspliced']
del adata.layers['Ms']
del adata.layers['Mu']


# In[17]:


adata.uns['regions']=list(adata.obs.region.cat.categories)
adata.uns['classes']=list(adata.obs['class'].cat.categories)


# In[18]:


adata.uns['class_colors']=adata.uns['supervised_name_colors']


# In[19]:


adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)]


# In[ ]:


adata.write_csvs('/wynton/home/ye/mschmitz1/MouseDevINcsvs/',skip_data=False)


# In[20]:


sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.scale(adata,max_value=10)


# In[ ]:


#adata.write_csvs('/wynton/home/ye/mschmitz1/MouseDevINcsvs/',skip_data=False)


# In[21]:


adata.write('/wynton/home/ye/mschmitz1/MouseDevInhibitoryNeurons.h5ad')


# In[2]:


adata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityDYNAMICALMouseWbGeAdult.h5ad')


# In[64]:


adata=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/KDCbVelocityMouseWbGeAdult.h5ad')


# In[3]:


adata


# In[4]:


adata


# In[5]:


adata.layers['unspliced'].sum(1)


# In[6]:


adata.obs['intronic_umi']=adata.layers['unspliced'].sum(1)
adata.obs['exonic_umi']=adata.layers['spliced'].sum(1)
adata.obs['umi_per_cell']=adata.obs['intronic_umi']+adata.obs['exonic_umi']
adata.obs['percent_introns']=adata.obs['intronic_umi']/adata.obs['umi_per_cell']
adata.obs['percent_exons']=adata.obs['exonic_umi']/adata.obs['umi_per_cell']


# In[7]:


del adata.layers['spliced']
del adata.layers['unspliced']
del adata.layers['Ms']
del adata.layers['Mu']


# In[8]:


adata.obs['agg_supervised_name']='nan'
adata.obs['old_leiden']='nan'
supercell=pd.read_csv('/wynton/group/ye/mtschmitz/macaquedevbrain/MouseAdultAggSupervised.txt')
adata=adata[~adata.obs.index.duplicated(),:]
supercell=supercell.loc[supercell['agg_supervised_name']!='nan',:]
ind=adata.obs.index[adata.obs.index.isin(supercell['full_cellname'])]
supercell.index=supercell['full_cellname']
supercell=supercell.loc[~supercell.index.duplicated(),:]
supercell=supercell.loc[ind,:]
adata.obs.loc[ind,'agg_supervised_name']=supercell['agg_supervised_name']
adata.obs.loc[ind,'old_leiden']=supercell['leiden']
adata.obs['agg_supervised_name']=adata.obs['agg_supervised_name'].astype(str)
adata.obs['old_leiden']=adata.obs['old_leiden'].astype(str)
adata=adata[adata.obs['agg_supervised_name']!='nan',:]


# In[9]:


adata.obs['leiden']=adata.obs['old_leiden']


# In[10]:


adata.obs['region']=[re.sub('srr|samn','',x) for x in adata.obs['region'].astype(str)]
adata.obs['region']=[re.sub('ctx','cortex',x) for x in adata.obs['region'].astype(str)]
adata.obs['region']=[re.sub('ssctx','cortex',x) for x in adata.obs['region'].astype(str)]
adata.obs['region']=adata.obs.region.astype(str)
adata.obs.loc[adata.obs.region.str.contains('ltx'),'region']='cortex'
adata.obs.loc[adata.obs.region.str.contains('(?i)sscortex'),'region']='cortex'
adata.obs.loc[adata.obs.region.str.contains('dentgyr|^ca$'),'region']='hippocampus'

adata.obs['timepoint']=adata.obs['timepoint'].astype(str)
adata.obs.loc[adata.obs['dataset_name']=='PRJNA498989_OB_mouse','timepoint']=84
adata.obs.loc[adata.obs['dataset_name']=='PRJNA515751_konopka_striatum','timepoint']=30
adata.obs.loc[adata.obs['batch_name'].str.contains('P5_',case=False),'timepoint']=26
adata.obs.loc[adata.obs['batch_name'].str.contains('p07_Cortex_SRR11947654',case=False),'timepoint']=28
adata.obs.loc[adata.obs['timepoint']=='nan','timepoint']=84
adata.obs.loc[adata.obs['timepoint'].astype(float)>100,'timepoint']=84

adata.obs.loc[adata.obs['dataset_name'].str.contains('dev_hypo'),'timepoint']=[sc_utils.tp_format_mouse(x) for x in adata.obs.loc[adata.obs['dataset_name'].str.contains('dev_hypo'),'batch_name']]
adata.obs['timepoint']=adata.obs['timepoint'].astype(float)


# In[ ]:


adata=adata[adata.obs['agg_supervised_name']!='nan',:]


# In[ ]:


adata.obs=adata.obs.rename(columns={'agg_supervised_name':'class'})
adata.obs=adata.obs.rename(columns={'old_leiden':'hires_leiden'})


# In[11]:


adata.obs


# In[ ]:


adata.obs=adata.obs.loc[:,['dataset_name','batch_name','timepoint','region','class','leiden','n_genes','latent_cell_probability','intronic_umi','exonic_umi','umi_per_cell','percent_introns','percent_exons']]


# In[ ]:


adata.obs.loc[:,['dataset_name','batch_name','timepoint','region','class','leiden','n_genes','latent_cell_probability','intronic_umi','exonic_umi','umi_per_cell','percent_introns','percent_exons']].to_csv(os.path.expanduser('~/MouseAdultCellMetadata.csv'))


# In[83]:


import json
filenames=[]
jsons=[]
for root, dirs, files in os.walk('/wynton/group/ye/mtschmitz/mousefastqpool/', topdown=False):
    for name in files:
        if name =='run_info.json':
            print(root)
            filenames.append(root)
            f=open(os.path.join(root, name))
            jsons.append(json.load(f))
            f.close()
jsonframe=pd.DataFrame(jsons)
jsonframe.index=[x.split('/')[-1] for x in filenames]
jsonframe=jsonframe.loc[jsonframe.index.isin((adata.obs['batch_name'].astype(str)+'_kOut').unique()),:]


# In[ ]:


jsonframe


# In[ ]:


jsonframe.drop(columns=['n_targets','n_bootstraps','start_time','call']).to_csv(os.path.expanduser('~/MouseSampleMetadata.csv'))


# In[ ]:


#adata.obs.drop('old_leiden',axis=1,inplace=True)


# In[84]:


adata.var=adata.var.loc[:,['feature_type-0-10','id-0-10','name-0-10']]
adata.var.index=list(adata.var['name-0-10'])
adata.raw.var.index=list(adata.raw.var['name-0-10'])


# In[85]:


adata


# In[86]:


get_ipython().run_line_magic('matplotlib', 'inline')
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
    'Str_LHX8/CHAT':'olivedrab',
    'Amy/Hypo_HAP1/PEG10':'darkslateblue',
    'GP_GBX1/GABRA1':'teal',
    'vSTR_HAP1/ZIC1':'darkslategray',
    'Excitatory':'whitesmoke',
    'Ctx_PROX1/LAMP5':'darkgoldenrod', 
    'Ctx_LHX6/LAMP5':'rosybrown', 
    'OB-GC_RPRM':'palegoldenrod', 
    'OB-GC_STXBP6/PENK':'olive', 
    'OB-PGC_FOXP2/CALB1':'aquamarine', 
    'OB-PGC_TH/SCGN':'orange', 
    'OB-PGC_ZIC':'saddlebrown', 
    'Ctx_LHX6/PVALB':'hotpink', 
    'Ctx_PROX1/SNCG':'lightslategray', 
    'Ctx_LHX6/SST':'salmon', 
    'Ctx/BN_SST/CHODL':'maroon',
    'Ctx_CCK/VIP':'tan',
    'Ctx_NR2F2/PAX6':'sienna',
    'Ctx_PVALB/VIPR2':'lightcoral',
    'BN-eSPN_FOXP2/TSHZ1':'springgreen', 
    'Str-dSPN_FOXP1/ISL1':'cyan', 
    'Str-iSPN_FOXP1/PENK':'navy',
    'vStr_DRD1/NPY1R':'violet',
    'Str-IN_CRABP1/MAF':'black',
    'Glia':'lightgoldenrodyellow',
    'OB-GC NR2F2/PENK':'brown',
    'Ctx_LAMP5/NDNF':'khaki',
    'Ctx_SST/NDNF':'peachpuff',
    'Ctx_CCK/DPY19L1':'dimgray'}

sc.pl.umap(adata,color=['class'])
my_pal=[paldict[x] for x in adata.obs['class'].cat.categories]
sc.pl.umap(adata,color=['class'],legend_loc='on data',legend_fontsize=4,palette=my_pal,save='agg_supervised_adult')


# In[87]:


adata.uns['regions']=list(adata.obs.region.cat.categories)
adata.uns['classes']=list(adata.obs['class'].cat.categories)


# In[88]:


#adata.uns['class_colors']=adata.uns['agg_supervised_name_colors']


# In[89]:


adata.X=adata.raw.X[:,adata.raw.var.index.isin(adata.var.index)]


# In[ ]:


adata.write_csvs('/wynton/home/ye/mschmitz1/MouseAdultINcsvs/',skip_data=False)


# In[90]:


sc.pp.normalize_per_cell(adata)
sc.pp.log1p(adata)
sc.pp.scale(adata,max_value=10)


# In[91]:


adata.write('/wynton/home/ye/mschmitz1/MouseAdultInhibitoryNeurons.h5ad')


# In[ ]:





# In[28]:





# In[17]:


region_color_dict=dict(zip(adata.obs.region.cat.categories, adata.uns['region_colors']))


# In[ ]:





# In[24]:


get_ipython().run_line_magic('matplotlib', 'inline')
sc.settings.figdir='/wynton/group/ye/mtschmitz/figures/mouseWbGeAdult2/'
scv.settings.figdir='/wynton/group/ye/mtschmitz/figures/mouseWbGeAdult2/'

adata.obs['agg_supervised_name']=adata.obs['class'].astype('category')
sc.tl.dendrogram(adata,groupby='agg_supervised_name')
sc.pl.dendrogram(adata,groupby='agg_supervised_name',save='agg_supervised_name')
df_plot = adata.obs.groupby(['agg_supervised_name', 'region']).size().reset_index().pivot(columns='region', index='agg_supervised_name', values=0).apply(lambda g: g / g.sum(),1)
df_plot=df_plot.loc[adata.uns["dendrogram_agg_supervised_name"]['categories_ordered'],:]
these_colors=[region_color_dict[x] for x in df_plot.columns]
ax = df_plot.plot(kind='bar', legend=False,stacked=True,color=these_colors,fontsize=7.5)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('proportion of cells in region')
ax.grid(False)
ax.get_figure().savefig(os.path.join(sc.settings.figdir,'supervisedNameRegionsStackBar.pdf'), bbox_inches="tight")


# In[ ]:





# In[ ]:





# In[ ]:


filepath='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_kallisto'
for f,b in zip(obs.file_name.unique(),obs.batch_name.unique()):
    print(f)
    outdir=re.sub('kOut','Out',os.path.join(filepath,f,'outs'))
    cellfile=re.sub('kOut','Out',os.path.join(outdir,'InhibitoryCells.txt'))
    try:
        pd.DataFrame(index=[re.split('_',x)[0]+'-1' for x in obs.loc[obs.file_name.isin([f]),:].index]).to_csv(cellfile,header=False)
        stream = os.popen('/wynton/home/ye/mschmitz1/utils/subset-bam_linux --out-bam '+os.path.join(outdir,b+'_inhibitory.bam')+' --bam '+ os.path.join(outdir,'possorted_genome_bam.bam')+' --cell-barcodes '+ cellfile)
        output = stream.read()
        print(output)
        print('first success')
        print('/wynton/home/ye/mschmitz1/utils/bamtofastq_linux ' + os.path.join(outdir,b+'_inhibitory.bam') + ' '+outdir+'/fastq/' + ' --reads-per-fastq=999999999999')
        stream = os.popen('/wynton/home/ye/mschmitz1/utils/bamtofastq_linux ' + os.path.join(outdir,b+'_inhibitory.bam') + ' '+outdir+'/fastq/' + ' --reads-per-fastq=999999999999')
        output = stream.read()
        print(output)
        print('double success')
        [os.rename(os.path.join(outdir+'/fastq/',y,x),os.path.join(outdir+'/fastq/',y,re.sub('bamtofastq',b,x))) for y in os.listdir(outdir+'/fastq/') for x in os.listdir(outdir+'/fastq/'+y)]
            except:
        pass


# In[90]:


filepath='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_kallisto'
for f,b in zip(obs.file_name.unique(),obs.batch_name.unique()):
    print(f)
    outdir=re.sub('kOut','Out',os.path.join(filepath,f,'outs'))
    cellfile=re.sub('kOut','Out',os.path.join(outdir,'InhibitoryCells.txt'))
    [os.rename(os.path.join(outdir+'/fastq/',y,x),os.path.join(outdir+'/fastq/',y,re.sub('bamtofastq',b,x))) for y in os.listdir(outdir+'/fastq/') for x in os.listdir(outdir+'/fastq/'+y)]


# In[ ]:


'Mac2_V1' in adata.obs.batch_name


# In[ ]:


list(adata.obs.batch_name.unique())


# In[ ]:


adata.obs


# In[ ]:


filedf=pd.DataFrame(list(adata.obs.groupby(['batch_name','file_name','region','timepoint']).sum().dropna().index.to_numpy()))


# In[ ]:


filedf.index=filedf[1]


# In[ ]:


filedf


# In[ ]:


l=[]
filepath='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_kallisto'
for x in os.listdir(filepath):
    if os.path.isdir(os.path.join(filepath,x)) & ('_Out' in x):
        try:
            for y in os.listdir(os.path.join(filepath,x,'outs','fastq')):
                fc=re.split('_',y)
                fc=fc[-1]
                print(fc)
                for z in os.listdir(os.path.join(filepath,x,'outs','fastq',y)):
                    print(os.path.join(filepath,x,'outs','fastq',y,z))
                    print(os.path.join(filepath,x,'outs','fastq',y,fc+'_'+z))
                    #os.rename(os.path.join(filepath,x,'outs','fastq',y,z),os.path.join(filepath,x,'outs','fastq',y,fc+'_'+z))
                    #os.rename(os.path.join(filepath,x,'outs','fastq',y,z),'/wynton/group/ye/mtschmitz/macaquedevbrain/INfastqs/'+z)
        except:
            print('FAIL',os.path.join(filepath,x,'outs','fastq'))


# In[ ]:


l=[]
filepath='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_kallisto'
for x in filedf[1]:
    origx=x
    x=re.sub('kOut','Out',x)
    if 'PEC' in x:
        continue
    for y in os.listdir(os.path.join(filepath,x,'outs','fastq')):
        print(os.path.join(filepath,x,'outs','fastq',y))
        print(os.listdir(os.path.join(filepath,x,'outs','fastq',y)))
        for z in os.listdir(os.path.join(filepath,x,'outs','fastq',y)):
            print(os.path.join(filepath,x,'outs','fastq',y,z))
            l.append([x,os.path.join(filepath,x,'outs','fastq',y,z)]+list(filedf.loc[origx,:]))


# In[ ]:


df=pd.DataFrame(l)
df


# In[182]:


df[1]=[re.split('/',x)[-1] for x in df[1]]


# In[2]:


df.columns


# In[ ]:


#unrelated mouse dataset taken from bams
filepath='/wynton/group/ye/mtschmitz/mousefastqpool/PRJNA637987_lamanno_devmouse'
#To run with kalliso use read schema 1,0,0:2,0,0:0,0,0
for x in os.listdir(filepath):
    try:
        if os.path.isdir(os.path.join(filepath,x)):
            for y in os.listdir(os.path.join(filepath,x)):
                print(os.path.join(filepath,x,y))
                l=re.split('_',y)
                for z in os.listdir(os.path.join(filepath,x,y)):
                    print(os.path.join(filepath,x,y,z))
                    print(os.path.join(filepath,re.sub('bamtofastq',x+'__'+l[-1],z)))
                    os.rename(os.path.join(filepath,x,y,z),os.path.join(filepath,re.sub('bamtofastq',x+'__'+l[-1],z)))
    except:
        pass
    
"""
filepath='/wynton/group/ye/mtschmitz/mousefastqpool/PRJNA637987_lamanno_devmouse'
#To run with kalliso use read schema 1,0,0:2,0,0:0,0,0
for x in os.listdir(filepath):
    if 'fastq.gz' in x:
        #print(x)
        #print(re.split('_',re.sub('.fastq.gz','',x)))
        l=re.split('_',re.sub('.fastq.gz','',x))
        print(os.path.join(filepath,x),os.path.join(filepath,re.sub('L[0-9]+','L'+l[-1],x)))
        os.rename(os.path.join(filepath,x),os.path.join(filepath,re.sub('L[0-9]+','L'+l[-1],x)))
"""


# In[2]:


f = open('/wynton/group/ye/mtschmitz/refdata2/mm10/Mus_musculus.GRCm38.dna.primary_assembly.fa', "r")
l=f.readlines()



# In[3]:


import tqdm
chrstring=[]
for x in tqdm.tqdm(l):
    if x[0] != '>':
        #if x != 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n':
        chrstring[-1].append(x[:-1])
    else:
        chrstring.append([])


# In[4]:


chrstringcat=[]
for x in chrstring:
    chrstringcat.append(''.join(x))


# In[5]:


import re
broken=[]
for x in tqdm.tqdm(chrstringcat):
    broken.append(re.split('GAGCTC',x))


# In[6]:


len(chrstring[0])


# In[7]:


flatlist=[len(item) for sublist in tqdm.tqdm(broken) for item in sublist]    


# In[8]:


print(len(flatlist))


# In[9]:


flatlist[2]


# In[11]:


get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn
seaborn.distplot(np.log10(np.array(flatlist)+1),kde=False)
plt.xlabel('log10 bases')
plt.ylabel('# fragments')


# In[28]:


logvals=np.log10(np.array(flatlist)+1)
h=np.histogram(logvals,bins=int((np.max(logvals)-np.min(logvals))/.1)+1)
10**h[1]


# In[36]:


h=np.histogram(flatlist,bins=[0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,11000,1000000,10000000])
pd.DataFrame(h[::-1]).T

