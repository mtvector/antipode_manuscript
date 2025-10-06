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
#Load my pipeline functions
import importlib
import importlib.util
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../utils/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)
figdir='/wynton/group/ye/mtschmitz/WBfigures/macaqueMultiseq_fixed/'
sc.settings.figdir=figdir
scv.settings.figdir=figdir
sc.settings.file_format_figs='pdf'
sc.settings.autosave=True
sc.settings.autoshow=False
import h5py
print(h5py.__version__)

def most_frequent(List): 
    return max(set(List), key = List.count)

def variance(X,axis=1):
    return( (X.power(2)).mean(axis=axis)-(np.power(X.mean(axis=axis),2)) )

filepath='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT_chang_kOut'
files=os.listdir(filepath)
files=[f for f in files if 'kOut' in f]
regionlist=['multi']#,'re-optic','poa']
fileList=[f for f in files if any([x in f.lower() for x in regionlist]) ]

newfile='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT_chang_h5ad/KDCbVelocityMultiseq.h5ad'
min_genes=700

if True:#not os.path.exists(newfile):
    adatas=[]
    filename='aem_cellbended_150_750_175e_V0.2'
    for f in fileList:
        print(f,flush=True)
        try:
            if os.path.exists(os.path.join(filepath,f,filename,'aem_cellbended_filtered.h5')):
                sadata = sc_utils.readCellbenderH5(os.path.join(filepath,f,filename,'aem_cellbended_filtered.h5'))
                sadata =sadata[sadata.obs['latent_cell_probability']>.99,:]
                sc.pp.filter_cells(sadata,min_genes=min_genes)
            else:
                sadata=sc_utils.loadPlainKallisto(os.path.join(filepath,f,'all_em'),min_genes=min_genes)
                sc.pp.filter_cells(sadata,min_genes=min_genes)
            sadata.obs.index=[re.sub("-1","",x) for x in sadata.obs.index]
            multiseq=pd.read_csv(os.path.join(filepath,f,'multiseq_calls.csv'),index_col=0)
            print(multiseq)
            #multiseq=multiseq.loc[multiseq['demux_type']=='singlet',:]
            print(multiseq,flush=True)
            print(sadata.obs.index,flush=True)
            sadata=sadata[sadata.obs.index.isin(multiseq.index),:]
            print(sadata,flush=True)
            multiseq=multiseq.loc[multiseq.index.isin(sadata.obs.index),:]
            multiseq=multiseq.loc[sadata.obs.index,:]
            print(multiseq,flush=True)
            sadata.obs['demux_type']=None
            sadata.obs['assignment']=None
            sadata.obs.loc[:,['demux_type','assignment']]=multiseq.loc[:,['demux_type','assignment']]
            sadata.obs.loc[sadata.obs['demux_type']=='singlet','region']=sadata.obs['assignment'][sadata.obs['demux_type']=='singlet']
            #sadata=sadata[sadata.obs['demux_type']!='doublet',:]
            sadata.obs['region']=[sc_utils.region_format_macaque(x) for x in sadata.obs['region']]
            sadata.obs['region']=[x.lower() for x in sadata.obs['region']]
            print(sadata.obs['region'])
            sadata.uns['name']=f
            sadata.obs['batch_name']=str(sadata.uns['name'])
            sadata.obs['timepoint']=sc_utils.tp_format_macaque(sadata.uns['name'])
            if not os.path.exists(os.path.join(filepath,f,'cbdoublets.txt')):
                doublets=sc_utils.doublescrub(sadata)
                print(doublets,flush=True)
                pd.DataFrame(doublets).to_csv(os.path.join(filepath,f,'cbdoublets.txt'),index=False, header=False)
            try:
                ddf=pd.read_csv(os.path.join(filepath,f,'cbdoublets.txt'),index_col=False, header=None)
                doublets=list(ddf[0])
            except:
                doublets=[]
            sadata=sadata[~sadata.obs.index.isin(doublets),:]
            pd.DataFrame(sadata.obs.index).to_csv(os.path.join(filepath,f,'cellbendedcells.txt'),index=False, header=False)
            if sadata.shape[0]>10:
                adatas.append(sadata)
        except Exception as e:
            print(e)
            print('fail')

    #adatas=adatas+[multi]
    adata=sc.AnnData.concatenate(*adatas)
    adata.var.columns = adata.var.columns.astype(str)
    adata.obs.columns = adata.obs.columns.astype(str)
    adata.obs['clean_cellname']=[re.sub('-[0-9]+','',x) for x in  adata.obs.index]
    adata.obs['full_cellname']=adata.obs['clean_cellname'].astype(str)+'_'+adata.obs['batch_name'].astype(str)
    adata.obs.index=list(adata.obs['full_cellname'])
    adata=sc_utils.sanitize_types_for_h5(adata)
    adata.obs['latent_cell_probability']=adata.obs['latent_cell_probability'].astype(float)
    adata.obs['latent_RT_efficiency']=adata.obs['latent_RT_efficiency'].astype(float)    
    print(adata.obs.dtypes)
    print(adata.var.dtypes)
    adata.raw=adata
    print(adata.obs['latent_cell_probability'])
    print(adata.obs['latent_RT_efficiency'])
    adata.write(newfile)        
    
adata=sc.read(newfile)
print("1",adata.obs)
supercell=pd.read_csv('/wynton/group/ye/mtschmitz/macaquedevbrain/MacaqueAllSupervisedCellNames.txt')
print(supercell)
ind=adata.obs.index[adata.obs['full_cellname'].isin(supercell['full_cellname'])]
supercell.index=supercell['full_cellname']
supercell=supercell.loc[ind,:]
adata.obs.loc[ind,'supervised_name']=supercell['supervised_name']
adata.obs['supervised_name']=adata.obs['supervised_name'].astype(str)
print("2",adata,flush=True)

adata.raw=adata
adata.var_names_make_unique()
ribo_genes=[name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL') ]
adata.obs['percent_ribo'] = np.sum(
adata[:, ribo_genes].X, axis=1) / np.sum(adata.X, axis=1)
mito_genes = [name for name in adata.var_names if name in ['ND1','ND2','ND4L','ND4','ND5','ND6','ATP6','ATP8','CYTB','COX1','COX2','COX3'] or 'chrM-' in name or 'MT-' in name]
adata.obs['percent_mito'] = np.sum(
adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
print(adata,flush=True)
#sc.pl.violin(adata,groupby='batch_name',keys=['percent_ribo','percent_mito'],rotation=45,save='ribomito')

adata=adata[adata.obs['percent_ribo'].to_numpy()<.4,:]
adata=adata[adata.obs['percent_mito'].to_numpy()<.15,:]
adata=adata[:,~adata.var.index.isin(mito_genes)]
print(adata)
sc.pp.filter_genes(adata,min_cells=10)
#sc.pp.filter_cells(adata,min_counts=1000)
sc.pp.filter_cells(adata,min_genes=min_genes)
sc.pl.highest_expr_genes(adata, n_top=20, )
print(adata)
sufficient_cells=adata.obs['batch_name'].value_counts().index[adata.obs['batch_name'].value_counts()>50]
adata=adata[adata.obs['batch_name'].isin(sufficient_cells),:]
adata.obs['batch_name'].cat.remove_unused_categories(inplace=True)

sc.pp.normalize_total(adata,exclude_highly_expressed=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=8000,batch_key='batch_name',subset=False)
sc.pp.scale(adata,max_value=10)
sc.pp.pca(adata,n_comps=100)
#sc.pp.neighbors(adata)
bbknn.bbknn(adata,batch_key='batch_name',n_pcs=100,neighbors_within_batch=3)
sc.tl.leiden(adata,resolution=4)
sc.tl.umap(adata,spread=2)

#adata=sc.AnnData(adata.raw.X,var=adata.raw.var,obsm=adata.obsm,uns=adata.uns,obs=adata.obs)
#adata.raw=adata
#sc.pp.normalize_per_cell(adata)
#sc.pp.log1p(adata)
#sc.pp.scale(adata,max_value=10)
sc.pl.umap(adata,color=['leiden'],save='leiden')
sc.pl.umap(adata,color=['batch_name'],save='batch_name')
sc.pl.umap(adata,color=['region'],save='region')
sc.pl.umap(adata,color=['timepoint'],save='timepoint')
easygenes= ['MKI67','SOX2','AQP4','EDNRB','IL33','PDGFRA','HOPX','ERBB4','CALB1','CALB2','GAD1','GAD2','GADD45G','LHX6','SST','NKX2-1','MAF','SP8','SP9','PROX1','VIP','CCK','NPY','LAMP5','HTR3A','NR2F1','NR2F2','TOX3','ETV1','SCGN','FOXG1','FOXP1','FOXP2','FOXP4','TH','DDC','SLC18A2','PAX6','MEIS2','SHTN1','ISL1','PENK','CRABP1','ZIC1','ZIC2','EBF1','DLX5','GSX2','RBFOX3','PPP1R17','EOMES','DCX','TUBB3','BCL11B','TLE4','FEZF2','SATB2','TBR1','AIF1','RELN','PECAM1','HBZ','ROBO1','ROBO2','ROBO3','ROBO4']
easygenes=[x for x in easygenes if x in adata.var.index]
sc.pl.umap(adata,color=easygenes,use_raw=False,save='FaveGenes')

sc.pl.umap(adata,color=['supervised_name'],save='supervised_name')

f = plt.figure()
df_plot = adata.obs.groupby(['supervised_name', 'batch_name']).size().reset_index().pivot(columns='batch_name', index='supervised_name', values=0).apply(lambda g: g / g.sum(),1)
ax = df_plot.plot(kind='bar', legend=False,stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('proportion')
ax.get_figure().savefig(os.path.join(sc.settings.figdir,'supervisedVbatchBar.pdf'), bbox_inches="tight")

f = plt.figure()
df_plot = adata.obs.groupby(['supervised_name', 'leiden']).size().reset_index().pivot(columns='leiden', index='supervised_name', values=0).apply(lambda g: g / g.sum(),1)
ax = df_plot.plot(kind='bar', legend=False,stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('proportion')
ax.get_figure().savefig(os.path.join(sc.settings.figdir,'supervisedVleidenBar.pdf'), bbox_inches="tight")


f = plt.figure()
df_plot = adata.obs.groupby(['supervised_name', 'timepoint']).size().reset_index().pivot(columns='timepoint', index='supervised_name', values=0).apply(lambda g: g / g.sum(),1)
ax = df_plot.plot(kind='bar', legend=False,stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('proportion')
ax.get_figure().savefig(os.path.join(sc.settings.figdir,'supervisedVtimepointBar.pdf'), bbox_inches="tight")


sc_utils.cell_cycle_score(adata)
sc_utils.marker_analysis(adata,variables=['leiden'],markerpath=os.path.expanduser('~/markers.txt')) 
sc_utils.log_reg_diff_exp(adata,obs_name='leiden')

gdata=adata.copy().T

import numba
@numba.njit()
def abscorrelation(x, y):
    mu_x = 0.0
    mu_y = 0.0
    norm_x = 0.0
    norm_y = 0.0
    dot_product = 0.0

    for i in range(x.shape[0]):
        mu_x += x[i]
        mu_y += y[i]

    mu_x /= x.shape[0]
    mu_y /= x.shape[0]

    for i in range(x.shape[0]):
        shifted_x = x[i] - mu_x
        shifted_y = y[i] - mu_y
        norm_x += shifted_x ** 2
        norm_y += shifted_y ** 2
        dot_product += shifted_x * shifted_y

    if norm_x == 0.0 and norm_y == 0.0:
        return 0.0
    elif dot_product == 0.0:
        return 1.0
    else:
        return 1.0 - np.absolute(dot_product / np.sqrt(norm_x * norm_y))
    
gdata.obsm['X_pca']=gdata.obsm['PCs']

#sc.pp.neighbors(gdata,n_neighbors=25)
sc.pp.neighbors(gdata,metric=abscorrelation,n_neighbors=25)

sc.tl.leiden(gdata,resolution=5)
sc.tl.umap(gdata,spread=5,min_dist=.2)

sc.pl.umap(gdata,color='leiden',legend_loc='on data',palette="Set2",save='GeneLeiden')
print(gdata[s_genes+g2m_genes,:].obs.leiden.value_counts())
vcs=gdata[s_genes+g2m_genes,:].obs.leiden.value_counts()
cc_mods=vcs.index[vcs>5]
cc_genes=gdata.obs.index[gdata.obs.leiden.isin(cc_mods)]
adata.var['leiden']=gdata.obs['leiden']
del gdata

adata=adata[(~adata.obs.region.isin(['ob','OB']) &  ~adata.obs.supervised_name.str.contains('"')) | adata.obs.supervised_name.str.contains('MEIS2'),:]

bbknn.bbknn(adata,batch_key='batch_name',n_pcs=100,neighbors_within_batch=3)

sc.pp.highly_variable_genes(adata,flavor='seurat_v3',layer='unspliced',n_top_genes=15000)
adata.var['highly_variable_rank_u']=adata.var['highly_variable_rank']
sc.pp.highly_variable_genes(adata,flavor='seurat_v3',layer='spliced',n_top_genes=15000)
adata.var['highly_variable_rank_s']=adata.var['highly_variable_rank']
print(adata,flush=True)
adata.var['velocity_genes']=adata.var.loc[:,['highly_variable_rank_s','highly_variable_rank_u']].mean(1,skipna=False).rank()<4000
print(adata.var['velocity_genes'].value_counts())

adata.uns['iroot'] = np.flatnonzero(adata.obs.index==sc_utils.get_median_cell(adata,'supervised_name','Transition'))[0]
ncells=20
mgenes=(adata.layers['unspliced']>0).sum(0).A1 < ncells
print(mgenes)
adata.var.loc[mgenes,'velocity_genes']=False

#mgenes=(adata.layers['unspliced']>0).sum(0) > ncells
#adata=adata[:,mgenes.A1]
#adata=adata[:,adata.var['highly_variable']]
#adata.var['velocity_genes']=mgenes.A1&adata.var['highly_variable']
print(list(adata.var['velocity_genes']))
scv.pp.normalize_per_cell(adata)
print(adata)
print(adata.layers['spliced'].shape)
#scv.pp.filter_genes_dispersion(adata)#,n_top_genes=min(6000,adata.layers['spliced'].shape[1]))
scv.pp.log1p(adata)
scv.pp.remove_duplicate_cells(adata)
scv.pp.moments(adata, n_neighbors=20)
print(adata,flush=True)
scv.tl.recover_dynamics(adata,var_names='velocity_genes',use_raw=False)
print('Recovered', flush=True)
scv.tl.velocity(adata,mode='dynamical',filter_genes=False)
print('Velocity done',flush=True)
adata.layers['linear_velocity']=adata.layers['velocity'].copy()
adata.layers['linear_velocity'][:,adata.var.index.isin(cc_genes)]=0
adata.layers['cc_velocity']=adata.layers['velocity'].copy()
adata.layers['cc_velocity'][:,~adata.var.index.isin(cc_genes)]=0

scv.tl.velocity_graph(adata,mode_neighbors='connectivities',vkey='cc_velocity',approx=False)
scv.tl.transition_matrix(adata,vkey='cc_velocity')
scv.tl.terminal_states(adata,vkey='cc_velocity')
scv.tl.recover_latent_time(adata,vkey='cc_velocity')
scv.tl.velocity_confidence(adata,vkey='cc_velocity')
scv.tl.velocity_embedding(adata, basis='umap',vkey='cc_velocity')
adata.obs['cc_latent_time']=adata.obs['latent_time']
scv.pl.velocity_embedding_grid(adata, basis='umap',vkey='cc_velocity',save='_cc_grid',color_map=matplotlib.cm.RdYlBu_r)
try:
    sc.pl.umap(adata,color=['cc_velocity_length', 'cc_velocity_confidence','latent_time','root_cells','end_points','cc_velocity_pseudotime'],save='cc_velocity_stats')
except:
    print('fail cc')

scv.tl.velocity_graph(adata,mode_neighbors='connectivities',vkey='linear_velocity',approx=False)
scv.tl.transition_matrix(adata,vkey='linear_velocity')
scv.tl.terminal_states(adata,vkey='linear_velocity')
scv.tl.recover_latent_time(adata,vkey='linear_velocity')
scv.tl.velocity_embedding(adata, basis='umap',vkey='linear_velocity')
scv.tl.velocity_confidence(adata,vkey='linear_velocity')
scv.pl.velocity_embedding_grid(adata,vkey='linear_velocity',basis='umap',color='latent_time',save='_linear_grid',color_map=matplotlib.cm.RdYlBu_r)
try:
    sc.pl.umap(adata,color=['linear_velocity_length', 'linear_velocity_confidence','latent_time','root_cells','end_points','linear_velocity_pseudotime'],save='linear_velocity_stats')
except:
    'fail linear'
scv.utils.cleanup(adata)
adata.write(re.sub('Velocity','VelocityDYNAMICAL',newfile))

sc.tl.paga(adata,groups='supervised_name')
sc.tl.paga(adata,use_rna_velocity=True,groups='supervised_name')
sc.pl.paga_compare(adata,legend_fontsize=4,arrowsize=10,edge_width_scale=.4,threshold=np.quantile(adata.uns['paga']['connectivities'].data,.9))
sc.pl.paga_compare(adata,legend_fontsize=4,arrowsize=10,edge_width_scale=.4,threshold=np.quantile(adata.uns['paga']['connectivities'].data,.9),save='connectivity')
sc.pl.paga_compare(adata,solid_edges='connectivities',transitions='transitions_confidence',legend_fontsize=4,arrowsize=10,threshold=np.quantile(adata.uns['paga']['transitions_confidence'].data,.9),save='DYNAMICALvelocity')
adata.write(re.sub('Velocity','VelocityDYNAMICAL',newfile))
