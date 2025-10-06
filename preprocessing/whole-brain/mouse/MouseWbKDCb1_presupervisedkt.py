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
import tqdm
import bbknn
#Load my pipeline functions
import importlib
import importlib.util
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../utils/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)
sc.settings.figdir='/wynton/group/ye/mtschmitz/WBfigures/WbMousePresupervisionCB202002/'
scv.settings.figdir='/wynton/group/ye/mtschmitz/WBfigures/WbMousePresupervisionCB202002/'
sc.settings.file_format_figs='pdf'
sc.settings.autosave=True
sc.settings.autoshow=False
print(dir(sc_utils))

pth='/wynton/group/ye/mtschmitz/mousefastqpool'
dirs=['GSE123335_cortex','10X_demonstration','PRJNA411878_10x_mouse','PRJNA637987_lamanno_devmouse','PRJNA547712_dev_hypothalamus']

newfile='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT_chang_h5ad/KDCbVelocityMouseWbPresupervisionKT.h5ad'

min_genes=300

if True:#not os.path.exists(newfile):
    adatas=[]
    filename='t_cellbended_150_750_175e_V0.2'
    for d in dirs:
        print(d,flush=True)
        filepath=os.path.join(pth,d)
        files=os.listdir(filepath)
        files=[f for f in files if 'ktOut' in f]
        print(files)
        for f in files:
            print(f,flush=True)
            if re.search('p10|p14|p45', f.lower()):
                continue
            try:
                if os.path.exists(os.path.join(filepath,f,filename,'t_cellbended_filtered.h5')):
                    sadata = sc_utils.readCellbenderH5(os.path.join(filepath,f,filename,'t_cellbended_filtered.h5'))
                    sadata =sadata[sadata.obs['latent_cell_probability']>.99,:]
                    sc.pp.filter_cells(sadata,min_genes=min_genes)
                else:
                    sadata=sc_utils.loadPlainKallisto(os.path.join(filepath,f,'t_em'),min_genes=min_genes+200)
                import json
                with open(os.path.join(filepath,f,'run_info.json')) as file:
                    data = json.load(file)
                
                sadata.obs["p_pseudoaligned"]=data["p_pseudoaligned"]
                sadata.obs["p_unique"]=data["p_unique"]

                sadata.obs.index=[re.sub("-1","",x) for x in sadata.obs.index]
                sadata.uns['name']=f
                sadata.uns['name']=sc_utils.mouse_process_irregular_names(sadata.uns['name'])
                sadata.obs['batch_name']=str(sadata.uns['name'])
                sadata.obs['dataset_name']=d
                sadata.obs['timepoint']=sc_utils.tp_format_mouse(sadata.uns['name'])
                regionstring=sc_utils.region_format_mouse(sadata.uns['name'])
                regionstring=regionstring.lower()
                sadata.obs['region']=regionstring
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
    print(adata.obs.index)
    print(adata.obs['clean_cellname'])
    #########REGION HANDLING
    from collections import OrderedDict
    reg_dict=OrderedDict({'cerebel|cbc|[^a-z]ce|vermis|cb':'Cb','vz|cortex|ctx|wt[0-9]+|Neurons_Sample_':'Ctx','ob|GSM3449':'OB','head|all|nan|vesicle|placode':'Head','subpall|sp':'FB','thal|pulv|lgn':'Thal',
     'stria|stiatum|putanum|putamen|caud|accumb|nac|basal|SAMN107673':'Str','gp':'GP','clau':'Ctx-Clau','amy':'Amy',
     'ge':'GE','lge':'LGE','mge':'MGE','somato|som':'Ctx-Somato','cge':'CGE','motor|mop|m1|dfc':'Ctx-MOp',
     'pfc|prefront':'Ctx-PFC','cing':'Ctx-Cing','v1|occ':'Ctx-Occ','central|par':'Ctx-Par','temp':'Ctx-Temp',
     'hippo|ca[0-4]|dent|hc':'Ctx-Hip','medull|myelenceph':'Medulla','choroid|cp':'Choroid',
     'sept':'Sept','pre-opt|poa|preoptic':'POA','pons':'Pons','hypo':'Hypo','drg':'DRG','org|cultured':'Cultured','insula':'Ctx-Ins',
     'forebrain|prosen|telenceph':'FB','diencepha':'DE','mesenceph|midbrain|md|[^a-z]mb':'MB','hind|rhombence|metenceph':'HB'
    })

    adata.obs['batch_name']=list(adata.obs['batch_name'].astype(str))
    adata.obs['region']=list(adata.obs['region'].astype(str))
    if 'msregion' in adata.obs.columns:
        adata.obs['msregion']=list(adata.obs['msregion'].astype(str))

    msd=adata.obs['batch_name'].str.contains('multi',case=False)
    adata.obs.loc[msd,'msregion']=adata.obs.loc[msd,'region']
    adata.obs['msregion']=adata.obs['msregion'].astype(str)
    for k in reg_dict.keys():
        print(k)
        adata.obs.loc[adata.obs['msregion'].str.contains(k,case=False),'region']=reg_dict[k]

    for k in reg_dict.keys():
        print(k)
        print(reg_dict[k])
        print(adata.obs.loc[adata.obs['batch_name'].str.contains(k,case=False),'batch_name'])
        adata.obs.loc[adata.obs['batch_name'].str.contains(k,case=False),'region']=reg_dict[k]

    region_groupings={'ctx': ['Ctx-PFC','Ctx-Cing','Ctx-Clau','Ctx-Ins','Ctx-Hip','Ctx-MOp', 'Ctx-Somato','Ctx','Ctx-Temp','Ctx-Par','Ctx-Occ'],
     'bn':[ 'POA','GP','Sept','Str','Amy'],
     'ob':['OB'],
     'cp':['Choroid'],
     'mb':['MB'],
     'hb':['Pons','Cb','HB','DRG','Medulla'],
     'ge':['MGE','GE','LGE','CGE'],
     'de':['Hypo','Thal','DE'],
     'h':['Head'],
     'xvctx':['ultured']}

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

    adata.obs['region']=adata.obs['region'].astype('category')
    regions=[]
    region_colors=[]
    for k in palettes.keys():
        regions+=region_groupings[k]
        region_colors+=palettes[k]

    region_colors=[col for reg,col in zip(regions,region_colors) if reg in adata.obs['region'].cat.categories]
    cats=[x for x in regions if x in adata.obs['region'].cat.categories]
    missings=list(set(adata.obs['region'].cat.categories)-set(cats))

    print(cats+missings)
    paldict=dict(zip(cats+missings,region_colors+['grey' for i in range(len(missings))]))

    region_colors=[paldict[x] for x in adata.obs['region'].cat.categories]
    #adata.obs['region'].cat.categories=cats+missings
    adata.uns['region_cats']=list(adata.obs['region'].cat.categories)
    adata.uns['region_colors']=region_colors
    
    general_region=[]
    for r in adata.obs['region']:
        for k in region_groupings.keys():
            if any([r in x for x in region_groupings[k]]):
                r=k
        general_region.append(r)
    adata.obs['general_region']=general_region

    #####END REGION HANDLING
    print(adata)
    print(adata.obs)
    adata.obs['timepoint']=adata.obs['timepoint'].astype(float)
    adata=sc_utils.sanitize_types_for_h5(adata)
    adata.raw=adata
    print(adata.obs)
    print(adata.var.dtypes)
    adata.write(newfile)        

adata=sc.read(newfile)
adata.var.index=[x.upper() for x in adata.var.index]
adata.raw=adata

adata=adata[adata.obs['batch_name']!='SRR6061129_P10_Cortex',:]
adata=adata[~adata.obs['batch_name'].str.contains('motor'),:]
#adata=adata[~adata.obs['batch_name'].str.contains('GSM3449'),:]
adata.obs['leiden']='nan'
adata.obs['supervised_name']='nan'

supercell=pd.read_csv('/wynton/group/ye/mtschmitz/macaquedevbrain/MouseAdultAggSupervised.txt')
supercell=supercell.loc[supercell['agg_supervised_name']!='nan',:]
ind=adata.obs.index[adata.obs['full_cellname'].isin(supercell['full_cellname'])]
supercell.index=supercell['full_cellname']
supercell=supercell.loc[ind,:]
if 'leiden' in adata.obs.columns:
    adata.obs['leiden']=adata.obs['leiden'].astype(str)
if 'agg_supervised_name' in adata.obs.columns:
    adata.obs['agg_supervised_name']=adata.obs['agg_supervised_name'].astype(str)
print(supercell)
adata.obs.loc[ind,'leiden']=supercell['leiden']
adata.obs.loc[ind,'supervised_name']=supercell['agg_supervised_name']
adata.obs['leiden']=[re.sub('\.0','',str(x)) for x in adata.obs['leiden']]
print(adata.obs.supervised_name.value_counts())
countmat=adata.obs.astype(str).groupby(['leiden', 'supervised_name']).size().reset_index().pivot(columns='supervised_name', index='leiden', values=0)
print(adata.obs['supervised_name'].unique())


adata.var_names_make_unique()
#adata.write(newfile)
sc.pp.filter_genes(adata,min_cells=10)
#sc.pp.filter_cells(adata,min_counts=1000)
sc.pp.filter_cells(adata,min_genes=min_genes)

ribo_genes=[name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL') ]
adata.obs['percent_ribo'] = np.sum(
adata[:, ribo_genes].X, axis=1) / np.sum(adata.X, axis=1)
mito_genes = [name for name in adata.var_names if name in ['ND1','ND2','ND4L','ND4','ND5','ND6','ATP6','ATP8','CYTB','COX1','COX2','COX3'] or name.startswith('CHRM-') or name.startswith('MT-')]

adata.obs['percent_mito'] = np.sum(
adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
print(adata)
sc.pl.violin(adata,groupby='batch_name',keys=['percent_ribo','percent_mito'],rotation=45,save='ribomito')

adata=adata[adata.obs['percent_ribo'].to_numpy()<.4,:]
adata=adata[adata.obs['percent_mito'].to_numpy()<.15,:]
adata=adata[:,~adata.var.index.isin(mito_genes)]
#adata._inplace_subset_obs(adata.obs['percent_ribo']<.4)
#adata._inplace_subset_obs(adata.obs['percent_mito']<.15)
#adata._inplace_subset_var(~adata.var.index.isin(mito_genes))
print(adata)
sc.pp.filter_genes(adata,min_cells=10)
#sc.pp.filter_cells(adata,min_counts=1000)
sc.pp.filter_cells(adata,min_genes=min_genes)
sc.pl.highest_expr_genes(adata, n_top=20, )
print(adata)
sc.pp.normalize_total(adata,exclude_highly_expressed=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata,n_top_genes=15000,batch_key='dataset_name',subset=False)#'dataset_name'
sc.pp.scale(adata,max_value=10)
sc.pp.pca(adata,n_comps=100)
#sc.pp.neighbors(adata)
bbknn.bbknn(adata,batch_key='batch_name',n_pcs=100,neighbors_within_batch=3)
sc.tl.leiden(adata,resolution=.7)
sc.tl.umap(adata,spread=2)
adata.write('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT_chang_h5ad/KDCbVelocityMouseWbPresupervisionProcessed.h5ad')

for i in ['ctx|ge|bn','de','mb','hb']:
    adata[adata.obs.general_region.str.contains(i),:].write(re.sub('\.h5ad',i+'_subset.h5ad',newfile))

sc.pl.umap(adata,color=['leiden'],save='leiden')
sc.pl.umap(adata,color=['batch_name'],save='batch_name')
sc.pl.umap(adata,color=['supervised_name'],save='supervised_name')
sc.pl.umap(adata,color=['timepoint'],save='timepoint')
sc.pl.umap(adata,color=['region'],save='region')
sc.pl.umap(adata,color=['general_region'],save='general_region')    
sc.pl.umap(adata,color=['dataset_name'],save='dataset_name')
easygenes= ['MKI67','SOX2','AQP4','EDNRB','IL33','PDGFRA','HOPX','ERBB4','CALB1','CALB2','GAD1','GAD2','GADD45G','LHX6','LHX7','LHX8','SST','NKX2-1','MAF','SP8','SP9','PROX1','VIP','CCK','NPY','LAMP5','HTR1A','HTR3A','NR2F1','NR2F2','TOX3','ETV1','SCGN','FOXP1','FOXP2','FOXP4','TH','DDC','SLC18A2','PAX6','MEIS2','SHTN1','ISL1','PENK','DRD1','ADORA2A','CHAT','ZNF503','STXBP6','CRABP1','ZIC1','ZIC2','TAC1','TAC2','TAC3','EBF1','DLX5','GSX2','RBFOX3','PPP1R17','EOMES','DCX','TUBB3','BCL11B','TLE4','FEZF2','SATB2','TBR1','AIF1','RELN','PECAM1','HBZ','ROBO1','ROBO2','ROBO3','ROBO4']
easygenes=[x for x in easygenes if x in adata.var.index]
sc.pl.umap(adata,color=easygenes,use_raw=False,save='FaveGenes')

f = plt.figure()
df_plot = adata.obs.groupby(['supervised_name', 'batch_name']).size().reset_index().pivot(columns='batch_name', index='supervised_name', values=0).apply(lambda g: g / g.sum(),1)
ax = df_plot.plot(kind='bar', legend=False,stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('proportion')
ax.get_figure().savefig(os.path.join(sc.settings.figdir,'supervisedVbatchBar.pdf'), bbox_inches="tight")

f = plt.figure()
df_plot = adata.obs.groupby(['supervised_name', 'timepoint']).size().reset_index().pivot(columns='timepoint', index='supervised_name', values=0).apply(lambda g: g / g.sum(),1)
ax = df_plot.plot(kind='bar', legend=False,stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('proportion')
ax.get_figure().savefig(os.path.join(sc.settings.figdir,'supervisedVtimepointBar.pdf'), bbox_inches="tight")

f = plt.figure()
df_plot = adata.obs.groupby(['supervised_name', 'region']).size().reset_index().pivot(columns='region', index='supervised_name', values=0).apply(lambda g: g / g.sum(),1)
ax = df_plot.plot(kind='bar', legend=False,stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('proportion')
ax.get_figure().savefig(os.path.join(sc.settings.figdir,'supervisedVregionBar.pdf'), bbox_inches="tight")

f = plt.figure()
df_plot = adata.obs.groupby(['supervised_name', 'leiden']).size().reset_index().pivot(columns='leiden', index='supervised_name', values=0).apply(lambda g: g / g.sum(),1)
ax = df_plot.plot(kind='bar', legend=False,stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('proportion')
ax.get_figure().savefig(os.path.join(sc.settings.figdir,'supervisedVleidenBar.pdf'), bbox_inches="tight")

f = plt.figure()
df_plot = adata.obs.groupby(['leiden', 'region']).size().reset_index().pivot(columns='region', index='leiden', values=0).apply(lambda g: g / g.sum(),1)
ax = df_plot.plot(kind='bar', legend=False,stacked=True)
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.ylabel('proportion')
ax.get_figure().savefig(os.path.join(sc.settings.figdir,'regionVleidenBar.pdf'), bbox_inches="tight")


sc_utils.cell_cycle_score(adata)
sc_utils.marker_analysis(adata,variables=['leiden'],markerpath=os.path.expanduser('~/markers.txt')) 
sc_utils.log_reg_diff_exp(adata)
sc_utils.log_reg_diff_exp(adata,obs_name='supervised_name')
adata.write('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT_chang_h5ad/KDCbVelocityMouseWbPresupervisionProcessed.h5ad')



