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
figdir='/wynton/group/ye/mtschmitz/WBfigures/WbHuman/'
sc.settings.figdir=figdir
scv.settings.figdir=figdir
sc.settings.file_format_figs='pdf'
sc.settings.autosave=True
sc.settings.autoshow=False


pth='/wynton/group/ye/mtschmitz/humanfastqpool/cathg38_gencodev33_kallisto'
dirs=['brainstem','kriegstein','GW16_19','pollen_2022','wang_hippocampus','wang_hypothalamus']

newfile='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT_chang_h5ad/KDCbVelocityPanHumanPresupervise.h5ad'
min_genes=700
if True:#not os.path.exists(newfile):
    adatas=[]
    filename='aem_cellbended_150_750_175e_V0.2'
    for d in dirs:
        print(d,flush=True)
        filepath=os.path.join(pth,d)
        files=os.listdir(filepath)
        files=[f for f in files if 'kOut' in f]
        print(files)
        for f in files:
            print(f,flush=True)
            try:
                if os.path.exists(os.path.join(filepath,f,filename,'aem_cellbended_filtered.h5')):
                    sadata = sc_utils.readCellbenderH5(os.path.join(filepath,f,filename,'aem_cellbended_filtered.h5'))
                    sadata =sadata[sadata.obs['latent_cell_probability']>.99,:]
                    sc.pp.filter_cells(sadata,min_genes=min_genes)
                else:
                    sadata=sc_utils.loadPlainKallisto(os.path.join(filepath,f,'all_em'),min_genes=min_genes)
                print(1,flush=True)
                sadata.obs.index=[re.sub("-1","",x) for x in sadata.obs.index]
                sadata.uns['name']=f
                sadata.obs['batch_name']=str(sadata.uns['name'])
                print(2,flush=True)
                print(2.01,flush=True)
                sadata.obs['timepoint']=sc_utils.tp_format_human(sadata.uns['name'])
                sadata.obs['dataset_name']=d
                print(2.1,flush=True)
                regionstring=sc_utils.region_format_human(sadata.uns['name'])
                print(2.2,flush=True)
                regionstring=regionstring.lower()
                print(3,flush=True)
                sadata.obs['region']=regionstring
                if not os.path.exists(os.path.join(filepath,f,'cbdoublets.txt')):
                    print(4,flush=True)
                    doublets=sc_utils.doublescrub(sadata)
                    print(doublets,flush=True)
                    pd.DataFrame(doublets).to_csv(os.path.join(filepath,f,'cbdoublets.txt'),index=False, header=False)
                try:
                    ddf=pd.read_csv(os.path.join(filepath,f,'cbdoublets.txt'),index_col=False, header=None)
                    doublets=list(ddf[0])
                except:
                    doublets=[]
                sadata=sadata[~sadata.obs.index.isin(doublets),:]
                print(5,flush=True)
                pd.DataFrame(sadata.obs.index).to_csv(os.path.join(filepath,f,'cellbendedcells.txt'),index=False, header=False)
                if sadata.shape[0]>10:
                    adatas.append(sadata)
            except Exception as e:
                print(e)
                print('fail')
    adata=sc.AnnData.concatenate(*adatas)
    adata.var.columns = adata.var.columns.astype(str)
    adata.obs.columns = adata.obs.columns.astype(str)
    adata.obs['clean_cellname']=[re.sub('-[0-9]+','',x) for x in  adata.obs.index]
    adata.obs['full_cellname']=adata.obs['clean_cellname'].astype(str)+'_'+adata.obs['batch_name'].astype(str)
    adata.obs.index=list(adata.obs['full_cellname'])

    gw20fix={'GW20V1':'GW20parietal',
     'GW20somato':'GW20PFC',
     'GW20cingulate-corpus':'GW20hypoventalmedial',
     'GW20claustrum':'GW20MGE',
     'GW20preoptic':'GW20caudalthalamus',
     'GW20motor':'GW20CGE',
     'GW20LGE':'GW20dorsalthalamus',
     'GW20putanum':'GW20ventralthalamus'}
    gw20fixnew={}
    for x in gw20fix.keys():
        gw20fixnew[x+'_kOut']=gw20fix[x]+'_kOut'
        gw20fixnew[gw20fix[x]+'_kOut']=x+'_kOut'
    adata.obs['batch_name']=adata.obs['batch_name'].replace(gw20fixnew)

    ########REGION HANDLING
    from collections import OrderedDict
    reg_dict=OrderedDict({'cerebel|cbc|[^a-z]ce|vermis|cb':'Cb','vz|cortex|ctx|wt[0-9]+|Neurons_Sample_':'Ctx','ob|GSM3449':'OB','head|all|nan|vesicle|placode':'Head','subpall|sp':'FB','thal|pulv|lgn':'Thal',
     'stria|stiatum|putanum|putamen|caud|accumb|nac|basal|SAMN107673':'Str','gp':'GP','clau':'Ctx-Clau','amy':'Amy',
     'ge':'GE','lge':'LGE','mge':'MGE','somato|som':'Ctx-Somato','cge':'CGE','motor|mop|m1|dfc':'Ctx-MOp',
     'pfc|prefront':'Ctx-PFC','cing':'Ctx-Cing','v1|occ':'Ctx-Occ','central|par':'Ctx-Par','temp':'Ctx-Temp',
     'hippo|ca[0-4]|dent|hc':'Ctx-Hip','medull|myelenceph':'Medulla','choroid|cp':'Choroid',
     'sept':'Sept','pre-opt|poa|preoptic':'POA','pons':'Pons','hypo':'Hypo','drg':'DRG','org|cultured':'Cultured','insula':'Ctx-Ins',
     'forebrain|prosen|telenceph':'FB','diencepha':'DE','mesenceph|midbrain|md$|[^a-z]mb':'MB','hind|rhombence|metenceph':'HB'
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
    adata=sc_utils.sanitize_types_for_h5(adata)
    adata.raw=adata

    adata.var.columns = adata.var.columns.astype(str)
    adata.obs.columns = adata.obs.columns.astype(str)
    adata.raw.var.columns = adata.raw.var.columns.astype(str)
    adata.obs['latent_cell_probability']=adata.obs['latent_cell_probability'].fillna(np.nan).astype(float)
    adata.obs['latent_RT_efficiency']=adata.obs['latent_RT_efficiency'].fillna(np.nan).astype(float)
    adata.write(newfile)        


adata=sc.read(newfile)
if True:#not os.path.exists(re.sub('\.h5ad','Processed.h5ad',newfile)):
    adata.raw=adata
    sc.pp.filter_genes(adata,min_cells=10)
    #sc.pp.filter_cells(adata,min_counts=1000)
    sc.pp.filter_cells(adata,min_genes=min_genes)

    adata.var_names_make_unique()
    ribo_genes=[name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL') ]
    adata.obs['percent_ribo'] = np.sum(
    adata[:, ribo_genes].X, axis=1) / np.sum(adata.X, axis=1)
    mito_genes = [name for name in adata.var_names if name in ['ND1','ND2','ND4L','ND4','ND5','ND6','ATP6','ATP8','CYTB','COX1','COX2','COX3'] or name.startswith('chrM-') or name.startswith('MT-')]
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
    #adata[np.random.choice(adata.obs.index,20000,replace=False),:].write('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT202002_h5ad/TinyHuman.h5ad')
    
    print(adata)
    sc.pl.highest_expr_genes(adata, n_top=20, )
    print(adata)
    sc.pp.normalize_total(adata,exclude_highly_expressed=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata,n_top_genes=8000,batch_key='batch_name',subset=False)
    sc.pp.scale(adata,max_value=10)
    sc.pp.pca(adata,n_comps=100)
    #sc.pp.neighbors(adata)
    bbknn.bbknn(adata,batch_key='batch_name',n_pcs=100,neighbors_within_batch=3)
    sc.tl.leiden(adata,resolution=.8)
    sc.tl.umap(adata,spread=2)
    adata.write(re.sub('\.h5ad','Processed.h5ad',newfile))
    #adata=sc.AnnData(adata.raw.X,var=adata.raw.var,obsm=adata.obsm,uns=adata.uns,obs=adata.obs)
    #adata.raw=adata
    #sc.pp.normalize_per_cell(adata)
    #sc.pp.log1p(adata)
    #sc.pp.scale(adata,max_value=10)

    for i in ['ctx|ge|bn','de','mb','hb']:
        adata[adata.obs.general_region.str.contains(i),:].write(re.sub('\.h5ad',i+'_subset.h5ad',newfile))

    sc.pl.umap(adata,color=['leiden'],save='leiden')
    sc.pl.umap(adata,color=['batch_name'],save='batch_name')
    sc.pl.umap(adata,color=['region'],save='region')
    sc.pl.umap(adata,color=['general_region'],save='general_region')    
    sc.pl.umap(adata,color=['timepoint'],save='timepoint')
    easygenes= ['MKI67','SOX2','AQP4','EDNRB','IL33','PDGFRA','HOPX','ERBB4','CALB1','CALB2','GAD1','GAD2','GADD45G','LHX6','SST','NKX2-1','MAF','SP8','SP9','PROX1','VIP','CCK','NPY','LAMP5','HTR3A','NR2F1','NR2F2','TOX3','ETV1','SCGN','FOXP1','FOXP2','FOXP4','TH','PAX6','MEIS2','ISL1','PENK','CRABP1','ZIC1','ZIC2','EBF1','DLX5','GSX2','RBFOX3','PPP1R17','EOMES','DCX','TUBB3','BCL11B','TLE4','FEZF2','SATB2','TBR1','AIF1','RELN','PECAM1','HBZ','ROBO1','ROBO2','ROBO3','ROBO4','FOXG1','EMX1','DLX1','DLX2','DLX5','POMC']
    easygenes=[x for x in easygenes if x in adata.var.index]
    sc.pl.umap(adata,color=easygenes,use_raw=False,save='FaveGenes')

    hierarchy_key='leiden'
    rgs=sc.tl.rank_genes_groups(adata,groupby=hierarchy_key,method='logreg',use_raw=False,copy=True).uns['rank_genes_groups']#,penalty='elasticnet',solver='saga')#or penalty='l1'
    result=rgs
    groups = result['names'].dtype.names
    df=pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores']})
    df.to_csv(os.path.join(sc.settings.figdir,"LogReg"+hierarchy_key+"Norm.csv"))
    topgenes=df.iloc[0:8,['_n' in x for x in df.columns]].T.values
    cols=df.columns[['_n' in x for x in df.columns]]
    cols=[re.sub('_n','',x) for x in cols]
    topdict=dict(zip(cols,topgenes))
    sc.tl.dendrogram(adata,groupby=hierarchy_key)
    var_dict=dict(zip(adata.uns["dendrogram_['"+hierarchy_key+"']"]['categories_ordered'],[topdict[x] for x in adata.uns["dendrogram_['"+hierarchy_key+"']"]['categories_ordered']]))
    sc.pl.matrixplot(adata,groupby=hierarchy_key,var_names=var_dict,cmap='RdBu_r',use_raw=False,dendrogram=True,save='leiden_marker')


    adata.write(re.sub('\.h5ad','Processed.h5ad',newfile))

 


