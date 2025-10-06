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
figdir='/wynton/group/ye/mtschmitz/WBfigures/WbMacaquePresupervise/'
sc.settings.figdir=figdir
scv.settings.figdir=figdir
sc.settings.file_format_figs='pdf'
sc.settings.autosave=True
sc.settings.autoshow=False


def most_frequent(List): 
    return max(set(List), key = List.count)

def variance(X,axis=1):
    return( (X.power(2)).mean(axis=axis)-(np.power(X.mean(axis=axis),2)) )

filepath='/wynton/group/ye/mtschmitz/macaquedevbrain/TOGA_kallisto'#CAT_fixed_kallisto'
files=os.listdir(filepath)
files=[f for f in files if 'ktOut' in f]
files=[f for f in files if '_ktOut' != f]
fileList=files
#fileList=fileList[0:5]
print(fileList)
newfile='/wynton/group/ye/mtschmitz/macaquedevbrain/CAT_chang_h5ad/KDCbVelocityMacaqueWbPresuperviseKT.h5ad'
min_genes=300
#fileList=[x for x in fileList if x != "E100insula_kOut"]
fileList=[x for x in fileList if 'ulti' not in x]

if True:#not os.path.exists(newfile):
    adatas=[]
    filename='t_cellbended_150_750_175e_V0.2'
    for f in fileList:
        print(f,flush=True)
        try:
            if os.path.exists(os.path.join(filepath,f,filename,'t_cellbended_filtered.h5')):
                sadata = sc_utils.readCellbenderH5(os.path.join(filepath,f,filename,'t_cellbended_filtered.h5'))
                sadata =sadata[sadata.obs['latent_cell_probability']>.99,:]
                sc.pp.filter_cells(sadata,min_genes=min_genes)
            else:
                #continue
                sadata=sc_utils.loadPlainKallisto(os.path.join(filepath,f,'t_em'),min_genes=min_genes+200)
            import json
            with open(os.path.join(filepath,f,'run_info.json')) as file:
                data = json.load(file)
            
            sadata.obs["p_pseudoaligned"]=data["p_pseudoaligned"]
            sadata.obs["p_unique"]=data["p_unique"]
                
            sadata.obs.index=[re.sub("-1","",x) for x in sadata.obs.index]
            sadata=sadata[~sadata.obs.index.str.contains('HUM'),:]
            sadata.uns['name']=f
            sadata.obs['file_name']=f
            sadata.uns['name']=sc_utils.macaque_process_irregular_names(sadata.uns['name'])
            sadata.obs['batch_name']=str(sadata.uns['name'])
            sadata.obs['timepoint']=sc_utils.tp_format_macaque(sadata.uns['name'])
            regionstring=sc_utils.region_format_macaque(sadata.uns['name'])
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
    multiseq=sc.read('/wynton/group/ye/mtschmitz/macaquedevbrain/CAT_chang_h5ad/KDCbVelocityMultiseqkt.h5ad')
    multiseq=multiseq[multiseq.obs.index,:]
    adata=sc.AnnData.concatenate(adata,multiseq,index_unique=None)
    ########REGION HANDLING
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
    adata=sc_utils.sanitize_types_for_h5(adata)
    adata.raw=adata

    adata.var.columns = adata.var.columns.astype(str)
    adata.obs.columns = adata.obs.columns.astype(str)
    adata.raw.var.columns = adata.raw.var.columns.astype(str)
    adata.obs['latent_cell_probability']=adata.obs['latent_cell_probability'].fillna(np.nan).astype(float)
    adata.obs['latent_RT_efficiency']=adata.obs['latent_RT_efficiency'].fillna(np.nan).astype(float)
    adata.write(newfile)        

    
if True:#not os.path.exists(re.sub('\.h5ad','Processed.h5ad',newfile)):
    adata=sc.read(newfile)
    print("1",adata.obs)
    supercell=pd.read_csv('/wynton/group/ye/mtschmitz/macaquedevbrain/MIP_annotations/MacaqueAllSupervisedCellNames.txt')
    ind=adata.obs.index[adata.obs['full_cellname'].isin(supercell['full_cellname'])]
    supercell.index=supercell['full_cellname']
    supercell=supercell.loc[ind,:]
    adata.obs.loc[ind,'supervised_name']=supercell['supervised_name']
    adata.obs['supervised_name']=adata.obs['supervised_name'].astype(str)
    
    supercell=pd.read_csv('/wynton/group/ye/mtschmitz/macaquedevbrain/MIP_annotations/MacaqueGEsupervisednames4.txt')
    ind=adata.obs.index[adata.obs['full_cellname'].isin(supercell['full_cellname'])]
    supercell.index=supercell['full_cellname']
    supercell=supercell.loc[ind,:]
    adata.obs.loc[ind,'supervised_name']=supercell['supervised_name']
    adata.obs['supervised_name']=adata.obs['supervised_name'].astype(str)
    
    print("2",adata)
    print(adata.obs['supervised_name'].unique())

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
    for i in ['ctx|ge|bn','de','mb','hb']:
        adata[adata.obs.general_region.str.contains(i),:].write(re.sub('\.h5ad',i+'_subset.h5ad',newfile))
    
    #adata=sc.AnnData(adata.raw.X,var=adata.raw.var,obsm=adata.obsm,uns=adata.uns,obs=adata.obs)
    #adata.raw=adata
    #sc.pp.normalize_per_cell(adata)
    #sc.pp.log1p(adata)
    #sc.pp.scale(adata,max_value=10)
    sc.pl.umap(adata,color=['leiden'],save='leiden',legend_loc='on data')
    sc.pl.umap(adata,color=['batch_name'],save='batch_name')
    sc.pl.umap(adata,color=['region'],save='region')
    sc.pl.umap(adata,color=['general_region'],save='general_region')    
    sc.pl.umap(adata,color=['timepoint'],save='timepoint')
    easygenes= ['MKI67','SOX2','AQP4','EDNRB','IL33','PDGFRA','HOPX','ERBB4','CALB1','CALB2','GAD1','GAD2','GADD45G','LHX6','SST','NKX2-1','MAF','SP8','SP9','PROX1','VIP','CCK','NPY','LAMP5','HTR3A','NR2F1','NR2F2','TOX3','ETV1','SCGN','FOXP1','FOXP2','FOXP4','TH','PAX6','MEIS2','SHTN1','ISL1','PENK','CRABP1','ZIC1','ZIC2','EBF1','DLX5','GSX2','RBFOX3','PPP1R17','EOMES','DCX','TUBB3','BCL11B','TLE4','FEZF2','SATB2','TBR1','AIF1','RELN','PECAM1','HBZ','ROBO1','ROBO2','ROBO3','ROBO4']
    easygenes=[x for x in easygenes if x in adata.var.index]
    sc.pl.umap(adata,color=easygenes,use_raw=False,save='FaveGenes')

    sc.pl.umap(adata,color=['supervised_name'],save='supervised_name')
    countmat=adata.obs.astype(str).groupby(['leiden', 'supervised_name']).size().reset_index().pivot(columns='supervised_name', index='leiden', values=0)
    if 'nan' in countmat.columns:
        countmat=countmat.drop('nan',axis=1)
    leidentosuper=dict(countmat.idxmax(1))
    adata.obs['supervised_name']=adata.obs['supervised_name'].astype(str)
    adata.obs['supervised_name']='nan'
    for i in adata.obs.loc[adata.obs['supervised_name']=='nan',:].index:
        l=adata.obs['leiden'][i]
        adata.obs.loc[i,'supervised_name']=leidentosuper[l]


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

    f = plt.figure()
    df_plot = adata.obs.groupby(['supervised_name', 'region']).size().reset_index().pivot(columns='region', index='supervised_name', values=0).apply(lambda g: g / g.sum(),1)
    ax = df_plot.plot(kind='bar', legend=False,stacked=True)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.ylabel('proportion')
    ax.get_figure().savefig(os.path.join(sc.settings.figdir,'supervisedVregionBar.pdf'), bbox_inches="tight")

    sc_utils.cell_cycle_score(adata)
    sc_utils.marker_analysis(adata,variables=['leiden'],markerpath=os.path.expanduser('~/markers.txt')) 
    sc_utils.log_reg_diff_exp(adata)
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


    adata.uns['iroot'] = np.flatnonzero(adata.obs.index==sc_utils.get_median_cell(adata,'supervised_name','S-phase Dividing Cell'))[0]
    sc.tl.diffmap(adata,n_comps=15)
    sc.tl.dpt(adata,n_branchings=0,n_dcs=15)
    #sc.pl.diffmap(adata,components='all',save='alldiffmap')

    sc.tl.paga(adata,groups='supervised_name')
    sc.pl.paga(adata,layout='eq_tree',save='normal_tree',legend_fontsize=4,edge_width_scale=.4,threshold=np.quantile(adata.uns['paga']['connectivities'].data,.9))
    sc.pl.paga_compare(adata,legend_fontsize=4,edge_width_scale=.4,threshold=np.quantile(adata.uns['paga']['connectivities'].data,.9),save='connectivity')

    adata.write(re.sub('\.h5ad','Processed.h5ad',newfile))

 


