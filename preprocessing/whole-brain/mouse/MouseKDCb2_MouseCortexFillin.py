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
#Load my pipeline functions
import importlib
import importlib.util
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../utils/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)
sc.settings.file_format_figs='pdf'
sc.settings.autosave=True
sc.settings.autoshow=False


pth='/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/HvQvM/extra_mouse_ctx'
dirs=['fastq']

newfile='/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/HvQvM/extra_mouse_ctx/mouse_ctx_fillin.h5ad'
min_genes=800
adatas=[]
filename='aem_cellbended_150_750_200e_V0.2'
for d in dirs:
    print(d,flush=True)
    filepath=os.path.join(pth,d)
    files=os.listdir(filepath)
    files=[f for f in files if 'kOut' in f]
    for f in files:
        print(os.path.join(filepath,f,filename,'aem_cellbended_filtered.h5'),flush=True)
        if os.path.exists(os.path.join(filepath,f,filename,'aem_cellbended_filtered.h5')):
            sadata = sc_utils.readCellbenderH5(os.path.join(filepath,f,filename,'aem_cellbended_filtered.h5'))
            sadata =sadata[sadata.obs['latent_cell_probability']>.99,:]
            sc.pp.filter_cells(sadata,min_genes=min_genes)
        else:
            # sadata=sc_utils.loadPlainKallisto(os.path.join(filepath,f,'all_em'),min_genes=min_genes)
            print('FAIL',flush=True)
            continue
        print(1,flush=True)
        sadata.obs.index=[re.sub("-1","",x) for x in sadata.obs.index]
        sadata.uns['name']=f
        sadata.obs['batch_name']=str(sadata.uns['name'])
        print(2,flush=True)
        print(2.01,flush=True)
        sadata.obs['timepoint']=sc_utils.tp_format_mouse(sadata.uns['name'])
        sadata.obs['dataset_name']='10.1038/s41586-021-03670-5'
        print(5,flush=True)
        pd.DataFrame(sadata.obs.index).to_csv(os.path.join(filepath,f,'cellbendedcells.txt'),index=False, header=False)
        if sadata.shape[0]>10:
            adatas.append(sadata)
        #except Exception as e:
        #    print(e)
        #    print('fail')
adata=sc.concat(adatas)
adata.var.columns = adata.var.columns.astype(str)
adata.obs.columns = adata.obs.columns.astype(str)
adata.obs['clean_cellname']=[re.sub('-[0-9]+','',x) for x in  adata.obs.index]
adata.obs['general_region'] = 'ctx'
adata.obs['region'] = 'Ctx'
adata.obs['female'] = 'mixed'

adata.var['mgi_symbol'] = list(adata.var.index)
mgi_tab = pd.read_csv('/home/matthew.schmitz/utils/HOM_AllOrganism.rpt',sep='\t')
mgi_tab.index = mgi_tab['DB Class Key']
join_tab = mgi_tab.loc[mgi_tab['Common Organism Name']=='mouse, laboratory',:].join(mgi_tab.loc[mgi_tab['Common Organism Name']=='human',:],rsuffix='_Human')
m2h = dict(zip(join_tab['Symbol'].str.upper(),join_tab['Symbol_Human'].str.upper()))
adata.var.index = list(adata.var['mgi_symbol'].str.upper().replace(m2h))
adata.obs.index = [str(x) for x in adata.obs.index]
adata.var.index = [str(x) for x in adata.var.index]
adata.obs_names_make_unique()
adata.var_names_make_unique()

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

div_markers={'G0':[],#'YWHAG','NNAT'
             'G1':['CCND1','CCND2','CCND3','CCNE1','CCNE2','MKI67','PCNA'],
'G1S':['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7',
      'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'RAD51AP1', 'GMNN', 
      'WDR76', 'SLBP', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 
      'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 
      'E2F8','H4C3','H1-3','H1-5','H2AC12','H2AC20','H1-2','H1-4','H1-1','FBXO5',
      'SPC25','FAM83D','HIST1H1A','HIST1H1B','HIST1H1C','HIST1H1D','HIST1H1E','CCND1','CCND2','CCND3'],#,'CCNE1', 'CCNE2'
'G2M':['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2', 
       'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 
       'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 
       'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 
       'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 
       'GAS2L3', 'CBX5','ASPM','CENPA','CCND1','CCND2','CCND3']
        }

adata.obs['n_genes'] = (adata.layers['spliced']>0).sum(1)
adata.obs['log10_n_genes'] = np.log10((adata.layers['spliced']).sum(1))
adata.obs['log10_n_counts'] = np.log10((adata.layers['spliced']>0).sum(1))

sc.tl.score_genes_cell_cycle(adata, s_genes=div_markers['G1S'], g2m_genes=div_markers['G2M'])
adata.obsm['phase_sex']=np.concatenate([adata.obs['S_score'].to_numpy().reshape(-1,1),adata.obs['G2M_score'].to_numpy().reshape(-1,1),adata.obs['log10_n_counts'].to_numpy().reshape(-1,1),np.zeros([adata.shape[0],2])],axis=1)
adata.obsm['phase_sex']=adata.obsm['phase_sex']/np.abs(adata.obsm['phase_sex']).max(0)#Scale so max value is magnitude 1
adata.obs.index = adata.obs['batch_name'] + '_' + adata.obs.index.to_series()

adata.write('/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/HvQvM/extra_mouse_ctx/mouse_ctx_fillin.h5ad')
