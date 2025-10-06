# various import statements
import os
import pandas as pd
import scanpy as sc
import tqdm
import numpy as np
from collections import defaultdict
import numpy as np
import re
import scipy

retina_path='/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/shekar_retina'
rfiles=os.listdir(retina_path)
adatas={}
for f in rfiles:
    print(f)
    if not f in adatas.keys() and '.csv' in f:
        adata=sc.read_csv(os.path.join(retina_path,f)).T
        adata.obs['file']=f
        adatas[f]=adata

species_names = [re.search(r'_(\w+)_count_mat', filename).group(1) for filename in adatas.keys()]
keydict=dict(zip(adatas.keys(),species_names))
for k in adatas.keys():
    adatas[k].obs['species']=keydict[k]

human_path='/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/shekar_retina/other'
rfiles=os.listdir(human_path)
for f in rfiles:
    print(f)
    if not f in adatas.keys():
        try:
            adata=sc.read_csv(os.path.join(human_path,f)).T
            adata.obs['file']=f
            adatas[f]=adata
        except:
            print("FAILFAL")

species_names = [re.search(r'_(\w+)_count_mat', filename).group(1) for filename in adatas.keys()]
keydict=dict(zip(adatas.keys(),species_names))
for k in adatas.keys():
    adatas[k].obs['species']=keydict[k]


df=pd.read_csv('/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/genomes/cleanome_genomes/all_gene_ids.csv',header=0)

mdict=dict(zip(df.loc[df['species']=='Mus_musculus','gene'],df.loc[df['species']=='Mus_musculus','ortholog_symbol']))

adatas['GSE133382_Mouse_count_mat.csv.gz'].var.index=adatas['GSE133382_Mouse_count_mat.csv.gz'].var.index.to_series().replace(mdict)

adatas['GSE133382_Mouse_count_mat.csv.gz'].var_names_make_unique()

adata=sc.concat(adatas.values())

adata

pattern = re.compile(r'(.+):')
adata.obs['batch']=[match.group(1) if (match := pattern.search(value)) else '' for value in adata.obs.index ]
adata.obs.loc[adata.obs['batch']=='','batch']=adata.obs.loc[adata.obs['batch']=='','file']

adata.layers['UMIs']=scipy.sparse.csr_matrix(adata.X.astype(int))

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

adata.write('/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/shekar_retina/retina.h5ad')

adatas={k:adatas[k] for k in adatas.keys() if not re.search('Lizard|Zebrafish',k)}

adata=sc.concat(adatas.values())

pattern = re.compile(r'(.+):')
adata.obs['batch']=[match.group(1) if (match := pattern.search(value)) else '' for value in adata.obs.index ]
adata.obs.loc[adata.obs['batch']=='','batch']=adata.obs.loc[adata.obs['batch']=='','file']

adata.layers['UMIs']=scipy.sparse.csr_matrix(adata.X.astype(int))

sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

adata.write('/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/shekar_retina/retina_mammals.h5ad')