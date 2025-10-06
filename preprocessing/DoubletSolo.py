import scanpy as sc
import os
import numpy as np
import seaborn
import matplotlib
import pandas as pd
import sys
import scvi
import torch
import scipy 
import tqdm

if torch.cuda.is_available():
    print("GPU is available")
    print("Number of GPUs:", torch.cuda.device_count())
    print("GPU Name:", torch.cuda.get_device_name(0))
else:
    print("GPU is not available")


sys.path.append('/home/matthew.schmitz/utils/mts-utils/')
from genomics import sc_analysis


adata=sc.read_h5ad('/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/HvQvM/HvQvMall_cere_clean.h5ad',backed='r')

#scvi.model.SCVI.setup_anndata(adata,batch_key="batch_name",layer='spliced')
#vae = scvi.model.SCVI(adata,n_hidden=2048,n_layers=3,n_latent=100)
#vae.train(max_epochs=100)
#vae.save('/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/extra/solo_outs/scvi_vae', overwrite=True) 

vae=scvi.model.SCVI.load('/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/extra/solo_outs/scvi_vae',adata=adata)

done_list=[]
for b in reversed(adata.obs['batch_name'].unique()):
    if (adata.obs['batch_name']==b).sum()>1000:
        if not b in done_list:
            try:
                try:
                    print(b)
                    solo_batch_1 = scvi.external.SOLO.from_scvi_model(vae, restrict_to_batch=b)
                    solo_batch_1.train(train_size=0.9,early_stopping=False)
                    predictions=solo_batch_1.predict()
                    pred_list=pd.DataFrame(scipy.special.softmax(predictions.to_numpy(),axis=-1),columns=predictions.columns,index=predictions.index)
                except:
                    solo_batch_1 = scvi.external.SOLO.from_scvi_model(vae, restrict_to_batch=b,)
                    solo_batch_1.train(train_size=0.8,early_stopping=False)
                    predictions=solo_batch_1.predict()
                    pred_list=pd.DataFrame(scipy.special.softmax(predictions.to_numpy(),axis=-1),columns=predictions.columns,index=predictions.index)
                pred_list.to_csv(os.path.join('/allen/programs/celltypes/workgroups/rnaseqanalysis/EvoGen/Team/Matthew/data/taxtest/extra/solo_outs',b+'.csv'))
            except:
                print('FAILFAILFAILFAILFAIL')
    
