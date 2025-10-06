#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import scvelo as scv
import anndata
import os
import re
import pandas as pd
import numpy as np
import scipy
import importlib.util
import seaborn as sns
import sklearn
import matplotlib.pyplot as plt
spec = importlib.util.spec_from_file_location("ScanpyUtilsMT", os.path.expanduser("~/code/pollye/MTsc/utils/ScanpyUtilsMT.py"))
sc_utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sc_utils)
sc.settings.autosave=False
sc.settings.autoshow=True


# In[2]:


headpath=os.path.expanduser("/home/mt/Downloads/")
adata=sc.read_h5ad(os.path.join(headpath,'MultiseqRaw.h5ad'))


# In[3]:


adata._inplace_subset_obs([x is not None for x in [re.search('[a-zA-Z]', x) for x in adata.obs['region']]])
adata.obs['region']=[re.sub('_A_|_B_','_',x) for x in adata.obs['region']]
adata.obs['simpleregion']=[re.sub('-1|-2','',x) for x in adata.obs['simpleregion']]
adata.obs['simpleregion']=[re.sub('Hippo$','Hippocampus',x) for x in adata.obs['simpleregion']]
adata.obs['simpleregion']=[re.sub('Hypo$','Hypothalamus',x) for x in adata.obs['simpleregion']]
adata.obs['simpleregion']=[re.sub('Thal$','Thalamus',x) for x in adata.obs['simpleregion']]
adata._inplace_subset_obs([x is None for x in [re.search('nan|Doublet|Negative', x) for x in adata.obs['region']]])


# In[4]:


#data._inplace_subset_obs(adata.obs['batch']=='0')
adata._inplace_subset_obs([x is not None for x in [re.search('[a-zA-Z]', x) for x in adata.obs['region']]])


# In[5]:


adata.var_names_make_unique()

mito_genes = [name for name in adata.var_names if name in ['ND1','ND2','ND4L','ND4','ND5','ND6','ATP6','ATP8','CYTB','COX1','COX2','COX3'] or 'MT-' in name]
ribo_genes = [name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL') ]

adata.obs['percent_mito'] = np.sum(
    adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
adata.obs['percent_ribo'] = np.sum(
    adata[:, ribo_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1
#scanpy.api.pp.filter_genes_dispersion(adata,n_top_genes=np.sum(np.sum(adata.X, axis=0)>0))
sc.pp.filter_genes(adata, min_cells=20,inplace=True)
sc.pp.filter_cells(adata,min_genes=400)
sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito'],
             jitter=0.4, multi_panel=True)

adata.raw = adata.copy() 

sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
varG=sc.pp.highly_variable_genes(adata, n_top_genes=10000,inplace=False)
adata._inplace_subset_var([x[0] for x in varG])
sc.pp.scale(adata,max_value=6)
sc.pp.pca(adata)
sc.pp.neighbors(adata,metric='manhattan')
sc.tl.umap(adata)
sc.tl.leiden(adata,resolution=5)
sc.pl.umap(adata,color='leiden')


# In[6]:


clusterpca=[]
for c in adata.obs['leiden'].unique():
    print(c)
    a=adata.copy()
    a._inplace_subset_obs(a.obs['leiden']==c)
    clusterpca.append(a.obsm['X_pca'].mean(0))
clusterpca=np.array(clusterpca)
print(clusterpca)


# In[7]:



adata.obs['leiden'].unique()[np.argsort([int(x)for x in adata.obs['leiden'].unique()])]


# In[8]:


from scipy.cluster.hierarchy import dendrogram, linkage
from matplotlib import pyplot as plt

clusterpca=clusterpca[np.argsort([int(x)for x in adata.obs['leiden'].unique()]),:]

linked = linkage(clusterpca, 'average')

labelList = np.sort([int(x)for x in adata.obs['leiden'].unique()])

plt.figure(figsize=(10, 7))
dendrogram(linked,
            orientation='top',
            labels=labelList,
            distance_sort='descending',
            show_leaf_counts=True)
plt.show()


# In[9]:


originaladata=adata.copy()


# In[116]:


adata=originaladata.copy()
tree=scipy.cluster.hierarchy.to_tree(linked, rd=False)
def getAllDescendants(root):
    leftdescend=[]
    rightdescend=[]
    queue = [] 
    queue.append(root.left) 
    while(len(queue) > 0): 
        node = queue.pop(0) 
        if node.left is None and node.right is None:
            leftdescend.append(node.id)
        if node.left is not None: 
            queue.append(node.left) 
        if node.right is not None: 
            queue.append(node.right) 
    
    queue.append(root.right) 
    while(len(queue) > 0): 
        node = queue.pop(0) 
        if node.left is None and node.right is None:
            rightdescend.append(node.id)
        if node.left is not None: 
            queue.append(node.left) 
        if node.right is not None: 
            queue.append(node.right)
    return((leftdescend,rightdescend))

hierarchy=[]
queue = [] 
queue.append(tree) 
while(len(queue) > 0): 
    node = queue.pop(0) 
    if node.left is not None and node.right is not None:
        hierarchy.append([node.id,getAllDescendants(node)])
    if node.left.left is not None and node.left.right is not None: 
        queue.append(node.left) 
    if node.right.left is not None and node.right.right is not None: 
        queue.append(node.right) 

maxgenes=8000
ncells=8000
cells=np.random.choice(adata.obs.index,ncells,replace=False)
adata._inplace_subset_obs(cells)
adata._inplace_subset_var((-adata.X).var(0).argsort()[0:min(maxgenes,adata.shape[1])])

from collections import defaultdict
'''des=defaultdict()
for h in hierarchy:
    includeclusters=h[1][0]+h[1][1]
    a=adata.copy()
    a._inplace_subset_obs([int(x) in includeclusters for x in a.obs['leiden']])
    a.obs['clusters']=['left' if int(x) in h[1][0] else 'right' for x in a.obs['leiden'] ]
    sc.pp.highly_variable_genes(a, min_disp=1e-4)
    if(len(a.obs.clusters.unique())>1):
        sc.tl.rank_genes_groups(a, 'clusters', method='t-test_overestim_var',groups=['left','right'])
        result = a.uns['rank_genes_groups']
        groups = result['names'].dtype.names
        df=pd.DataFrame({group + '_' + key[:1]: result[key][group] for group in groups for key in ['names', 'scores']})
        des[str(h[0])]=(df['left_n'][0:6].tolist(),df['right_n'][0:6].tolist())
'''

des=defaultdict()
for h in hierarchy:
    print(h[0])
    des[str(h[0])]=([],[])
    includeclusters=h[1][0]+h[1][1]
    a=adata.copy()
    a._inplace_subset_obs([int(x) in includeclusters for x in a.obs['leiden']])
    a.obs['clusters']=['left' if int(x) in h[1][0] else 'right' for x in a.obs['leiden'] ]
    sc.pp.highly_variable_genes(a, min_disp=1e-4)
    a=a.T
    print('pca calcd')
    aleft=a.copy()
    aleft._inplace_subset_var(a.var['clusters']=='left')
    aright=a.copy()
    aright._inplace_subset_var(a.var['clusters']=='right')
    if aright.shape[1]>20 and aleft.shape[1]>20:
        sc.pp.pca(aleft)
        sc.pp.neighbors(aleft,metric=abscorrelation)
        sc.pp.pca(aright)
        sc.pp.neighbors(aright,metric=abscorrelation)
        print('neighbors calcd')
        def inverseTopNNoncor(l,r,N=10):
            cormat=l.uns['neighbors']['distances']
            nz=cormat.nonzero()
            d=scipy.sparse.csr_matrix(cormat.shape)
            for x,y in zip(nz[0],nz[1]):
                d[x,y]=abscorrelation(r.obsm['X_pca'][x,:],r.obsm['X_pca'][y,:])
            leftdiff=d-cormat
            sortinds=np.argsort(leftdiff.data)
            corpairs=["+".join(sorted([a1,a2])) for a1,a2 in zip(l.obs.index[nz[0][sortinds]],l.obs.index[nz[1][sortinds]])]
            return(list(set(corpairs))[0:N])
        corL=inverseTopNNoncor(aleft,aright)
        corR=inverseTopNNoncor(aright,aleft)
        des[str(h[0])]=(corL,corR)


# In[ ]:


des


# In[118]:


import ete3
def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = ")%s:%.2f%s" % (node.id,parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick
    
import ete3
print(getNewick(tree, "", tree.dist, labelList))
t= ete3.Tree(getNewick(tree, "", tree.dist, labelList), format=1)
ts = ete3.TreeStyle()
ts.show_leaf_name = True
ts.show_branch_length = False
ts.show_branch_support = False
#t.show(tree_style=ts)


# In[119]:


queue = [] 
queue.append(t) 
while(len(queue) > 0): 
    node = queue.pop(0)
    if node.name in des.keys():    
        node.children[1].add_face(ete3.TextFace(' '.join(des[node.name][0])), column=0, position="branch-bottom")
        node.children[0].add_face(ete3.TextFace(' '.join(des[node.name][1])), column=0, position="branch-top")
    if node.children is not None:
        for c in node.children:
            queue.append(c)


# In[120]:


t.show(tree_style=ts)


# In[121]:


t.render("/home/mt/Downloads/HierarchicalCorrelationSmalltest.pdf", w=600, dpi=300, tree_style=ts)


# In[6]:


adata.X=adata.raw[:,[x[0] for x in varG]].X


# In[7]:


gdata=anndata.AnnData(adata.X.T,var=adata.obs,obs=adata.var)
sc.pp.pca(gdata, n_comps=50)


# In[21]:


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
    
sc.pp.neighbors(gdata,metric=abscorrelation)


# In[9]:


sc.tl.umap(gdata)
sc.tl.leiden(gdata,resolution=5)


# In[10]:


sc.pl.umap(gdata,color=['n_cells','leiden'])


# In[11]:


nz=gdata.uns['neighbors']['distances'].nonzero()
newdist=[]

for x,y in zip(nz[0],nz[1]):
    newdist.append(abscorrelation(gdata.obsm['X_pca'][x,:],gdata.obsm['X_pca'][y,:]))
d=scipy.sparse.csr_matrix((newdist,(nz[0],nz[1])))


# In[12]:


sns.distplot(d.data)
sns.distplot(gdata.uns['neighbors']['distances'].data)


# In[13]:


cat_columns = gdata.obs.select_dtypes(['category']).columns
gdata.obs[cat_columns] = gdata.obs[cat_columns].astype(str)
del cat_columns
modulesvd={}
import sklearn
for c in gdata.obs['leiden'].unique():
    svdinstance=sklearn.decomposition.TruncatedSVD(n_components=1)
    svd=svdinstance.fit(gdata[gdata.obs['leiden']==c,:].X)
    modulesvd[c]=svd.components_


# In[14]:


for c in modulesvd.keys():
    adata.obs['eigengene'+c]=modulesvd[c].flatten()


# In[15]:


sc.pl.umap(adata, color=[c for c in adata.obs.columns if 'eigengene' in c])


# In[40]:


print(list(gdata.obs.index[gdata.obs['leiden']=='50']))


# In[16]:


melist=[]
for c in modulesvd.keys():
    melist.append(modulesvd[c].flatten())
sns.clustermap(scipy.stats.spearmanr(np.array(melist).T)[0])


# In[17]:


for c in gdata.obs['leiden'].unique():
    gsub=gdata[gdata.obs['leiden']==c,:].X
    cors=[]
    for i in range(gsub.shape[0]):
        cors.append(np.corrcoef(modulesvd[c][0,:],gsub[i,:].toarray())[0,1])
    sns.distplot(cors)


# In[19]:


gsub=gdata[gdata.obs['leiden']=='20',:].X

cors=[]
for i in range(gsub.shape[0]):
    cors.append(np.corrcoef(modulesvd[c][0,:],gsub[i,:])[0,1])
    print(np.corrcoef(modulesvd[c][0,:],gsub[i,:])[0,1])
    sns.scatterplot(modulesvd[c][0,:],gsub[i,:])
    plt.show()
sns.distplot(cors)


# In[58]:


@numba.njit()
def correlation_vector(x, y):
    mu_x = 0.0
    mu_y = 0.0
    norm_x = 0.0
    norm_y = 0.0
    cov=[]
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
        covi=shifted_x * shifted_y
        cov.append(covi)
        dot_product += covi

    norm_dot=dot_product/float(x.shape[0])
    normprod=np.sqrt(norm_x * norm_y)
    return [(x - norm_dot)/normprod for x in cov]

correlation_vector(np.array([1,2,3]),np.array([2,3,12]))


# In[61]:


for c in gdata.obs['leiden'].unique():
    gsub=gdata[gdata.obs['leiden']==c,:].X
    y=modulesvd[c]
    for i in range(gsub.shape[0]):
        x=gsub[i,:].toarray()
        #print(scipy.stats.pearsonr(x,y.flatten()))
        #sns.scatterplot(x.flatten(),y.flatten())
        #plt.show()
        #The difference between covariance of gene i in cells
        sns.distplot(correlation_vector(x.flatten(),y.flatten()))
        #plt.show()
    plt.show()


# In[ ]:


for i in range(adata.shape[0]):
    print(i)
    nz=adata.uns['neighbors']['distances'][i,:].nonzero()[1]
    g=gdata[:,nz].copy()
    print(nz)
    print(g)


# In[ ]:


distdict={}
pcadict={}
for r in gdata.var['simpleregion'].unique():
    print(r)
    g=gdata[:,gdata.var['simpleregion']==r].copy()
    del g.uns['neighbors']
    del g.obsm['X_pca']
    sc.pp.pca(g)
    sc.pp.neighbors(g,metric='correlation',n_neighbors=10,use_rep='X')
    #nz=g.uns['neighbors']['distances'].nonzero()
    #newdist=[]
    #for x,y in zip(nz[0],nz[1]):
    #    newdist.append(scipy.spatial.distance.correlation(gdata.obsm['X_pca'][x,:],gdata.obsm['X_pca'][y,:]))
    #d=scipy.sparse.csr_matrix((newdist,(nz[0],nz[1])))
    #sns.distplot(d.data)
    #sns.distplot(gdata.uns['neighbors']['distances'].data)
    distdict[r]=g.uns['neighbors']['distances']
    pcadict[r]=g.obsm['X_pca']
    #sns.distplot(d.data-g.uns['neighbors']['distances'].data)
    #plt.show()


# In[23]:


def sparsecordistcontrasts(distdict,pcadict,a,b):
    acor=distdict[a]
    bcor=distdict[b]
    apca=pcadict[a]
    bpca=pcadict[b]
    nz=acor.nonzero()
    newdist=[]
    for x,y in zip(nz[0],nz[1]):
        newdist.append(abscorrelation(bpca[x,:],bpca[y,:]))
    return(acor-scipy.sparse.csr_matrix((newdist,(nz[0],nz[1]))))


# In[24]:


sns.distplot(sparsecordistcontrasts(distdict,pcadict,'MGE','CGE').data)


# In[25]:


sns.distplot(sparsecordistcontrasts(distdict,pcadict,'CGE','MGE').data)


# In[48]:


region1='Thalamus'
region2='Hypothalamus'
a=sparsecordistcontrasts(distdict,pcadict,region1,region2)


# In[51]:


sc.pl.umap(adata,color='simpleregion')
sc.pl.umap(adata[(adata.obs['simpleregion']==region1)| (adata.obs['simpleregion']==region2),:],color='simpleregion')


# In[49]:


sc.set_figure_params(color_map="jet")
ranks=np.argsort(a.data)#scipy.stats.rankdata(a.data,'ordinal')-1
nz=a.nonzero()
sns.distplot([np.corrcoef(gdata.X[nz[0][ranks][i],:], gdata.X[nz[1][ranks][i],:])[0,1] for i in range(100)])
print([a.data[ranks][i] for i in range(100)])
topgenes=[val for pair in zip(gdata.obs.index[nz[0][ranks]], gdata.obs.index[nz[1][ranks]]) for val in pair][0:40]
sc.pl.umap(adata[(adata.obs['simpleregion']==region1)| (adata.obs['simpleregion']==region2),:],color=topgenes)


# In[50]:


ranks=np.argsort(-a.data)#scipy.stats.rankdata(-(a.data),'ordinal')-1
nz=a.nonzero()
print([a.data[ranks][i] for i in range(100)])
sns.distplot([np.corrcoef(gdata.X[nz[0][ranks][i],:], gdata.X[nz[1][ranks][i],:])[0,1] for i in range(100)])
topgenes=[val for pair in zip(gdata.obs.index[nz[0][ranks]], gdata.obs.index[nz[1][ranks]]) for val in pair][0:40]
sc.pl.umap(adata[(adata.obs['simpleregion']==region1)| (adata.obs['simpleregion']==region2),:],color=topgenes)


# In[56]:


clustersvd={}
import sklearn
for c in adata.obs['leiden'].unique():
    svdinstance=sklearn.decomposition.TruncatedSVD(n_components=1)
    svd=svdinstance.fit(adata[adata.obs['leiden']==c,:].X)
    clustersvd[c]=svd.components_


# In[63]:


sklearn np.array(list(clustersvd.values()))[:,0,:]


# In[ ]:


#for c in gdata.obs['leiden'].unique():
#    print(list(gdata.obs.index[gdata.obs['leiden']==c]))


# In[39]:






#sc_utils.cell_cycle_score(adata)
#sc.pl.umap(adata, color=['leiden','n_counts','region','phase'])
sc.pl.umap(adata, color=['region','simpleregion','batch','leiden'],save='_MultiseqSummary')
sc.pl.umap(adata, color=['simpleregion'],save='_region')
sc.pl.umap(adata, color=['leiden'],save='_leiden')
sc.pl.umap(adata, color=['tp'],save='_tp')
sc.pl.umap(adata, color=['ActualRegion'],save='_ActualRegion')


# In[ ]:





# In[43]:


import re
adata.obs['subclassname']=[re.sub('mkrscore','',x) for x in adata.obs.loc[:,['mkrscore' in x for x in adata.obs.columns]].astype('float').idxmax(axis=1)]
normalizedadata.obs['subclassname']=[re.sub('mkrscore','',x) for x in adata.obs.loc[:,['mkrscore' in x for x in adata.obs.columns]].astype('float').idxmax(axis=1)]
sc.pl.umap(adata,color=['subclassname'])


# In[11]:


def most_frequent(List): 
    return max(set(List), key = List.count) 
classlist=[]
for c in adata.obs['subclassname']:
    fullbool=[c in x for x in adata.uns['markers']['fullclass']]
    flatclass=[item for sublist in adata.uns['markers'].loc[fullbool,'type'] for item in sublist]
    classlist.append(most_frequent(flatclass))
adata.obs['classname']=classlist
normalizedadata.obs['classname']=classlist
sc.pl.umap(adata,color=['classname'])


# In[12]:


adata.obs['region']=[re.sub('_A_|_B_','_',x) for x in adata.obs['region']]
normalizedadata.obs['region']=[re.sub('_A_|_B_','_',x) for x in adata.obs['region']]


# In[69]:


def variance(X,axis=1):
    return( (X.power(2)).mean(axis=axis)-(np.power(X.mean(axis=axis),2)) )

