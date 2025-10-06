#!/usr/bin/env python
# coding: utf-8

# In[1]:



import os
import sys
from collections import defaultdict
import gzip
import pandas as pd
import re
import csv
import seaborn
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')


def dataframe(filename):
    """Open an optionally gzipped GTF file and return a pandas.DataFrame.
    """
    # Each column is a list stored as a value in this dict.
    result = defaultdict(list)

    for i, line in enumerate(lines(filename)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                result[key] = [None] * i

        # Ensure this row has some value for each column.
        for key in result.keys():
            result[key].append(line.get(key, None))

    return pd.DataFrame(result)


def lines(filename):
    """Open an optionally gzipped GTF file and generate a dict for each line.
    """
    fn_open = gzip.open if filename.endswith('.gz') else open

    with fn_open(filename) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield parse(line)


def parse(line):
    """Parse a single GTF line and return a dict.
    """
    result = {}

    fields = line.rstrip().split('\t')

    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

    for i, info in enumerate(infos, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            result[key] = _get_value(value)

    return result


def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value


# In[2]:


species={'mouse':'/wynton/group/ye/mtschmitz/refdata2/mm10/Mus_musculus.GRCm38.100.gtf',
        'human':'/wynton/group/ye/mtschmitz/refdata2/hg38/gencodev33/gencode.v33.annotation.gtf',
        'macaque':'/wynton/group/ye/mtschmitz/refdata2/rhemac10/CAT_chang/Rhesus.gtf',
         'cow':'/wynton/group/ye/mtschmitz/refdata2/cetaceans/bosTau/genome/exon_sorted.gtf'}


# In[3]:


feature_names={'mouse':'exon',
        'human':'exon',
        'macaque':'transcript',
         'cow':'exon'}


# In[4]:


gtfs={}
for k in species.keys():
    print(k)
    gtfs[k] = dataframe(species[k])


# In[5]:


gtfs['human'] = dataframe('/wynton/group/ye/mtschmitz/refdata2/hg38/gencodev33/gencode.v33.annotation.gtf')


# In[6]:


gtfs['human']


# In[7]:


orthos=pd.read_csv('/wynton/home/ye/mschmitz1/utils/HOM_AllOrganism.rpt',sep='\t')
orthos=orthos.loc[orthos['NCBI Taxon ID'].isin([10090,9606]),:]
classcounts=orthos['DB Class Key'].value_counts()
one2one=classcounts.index[list(classcounts==2)]
orthos=orthos.loc[orthos['DB Class Key'].isin(one2one),:]

htab=orthos.loc[orthos['NCBI Taxon ID']==9606,:]
mtab=orthos.loc[orthos['NCBI Taxon ID']==10090,:]
genemapping=dict(zip([x.upper() for x in mtab['Symbol']],htab['Symbol']))


# In[8]:


for k in gtfs.keys():
    print(gtfs[k].columns)


# In[9]:


transcript_gtfs={}
for k in gtfs.keys():
    transcript_gtfs[k]=gtfs[k].loc[gtfs[k]['feature']==feature_names[k],['seqname','start','feature','gene_name']]


# In[10]:


for k in gtfs.keys():
    transcript_gtfs[k]=transcript_gtfs[k].loc[~transcript_gtfs[k]['gene_name'].duplicated(),:]


# In[11]:


#convert mouse ids to standard ids
transcript_gtfs['mouse']['gene_name']=transcript_gtfs['mouse']['gene_name'].str.upper().replace(genemapping)


# In[12]:


for k in transcript_gtfs.keys():
    print(transcript_gtfs[k])


# In[13]:


for k in transcript_gtfs.keys():
    print(transcript_gtfs[k]['gene_name'].unique())


# In[14]:


for k in transcript_gtfs.keys():
    transcript_gtfs[k]['start']=transcript_gtfs[k]['start'].astype(int)


# In[15]:


for k in transcript_gtfs.keys():
    transcript_gtfs[k]=transcript_gtfs[k].sort_values(['seqname','start']).reset_index()


# In[16]:


for k in transcript_gtfs.keys():
    transcript_gtfs[k]['linear_rank']=transcript_gtfs[k].index.astype(int)


# In[17]:


for k in transcript_gtfs.keys():
    transcript_gtfs[k].index=transcript_gtfs[k]['gene_name']


# In[18]:


transcript_gtfs


# In[19]:


set(transcript_gtfs['human']['gene_name']).intersection(transcript_gtfs['mouse']['gene_name'])


# In[24]:


species_pairs=[['human','humanscramble'],
               ['human','macaque'],
               ['human','mouse'],
               ['human','cow'],
               ['macaque','cow']]


# In[25]:


for pairs in species_pairs:
    print(pairs)
    if pairs[1]=='humanscramble':
        scrambled_gtf=transcript_gtfs['human'].copy(deep=True)
        scrambled_gtf['linear_rank']=list(np.random.choice(scrambled_gtf['linear_rank'],size=len(scrambled_gtf['linear_rank']),replace=False))
        merged=pd.merge(transcript_gtfs['human'], scrambled_gtf, left_index=True, right_index=True)
    else:
        merged=pd.merge(transcript_gtfs[pairs[0]], transcript_gtfs[pairs[1]], left_index=True, right_index=True)
    seaborn.scatterplot(x='linear_rank_x',y='linear_rank_y',data=merged,alpha=.2,color='slategray')
    plt.title(pairs[0]+' vs ' +pairs[1])
    plt.show()


# In[23]:


transcript_gtfs


# In[45]:


merged.groupby('seqname_y')['seqname_x'].value_counts().unstack().fillna(0)


# In[49]:


seaborn.heatmap(merged.groupby('seqname_y')['seqname_x'].value_counts(normalize=True).unstack().T.fillna(0),cmap='coolwarm')
plt.show()
seaborn.heatmap(merged.groupby('seqname_x')['seqname_y'].value_counts(normalize=True).unstack().fillna(0),cmap='coolwarm')


# In[41]:


merged['seqname_x']=[x[0] for x in merged['seqname_x'].str.split('_')]


# In[65]:


adj_xy=merged.groupby('seqname_y')['seqname_x'].value_counts(normalize=True).T.fillna(0).unstack().stack().reset_index()
adj_yx=merged.groupby('seqname_x')['seqname_y'].value_counts(normalize=True).fillna(0).unstack().stack().reset_index()


# In[77]:


import networkx as nx
G=nx.from_pandas_edgelist(pd.concat([adj_xy,adj_yx],axis=0),source='seqname_y',target='seqname_x')


# In[78]:


top = nx.bipartite.sets(G)[0]

pos = nx.bipartite_layout(G, top)


# In[63]:


adj_xy.stack().reset_index()


# In[ ]:




