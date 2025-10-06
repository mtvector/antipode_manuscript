import re
import scanpy
import fuzzywuzzy
import pandas as pd
import numpy as np
#https://www.biostars.org/p/317524/
import csv
import sys
import os
import gzip

def processfq(lines=None):
    ks = ['name', 'sequence', 'optional', 'quality']
    return {k: v for k, v in zip(ks, lines)}

#fq1='/wynton/group/ye/mtschmitz/macaquedevbrain/MULTISEQmacaque/pollena-MULTIseq_E65-1/MULTIseq_E65-1_S99_L002_R1_001.fastq.gz'
#fq2='/wynton/group/ye/mtschmitz/macaquedevbrain/MULTISEQmacaque/pollena-MULTIseq_E65-1/MULTIseq_E65-1_S99_L002_R2_001.fastq.gz'
fq1=sys.argv[1]
fq2=re.sub('_R1_','_R2_',fq1)
#indexfile='/wynton/group/ye/mtschmitz/macaquedevbrain/MULTISEQmacaque/MultiseqIndices.txt'
indexfile=sys.argv[2]
#cellfile='/wynton/scratch/mtschmitz/fastqpool/E80-2019_Multi-seq_kOut/cellbendedcells.txt'
cellfile=sys.argv[3]
Indices=pd.read_csv(indexfile,sep='\t')
Indices.index=Indices.index+1

cells=[]
with open(cellfile) as f:
    for l in f.readlines():
        if len(l)>2:
            cells.append(re.sub('\n','',l))

cellset=list(set(cells))
inds=list(Indices['Barcode Sequence'])
inddictrev=dict(enumerate(Indices['Barcode Sequence']))
inddict=dict(zip(inddict.values(),inddict.keys()))

celldictrev=dict(enumerate(cellset))
celldict=dict(zip(celldictrev.values(),celldictrev.keys()))

mat=np.zeros((len(cellset),len(inds)),dtype=np.int32)

#Slim processing
#Slim processing
n = 4
readpairs=[]
unassigned=[]
i=0
with gzip.open(fq1, 'r') as fh:
    with gzip.open(fq2, 'r') as rh:
        flines = []
        rlines = []
        linef = fh.readline()
        liner = rh.readline()
        while linef and liner:
            i+=1
            flines.append(linef.rstrip())
            rlines.append(liner.rstrip())
            if (len(flines) == n) and (len(rlines) == n):
                recordf = processfq(flines)
                recordr = processfq(rlines)
                indexseq=recordr['sequence'][0:8].decode('utf-8')
                ind=process.extractOne(indexseq,inds,score_cutoff=75)
                if ind is not None:
                    cbseq=recordf['sequence'][0:16].decode('utf-8')
                    cb=process.extractOne(cbseq,cellset,score_cutoff=80)
                    if cb is not None:
                        mat[celldict[cb[0]],inddict[ind[0]]]+=1
                    else:
                        unassigned.append([cbseq,indexseq])
                flines = []
                rlines = []
            if i%100000 == 0:
                print(i)
            linef = fh.readline()
            liner = rh.readline()

pd.DataFrame(mat,index=celldict.keys(),columns=[int(x)+1 for x in inddict.keys()]).to_csv(outfile,header=True,index=True,sep='\t',quoting=csv.QUOTE_NONE)