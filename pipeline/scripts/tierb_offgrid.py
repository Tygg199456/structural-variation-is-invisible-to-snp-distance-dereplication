#!/usr/bin/env python3
"""Tier B: evaluate the OFF-GRID integers t=18,19 next to the grid points."""
import csv, numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import re
D=np.load('results/dedup/snp_distance_raw_tierB.npy').astype(np.int64)
order=[l.strip() for l in open('results/dedup/snp_distance_raw_tierB.order.txt')]
np.fill_diagonal(D,0)
ss=list(csv.DictReader(open('results/figures/fig1b_same_strain.csv')))
idx={a:i for i,a in enumerate(order)}
Z=linkage(squareform(D.astype(np.float64),checks=False),method='complete')
print("Tier B n=%d pairs=%d"%(len(order), D.shape[0]*(D.shape[0]-1)//2))
for t in (15,16,17,18,19,20,25):
    lab=fcluster(Z,t=t,criterion='distance'); sizes=np.bincount(lab)[1:]
    m=sum(1 for r in ss if r['acc_a'] in idx and r['acc_b'] in idx and lab[idx[r['acc_a']]]==lab[idx[r['acc_b']]])
    diam=0
    top=sorted(sizes,reverse=True)[:1]
    print("  t=%2d  merged=%d/5  max_size=%d  diameter>=%d"%(t,m,int(sizes.max()),t))
