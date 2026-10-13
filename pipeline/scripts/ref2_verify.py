#!/usr/bin/env python3
"""Like-with-like re-evaluation on the SAME sweep grid, plus the decomposition the
reviewer asked for."""
import csv, gzip, json, os, re, sys
from pathlib import Path
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

SNP="results_ref2/snippy"; QC="results/qc"
CODE={"A":1,"C":2,"G":3,"T":4}; MISS=set("Nn-.*"); SENT=np.int32(10**6)
GRID=[0,3,5,10,15,20,25,30,50,75,100]
norm=lambda a: re.sub(r"\.\d+$","",(a or "").strip())
def omg(p): return gzip.open(p,"rt") if str(p).endswith(".gz") else open(p)
def fa(path):
    seqs,name,buf={},None,[]
    for line in omg(path):
        if line.startswith(">"):
            if name is not None: seqs[name]="".join(buf)
            name,buf=line[1:].split()[0],[]
        else: buf.append(line.strip())
    if name is not None: seqs[name]="".join(buf)
    return seqs
def blocks(seqs,minbp=10000,binbp=5000):
    out=[]
    for chrom,seq in seqs.items():
        n=len(seq); i=0
        while i<n:
            if seq[i] in MISS:
                j=i
                while j<n and seq[j] in MISS: j+=1
                if j-i>=minbp:
                    lo=(i//binbp)*binbp; hi=((j+binbp-1)//binbp)*binbp
                    out.append("%s:%d:%d"%(chrom,lo,hi-lo))
                i=j
            else: i+=1
    return out
accs=sorted(d.name for d in Path(SNP).iterdir() if (d/"snps.vcf").exists())
sites=set()
for a in accs:
    for line in omg(Path(SNP)/a/"snps.vcf"):
        if line.startswith("#"): continue
        f=line.split("\t",3)
        if len(f)>=3: sites.add((f[0],int(f[1])))
sites=sorted(sites); by={}
for col,(c,p) in enumerate(sites): by.setdefault(c,[]).append((p,col))
G=np.zeros((len(accs),len(sites)),dtype=np.int8); blk=[]
for si,a in enumerate(accs):
    seqs=fa(Path(SNP)/a/"snps.aligned.fa"); row=G[si]
    for chrom,poss in by.items():
        s=seqs.get(chrom)
        if s is None: continue
        L=len(s)
        for pos,col in poss:
            if pos<=L:
                ch=s[pos-1]
                if ch not in MISS: row[col]=CODE.get(ch.upper(),0)
    blk.append(blocks(seqs))
cm=(G!=0).astype(np.float32); C=cm@cm.T
oh=np.zeros((G.shape[0],G.shape[1],4),dtype=np.float32)
for k in range(4): oh[:,:,k]=(G==k+1)
ag=np.zeros_like(C)
for k in range(4): ag+=oh[:,:,k]@oh[:,:,k].T
D=np.rint(C-ag).astype(np.int32); np.fill_diagonal(D,0)
raw=D.copy(); D=D.copy(); D[(C/float(len(sites)))<0.70]=SENT

prof={norm(r["accession"]):(r.get("replicons_absent") or "",r.get("big_blocks_5kb") or "")
      for r in csv.DictReader(open(f"{QC}/deletion_profile_all.tsv"),delimiter="\t")}
ts={norm(r["accession"]):(r["pCD1_state"],r["pPCP1_state"],r["pMT1_state"])
    for r in csv.DictReader(open(f"{QC}/replicon_three_state_seqreport.tsv"),delimiter="\t")}
T=[ts.get(norm(a)) for a in accs]
B=[frozenset(blk[i]) for i in range(len(accs))]
PL=[frozenset((prof.get(norm(a),[None,''])[0] or '').split(';')) for a in accs]

def cdiff(i,j): return B[i]!=B[j]
def pclass(i,j):
    ta,tb=T[i],T[j]
    if ta is None or tb is None: return None
    if ta==tb: return "same"
    if any(ta[k]=="UNRESOLVED" or tb[k]=="UNRESOLVED" for k in range(3)
           if ta[k]!=tb[k]): return "unres"
    return "res"
def blocked(i,j,mode):
    if mode=="full": return (B[i]!=B[j]) or (PL[i]!=PL[j]) or (T[i]!=T[j])
    if cdiff(i,j): return True
    pc=pclass(i,j)
    if mode=="threestate": return pc=="res"
    return pc in ("res","unres")

ss=list(csv.DictReader(open("results/figures/fig1b_same_strain.csv")))
idx={a:i for i,a in enumerate(accs)}
out={"reference":"GCF_000013825.1","n":len(accs),"sites":len(sites),"grid":GRID,"rule":[]}
Z=None
for mode in ("off","full","threestate","conservative"):
    Dm=D.copy()
    if mode!="off":
        for i in range(len(accs)):
            for j in range(i+1,len(accs)):
                if Dm[i,j]<SENT and blocked(i,j,mode): Dm[i,j]=Dm[j,i]=SENT
    Zm=linkage(squareform(np.where(Dm>=SENT,1e6,Dm).astype(np.float64),checks=False),method="complete")
    if mode=="off": Z=Zm
    for t in GRID:
        lab=fcluster(Zm,t=t,criterion="distance"); sizes=np.bincount(lab)[1:]
        m=sum(1 for r in ss if r["acc_a"] in idx and r["acc_b"] in idx and lab[idx[r["acc_a"]]]==lab[idx[r["acc_b"]]])
        diam=0
        for c in np.unique(lab):
            mm=np.where(lab==c)[0]
            if mm.size>1: diam=max(diam,int(raw[np.ix_(mm,mm)].max()))
        ok=(m==5) and (diam<=2*t) and (sizes.max()<=25)
        if mode=="off":
            out["rule"].append({"t":t,"merged":m,"diameter":diam,"max_size":int(sizes.max()),"pass":bool(ok)})
    if mode=="off":
        lab=fcluster(Zm,t=20,criterion="distance")
        import collections
        sizes=collections.Counter(lab); out["size_distribution_t20"]=dict(sorted(collections.Counter(sizes.values()).items()))
        # decomposition of multi-member clusters under the RAW profile rule
        comp=collections.Counter()
        for c in np.unique(lab):
            mm=np.where(lab==c)[0]
            if mm.size<2: continue
            kind="not_mixed"
            if any(cdiff(int(mm[x]),int(mm[y])) for x in range(mm.size) for y in range(x+1,mm.size)): kind="a_chromosome"
            elif any(pclass(int(mm[x]),int(mm[y]))=="res" for x in range(mm.size) for y in range(x+1,mm.size)): kind="c_plasmid_resolved"
            elif any(pclass(int(mm[x]),int(mm[y]))=="unres" for x in range(mm.size) for y in range(x+1,mm.size)): kind="b_plasmid_unresolved"
            comp[kind]+=1
        out["off_decomposition_t20"]=dict(comp)
        out["off_multi_total"]=sum(v for k,v in comp.items() if k!="not_mixed")
    out[mode+"_retained_t20"]=int(len(set(fcluster(Zm,t=20,criterion="distance"))))
    if mode=="conservative":
        el=[(i,j) for i in range(len(accs)) for j in range(i+1,len(accs)) if raw[i,j]<=20]
        out["conservative_eligible"]=len(el)
        out["conservative_blocked_eligible"]=sum(1 for i,j in el if blocked(i,j,"conservative"))
    out[mode+"_blocked_total"]=int(sum(1 for i in range(len(accs)) for j in range(i+1,len(accs)) if Dm[i,j]>=SENT and raw[i,j]<SENT))
u=raw[np.triu_indices_from(raw,k=1)]
out["distance_p50"]=float(np.percentile(u,50)); out["distance_p95"]=float(np.percentile(u,95))
print(json.dumps(out,indent=2,sort_keys=True))
open("results_ref2/verify.json","w").write(json.dumps(out,indent=2,sort_keys=True))
