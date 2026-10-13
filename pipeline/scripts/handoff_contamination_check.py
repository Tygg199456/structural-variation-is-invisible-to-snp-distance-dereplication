import os, csv, re, collections, glob
norm=lambda a: re.sub(r"\.\d+$","",a.strip())
prof={norm(r["accession"]):r for r in csv.DictReader(open("results/qc/deletion_profile_all.tsv"),delimiter="\t")}
meta={}
for r in csv.DictReader(open("results/qc/assembly_metadata.tsv"),delimiter="\t"):
    meta.setdefault(norm(r["accession"]),r)
def nseq(k):
    g=glob.glob("data/genomes/%s*/*.fna"%k)
    if not g: return None
    return sum(1 for l in open(g[0]) if l.startswith(">"))
lv=lambda k:(meta.get(k,{}).get("assembly_level") or "")
def run(name,S):
    c=collections.Counter(); miss=0
    for k in S:
        n=nseq(k)
        if n is None: miss+=1; continue
        ab=bool((prof[k].get("replicons_absent") or "").strip())
        if ab and n<=1: c["claims plasmid loss, but n_seq<=1"]+=1
        elif ab and n>1: c["claims plasmid loss, n_seq>1"]+=1
        elif not ab and n<=1: c["no plasmid loss claimed, n_seq<=1"]+=1
        else: c["no plasmid loss, n_seq>1"]+=1
    print("  %-8s n=%-5d no_fasta=%-4d" % (name,len(S),miss))
    for k,v in c.most_common(): print("       %-38s %5d  (%.1f%%)" % (k,v,100*v/(len(S)-miss)))
print("="*94)
print("CONFLATION TEST: does replicons_absent track plasmid SEQUENCES or plasmid CLAIMS?")
print("="*94)
run("P3441", list(prof))
P258=[k for k in prof if lv(k).lower().startswith(("complete","chromosome"))]
run("P258", P258)
print()
for rep,nm in (("NC_003134","pMT1"),("NC_003131","pCD1"),("NC_003132","pPCP1")):
    cls=[k for k in P258 if rep in (prof[k].get("replicons_absent") or "")]
    n1=sum(1 for k in cls if (nseq(k) or 9)<=1)
    print("  P258 %-5s class n=%-3d   n_seq<=1: %d (%.0f%%)" % (nm,len(cls),n1,100*n1/len(cls) if cls else 0))
print()
print("  and the reverse question: how many genomes WITHOUT replicons_absent have n_seq<=1?")
no=[k for k in P258 if not (prof[k].get("replicons_absent") or "").strip()]
n1=sum(1 for k in no if (nseq(k) or 9)<=1)
print("     %d / %d (%.0f%%)" % (n1,len(no),100*n1/len(no) if no else 0))
