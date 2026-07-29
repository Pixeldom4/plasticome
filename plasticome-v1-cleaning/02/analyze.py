import csv, re
from Bio import Align
from Bio.Align import substitution_matrices

def read_fasta(path):
    d={};h=None
    for l in open(path):
        l=l.rstrip()
        if l.startswith('>'):h=l[1:];d[h]=''
        elif h:d[h]+=l
    return d

src=read_fasta('plasticome-v1-cleaning/02/source_sequences.fasta')
# map: acc -> (enzyme, like, seq)
meta={}
for hdr,seq in src.items():
    acc=hdr.split()[0]
    enz=re.search(r'enzyme=(\S+)',hdr).group(1)
    like=re.search(r'like=(\d)',hdr).group(1)=='1'
    meta[acc]=(enz,like,seq)

# D1 reported
rep={}
with open('D1-erickson.csv',newline='',encoding='utf-8-sig') as f:
    data=list(csv.reader(f))
for i,r in enumerate(data):
    if r and r[0].strip()=='Enzyme ID':
        for r2 in data[i+1:]:
            if r2 and r2[0].strip():rep[r2[0].strip()]=r2[1].strip().upper()
        break

def strip_tag(s):
    m=re.search(r'(LE)?H{5,}\*?$',s);return s[:m.start()] if m else s

aln=Align.PairwiseAligner();aln.substitution_matrix=substitution_matrices.load("BLOSUM62")
aln.mode="local";aln.open_gap_score=-11;aln.extend_gap_score=-1

def analyze(n,rt):
    if rt in n:
        off=n.index(rt);return "IDENTICAL",100.0,100.0,f"signal/N-term trim {off}; C-term trim {len(n)-off-len(rt)}"
    if rt and rt[0]=='M' and rt[1:] in n:
        off=n.index(rt[1:]);return "IDENTICAL",100.0,100.0,f"+start-Met; N-term trim {off}; C-term trim {len(n)-off-len(rt[1:])}"
    a=aln.align(n,rt)[0];c=a.counts();alen=c.identities+c.mismatches+c.gaps
    idp=100*c.identities/alen
    rcov=sum(e-s for s,e in a.aligned[1]);covp=100*rcov/len(rt)
    cls="NEAR-IDENTICAL" if idp>=99 else "HIGH-SIMILARITY" if idp>=90 else "HOMOLOG" if idp>=40 else "DIVERGENT"
    return cls,round(idp,1),round(covp,1),f"{c.mismatches} mismatch, {c.gaps} gap"

# cross-validate overlaps vs set 01 fetched NCBI
prev=read_fasta('plasticome-v1-cleaning/ipet_ncbi_sequences.fasta')
prev={k.split()[0]:v for k,v in prev.items()}
print("=== transcription cross-check vs set-01 NCBI fetch (overlapping accessions) ===")
mismatch_found=False
for acc in meta:
    if acc in prev:
        same = meta[acc][2]==prev[acc]
        if not same:
            mismatch_found=True
            print(f"  DIFF {acc}: 02-len={len(meta[acc][2])} vs set01-len={len(prev[acc])}")
print("  all overlapping sequences identical to set-01 fetch" if not mismatch_found else "  ^ review diffs above")

order=sorted(meta.items(), key=lambda kv:(int(re.sub(r'\D','',kv[1][0])), kv[0]))
rows=[]
print("\n"+f"{'accession':16}{'enz':7}{'lbl':6}{'src':6}{'rep':6}{'class':16}{'id%':7}{'cov%':7} notes")
print("-"*104)
for acc,(enz,like,seq) in order:
    enz_num=re.sub(r'\D','',enz)
    r=rep.get(enz_num)
    if r is None:
        print(f"{acc:16}{enz:7}{'like' if like else '':6}{len(seq):<6}{'--':6}{'NO-D1-ENTRY':16}");continue
    rt=strip_tag(r)
    cls,idp,covp,note=analyze(seq,rt)
    lbl='like' if like else ''
    print(f"{acc:16}{enz:7}{lbl:6}{len(seq):<6}{len(rt):<6}{cls:16}{idp:<7}{covp:<7} {note}")
    rows.append(dict(accession=acc,enzyme_id=enz_num,table_label=(f"Enzyme {enz_num} like" if like else f"Enzyme {enz_num}"),
                     is_like=like,source_len=len(seq),reported_mature_len=len(rt),
                     classification=cls,pct_identity=idp,coverage_pct=covp,notes=note))

out='plasticome-v1-cleaning/02/source_vs_D1_diff.csv'
with open(out,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)
from collections import Counter
print("\nsummary:",dict(Counter(r['classification'] for r in rows)))
print("wrote",out)
