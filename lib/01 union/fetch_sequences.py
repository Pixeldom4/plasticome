#!/usr/bin/env python3
"""Fetch amino-acid sequences for union rows whose aa_sequence is blank but which
have an accession. Cascade: NCBI protein efetch -> UniProt REST -> RCSB (PDB) ->
UniParc. Fills the union table in place and prints a report.

Stage 2 of the union build: build_union.py deliberately leaves every non-manual
v1.1-only row blank so this script fills it from the accession. Rows that already
have a sequence are never touched, so manual assignments and v260701's retrieved
sequences pass through untouched. Idempotent.

Usage
-----
  python3 fetch_sequences.py [union.tsv]

Provenance for the accessions resolved this run is written to
fetched_sequences.json next to the union table."""
import argparse, csv, os, sys, time, json, urllib.request, urllib.parse

DEFAULT_UNION = "plasticome_v1.v260701-union.tsv"
EMAIL = "pixeldom04@gmail.com"
KEY = os.environ.get("NCBI_API_KEY", "")

def http(url, data=None, timeout=60):
    req = urllib.request.Request(url, data=data,
        headers={"User-Agent": f"plasticome-seqfetch/1.0 ({EMAIL})"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")

def parse_fasta(text):
    """-> {header_first_token: sequence}"""
    out = {}; acc = None; buf = []
    for line in text.splitlines():
        if line.startswith(">"):
            if acc is not None: out[acc] = "".join(buf)
            hdr = line[1:].strip()
            acc = hdr.split()[0] if hdr else None
            buf = []
        elif acc is not None:
            buf.append(line.strip())
    if acc is not None: out[acc] = "".join(buf)
    return out

def base(a):  # strip version
    return a.split(".")[0]

def ncbi_batch(accs):
    """efetch protein FASTA for a list; return {input_acc: seq}."""
    got = {}
    for i in range(0, len(accs), 50):
        chunk = accs[i:i+50]
        params = {"db": "protein", "id": ",".join(chunk),
                  "rettype": "fasta", "retmode": "text", "email": EMAIL}
        if KEY: params["api_key"] = KEY
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(params)
        try:
            fa = parse_fasta(http(url))
        except Exception as e:
            sys.stderr.write(f"  NCBI chunk err: {e}\n"); fa = {}
        # map returned headers back to inputs by base accession
        by_base = {}
        for hdr, seq in fa.items():
            by_base.setdefault(base(hdr), seq)
        for a in chunk:
            if a in fa: got[a] = fa[a]
            elif base(a) in by_base: got[a] = by_base[base(a)]
        time.sleep(0.34 if not KEY else 0.11)
    return got

def uniprot_one(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{base(acc)}.fasta"
    fa = parse_fasta(http(url))
    return next(iter(fa.values()), None) if fa else None

def rcsb_one(acc):
    url = f"https://www.rcsb.org/fasta/entry/{base(acc).upper()}"
    fa = parse_fasta(http(url))
    # pick the longest chain sequence
    return max(fa.values(), key=len) if fa else None

def uniparc_one(acc):
    """Fallback for accessions DELETED from UniProtKB but archived in UniParc
    (e.g. entries dropped from a reference proteome). The sequence is immutable
    in UniParc, so we take the top cross-referenced UniParc record's FASTA."""
    url = "https://rest.uniprot.org/uniparc/search?" + urllib.parse.urlencode(
        {"query": base(acc), "format": "fasta", "size": 1})
    fa = parse_fasta(http(url))
    return next(iter(fa.values()), None) if fa else None

def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("union", nargs="?", default=DEFAULT_UNION,
                   help=f"union TSV to fill in place (default: ./{DEFAULT_UNION})")
    union = p.parse_args().union

    with open(union) as f:
        rows = list(csv.reader(f, delimiter="\t"))
    h = rows[0]; d = rows[1:]
    ai = h.index("aa_sequence"); acci = h.index("accession")
    targets = sorted({r[acci].strip() for r in d if not r[ai].strip() and r[acci].strip()})
    print(f"targets (missing seq, has accession): {len(targets)}")

    seqs = {}
    # 1) NCBI for everything
    seqs.update({k: v for k, v in ncbi_batch(targets).items() if v})
    print(f"  after NCBI: {len(seqs)}")
    # 2) UniProt for the rest
    for a in [x for x in targets if x not in seqs]:
        try:
            s = uniprot_one(a)
            if s: seqs[a] = s
        except Exception: pass
        time.sleep(0.1)
    print(f"  after UniProt: {len(seqs)}")
    # 3) RCSB (PDB) for the rest
    for a in [x for x in targets if x not in seqs]:
        try:
            s = rcsb_one(a)
            if s: seqs[a] = s
        except Exception: pass
        time.sleep(0.1)
    print(f"  after RCSB: {len(seqs)}")
    # 4) UniParc archive for the rest (deleted-from-UniProtKB accessions)
    for a in [x for x in targets if x not in seqs]:
        try:
            s = uniparc_one(a)
            if s: seqs[a] = s
        except Exception: pass
        time.sleep(0.1)
    print(f"  after UniParc: {len(seqs)}")

    missing = [a for a in targets if a not in seqs]
    print(f"  still unresolved: {len(missing)} -> {missing}")

    # write back -- only ever into blanks, so manual assignments survive
    filled = 0
    for r in d:
        if not r[ai].strip() and r[acci].strip() in seqs:
            r[ai] = seqs[r[acci].strip()]; filled += 1
    with open(union, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(h); w.writerows(d)
    print(f"rows filled: {filled}")
    still_blank = sum(1 for r in d if not r[ai].strip())
    print(f"rows still blank: {still_blank}")
    # persist a small provenance file next to the union table
    prov = os.path.join(os.path.dirname(os.path.abspath(union)), "fetched_sequences.json")
    with open(prov, "w") as f:
        json.dump({"resolved": {a: len(s) for a, s in seqs.items()},
                   "unresolved": missing}, f, indent=2)

if __name__ == "__main__":
    main()
