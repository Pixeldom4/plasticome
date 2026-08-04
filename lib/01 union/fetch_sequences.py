#!/usr/bin/env python3
"""Fetch amino-acid sequences for union rows whose aa_sequence is blank but which
have an accession. Cascade: NCBI protein efetch (twice -- see below) -> UniProt
REST -> RCSB (PDB) -> NCBI nuccore CDS translations -> MGnify -> UniParc. Fills
the union table in place and prints a report.

Accession shapes that need more than a plain protein efetch:

  1ETH_A            a PDB entry+chain. NCBI has it, but answers with the header
                    `>pdb|1ETH|A`, so the record has to be mapped back to the
                    `1ETH_A` we asked for (ncbi_keys). RCSB likewise wants the
                    bare entry -- `/fasta/entry/1ETH_A` returns HTTP 200 with the
                    body "No valid PDB IDs were submitted.", which parses to
                    nothing at all rather than raising.
  OP972509          a GenBank *nucleotide* accession. db=protein errors on it
                    outright; the protein is a CDS translation in db=nuccore.
  MGYP000221121644  an MGnify metagenomic prediction, in none of the sequence
                    databases above. MGnify has no fetch API for these, but its
                    landing page is server-rendered and carries the sequence.

The NCBI pass runs twice because a chunk occasionally comes back short a record
or two without raising -- the second pass asks only for what is still missing.

Stage 2 of the union build: build_union.py deliberately leaves every non-manual
v1.1-only row blank so this script fills it from the accession. Rows that already
have a sequence are never touched, so manual assignments and v260701's retrieved
sequences pass through untouched. Idempotent.

Usage
-----
  python3 fetch_sequences.py [union.tsv]

Provenance for the accessions resolved this run is written to
fetched_sequences.json next to the union table."""
import argparse, csv, os, re, sys, time, json, urllib.request, urllib.parse

DEFAULT_UNION = "plasticome_v1.v260701-union.tsv"
EMAIL = "pixeldom04@gmail.com"
KEY = os.environ.get("NCBI_API_KEY", "")

def http(url, data=None, timeout=60):
    req = urllib.request.Request(url, data=data,
        headers={"User-Agent": f"plasticome-seqfetch/1.0 ({EMAIL})"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")

def fasta_records(text):
    """-> [(full_header_line, sequence)] in file order.

    The whole header is kept because RCSB puts spaces inside its pipe-delimited
    fields (`>1ETH_2|Chains B, D|COLIPASE|...`), so the first token alone loses
    the chain list that rcsb_one has to read.
    """
    out = []; hdr = None; buf = []
    for line in text.splitlines():
        if line.startswith(">"):
            if hdr is not None: out.append((hdr, "".join(buf)))
            hdr = line[1:].strip(); buf = []
        elif hdr is not None:
            buf.append(line.strip())
    if hdr is not None: out.append((hdr, "".join(buf)))
    return out

def parse_fasta(text):
    """-> {header_first_token: sequence}"""
    out = {}
    for hdr, seq in fasta_records(text):
        tok = hdr.split()[0] if hdr else None
        if tok: out[tok] = seq
    return out

def base(a):  # strip version
    return a.split(".")[0]

PDB_ACC = re.compile(r"^([0-9][A-Za-z0-9]{3})_([A-Za-z0-9]+)$")

def pdb_parts(acc):
    """('1ETH', 'A') for a PDB entry+chain accession, else None."""
    m = PDB_ACC.match(acc.strip())
    return (m.group(1), m.group(2)) if m else None

def ncbi_keys(hdr):
    """Every accession form a returned FASTA header should match, upper-cased.

    The header is not always the string we asked for: PDB records come back as
    `pdb|1ETH|A` when the request said `1ETH_A`, so without the reconstructed
    entry_chain form the record is fetched and then silently dropped.
    """
    keys = {hdr, base(hdr)}
    parts = hdr.split("|")
    if len(parts) == 3 and parts[0] == "pdb" and parts[1]:
        keys.add(f"{parts[1]}_{parts[2]}" if parts[2] else parts[1])
    return {k.upper() for k in keys if k}

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
        # map returned headers back to inputs: exact, then version-stripped, then
        # the reconstructed PDB entry_chain form
        by_key = {}
        for hdr, seq in fa.items():
            for k in ncbi_keys(hdr):
                by_key.setdefault(k, seq)
        for a in chunk:
            for k in (a, base(a)):
                if k.upper() in by_key:
                    got[a] = by_key[k.upper()]
                    break
        time.sleep(0.34 if not KEY else 0.11)
    return got

def ncbi_nuccore_cds(acc):
    """Protein translation of a GenBank *nucleotide* accession.

    PAZy sometimes records the nucleotide deposit (OP972509) rather than the
    protein it encodes, and db=protein errors on those rather than redirecting.
    fasta_cds_aa returns one record per CDS; these are single-enzyme deposits, so
    take the longest translation.
    """
    params = {"db": "nuccore", "id": acc, "rettype": "fasta_cds_aa",
              "retmode": "text", "email": EMAIL}
    if KEY: params["api_key"] = KEY
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(params)
    fa = parse_fasta(http(url))
    return max(fa.values(), key=len) if fa else None

def uniprot_one(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{base(acc)}.fasta"
    fa = parse_fasta(http(url))
    return next(iter(fa.values()), None) if fa else None

def rcsb_one(acc):
    """FASTA for a PDB accession, honouring an entry_chain suffix.

    The endpoint takes the 4-character entry ID only: `1ETH_A` gets a 200 whose
    body is "No valid PDB IDs were submitted.", which parse_fasta turns into {}
    -- a silent miss, not an error. Records are headed `>1ETH_2|Chains B, D|...`,
    so the chain we asked for is picked out of that field; without a suffix, or
    if the chain is not listed, fall back to the longest chain.
    """
    parts = pdb_parts(acc)
    entry, chain = parts if parts else (base(acc), "")
    recs = fasta_records(http(f"https://www.rcsb.org/fasta/entry/{entry.upper()}"))
    if not recs:
        return None
    if chain:
        for hdr, seq in recs:
            fields = hdr.split("|")
            listed = fields[1].upper() if len(fields) > 1 else ""
            listed = listed.replace("CHAINS", " ").replace("CHAIN", " ")
            if chain.upper() in [c for c in re.split(r"[\s,\[\]]+", listed) if c]:
                return seq
    return max((seq for _, seq in recs), key=len)

MGYP_ACC = re.compile(r"^MGYP\d{12}$", re.I)
MGNIFY_SEQ = re.compile(r'id="proteinSequenceContainer"[^>]*>([A-Za-z\s]*?)</div>', re.S)
MGNIFY_LEN = re.compile(r'<nightingale-sequence[^>]*\blength="(\d+)"')

def mgnify_one(acc):
    """Sequence for an MGnify metagenomic prediction (MGYP…).

    MGnify exposes no fetch API for these -- /metagenomics/api/v1/proteins/<acc>
    is 404 and /metagenomics/proteins-api/ serves the landing page itself -- but
    that landing page is server-rendered and carries the sequence in
    #proteinSequenceContainer. Because this is scraping and not an API, the
    sequence viewer's declared length is checked against what was extracted, so a
    page-layout change is reported rather than silently returning a truncation.

    Only MGYP accessions are attempted; anything else returns None untouched.
    """
    acc = acc.strip()
    if not MGYP_ACC.match(acc):
        return None
    html = http(f"https://www.ebi.ac.uk/metagenomics/proteins/{acc.upper()}/")
    m = MGNIFY_SEQ.search(html)
    if not m:
        return None
    seq = "".join(m.group(1).split())
    declared = MGNIFY_LEN.search(html)
    if declared and int(declared.group(1)) != len(seq):
        sys.stderr.write(f"  MGnify {acc}: extracted {len(seq)} aa but page declares "
                         f"{declared.group(1)} -- page layout changed, skipping\n")
        return None
    return seq or None

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
    src = {}  # accession -> the stage that resolved it, for the provenance file

    def record(stage, hits):
        for a, s in hits.items():
            if s and a not in seqs:
                seqs[a] = s; src[a] = stage
        print(f"  after {stage}: {len(seqs)}")

    def one_by_one(stage, fn):
        """Run a single-accession lookup over whatever is still missing."""
        hits = {}
        for a in [x for x in targets if x not in seqs]:
            try:
                hits[a] = fn(a)
            except Exception:
                pass
            time.sleep(0.1)
        record(stage, hits)

    # 1) NCBI protein for everything, then again for the stragglers: a chunk can
    #    come back short a record or two without raising, and the retry costs one
    #    request against a handful of ids.
    record("NCBI", ncbi_batch(targets))
    retry = [x for x in targets if x not in seqs]
    if retry:
        record("NCBI retry", ncbi_batch(retry))
    # 2) UniProt, 3) RCSB (PDB entry+chain), 4) nucleotide CDS translations,
    # 5) MGnify (MGYP… only), 6) UniParc archive (deleted from UniProtKB).
    one_by_one("UniProt", uniprot_one)
    one_by_one("RCSB", rcsb_one)
    one_by_one("nuccore CDS", ncbi_nuccore_cds)
    one_by_one("MGnify", mgnify_one)
    one_by_one("UniParc", uniparc_one)

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
    by_stage = {}
    for a in seqs:
        by_stage[src[a]] = by_stage.get(src[a], 0) + 1
    print("resolved by: " + ", ".join(f"{k} {v}" for k, v in by_stage.items()))
    # persist a small provenance file next to the union table
    prov = os.path.join(os.path.dirname(os.path.abspath(union)), "fetched_sequences.json")
    with open(prov, "w") as f:
        json.dump({"resolved": {a: len(s) for a, s in seqs.items()},
                   "resolved_by": src,
                   "by_stage": by_stage,
                   "unresolved": missing}, f, indent=2)

if __name__ == "__main__":
    main()
