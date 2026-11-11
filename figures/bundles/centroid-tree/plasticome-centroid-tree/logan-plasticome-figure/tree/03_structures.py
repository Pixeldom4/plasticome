#!/usr/bin/env python3
"""
Step 3: one representative structure per component.

Components in different SSN components share under 30 percent identity, which
is where sequence alignment stops being trustworthy, so the comparison between
components is done on structure. Each component is represented by one
AlphaFold DB model.

The representative is the most central tip that has a usable model: walk the
medoid ranking from step 1 and take the first tip whose AlphaFold DB entry
carries the same sequence. UniProt accessions come from the annotation table,
and for tips without one, from a UniParc lookup on the sequence checksum.

A model is usable when its sequence equals the tip's, when one contains the
other with at least 90 percent of the longer one covered (a mature chain
against its precursor), or when the two are at least 90 percent identical over
the longer one (a crystallised variant against the UniProt entry).

If no tip of a component has a usable model, the medoid is folded with ESMFold
through the ESM Atlas API, which takes sequences up to 400 residues. That sends
the sequence to api.esmatlas.com. Pass --no-esmfold to keep everything local.

Components that still have no structure are listed, and step 5 adds them to
the backbone by sequence instead.

Writes work/structures/cNN.pdb and work/structures.tsv. Downloads are cached
in work/structures/cache, so a rerun makes no network calls.
"""

import argparse
import hashlib
import json
import os
import time
import urllib.error
import urllib.request

from Bio import Align
from Bio.Align import substitution_matrices

from common import WORK, WORK_ALL, cdir, read_fasta, read_tsv, write_tsv

SDIR = os.path.join(WORK, "structures")
CACHE = os.path.join(WORK_ALL, "structures", "cache")    # one cache for every set
MIN_COVER = 0.9
MIN_IDENT = 0.9
ESM_MAX = 400
ESM_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"

ALIGNER = Align.PairwiseAligner(mode="local", open_gap_score=-11,
                                extend_gap_score=-1)
ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
MAX_TRIES = 25      # tips tried per component before giving up


def fetch(url, cache_name, binary=False):
    """GET with an on-disk cache. A 404 is cached as an empty file."""
    path = os.path.join(CACHE, cache_name)
    if not os.path.exists(path):
        for attempt in range(5):
            try:
                with urllib.request.urlopen(url, timeout=60) as r:
                    data = r.read()
                break
            except urllib.error.HTTPError as e:
                if e.code in (400, 404):
                    data = b""
                    break
                if attempt == 4:
                    raise
            except (urllib.error.URLError, TimeoutError):
                if attempt == 4:
                    raise
            time.sleep(2 ** attempt)
        with open(path, "wb") as fh:
            fh.write(data)
    with open(path, "rb") as fh:
        data = fh.read()
    return data if binary else data.decode()


def accessions(tip_row, seq):
    """UniProt accessions to try: the annotated one, then UniParc by checksum."""
    out = [a.strip() for a in tip_row["uniprot_accession"].split(";")
           if a.strip() not in ("", "NA", "-")]
    yield from out
    md5 = hashlib.md5(seq.encode()).hexdigest().upper()
    txt = fetch("https://rest.uniprot.org/uniparc/search?query=checksum:"
                f"{md5}&fields=upi,accession&format=json&size=5",
                f"uniparc_{md5}.json")
    if txt:
        for r in json.loads(txt).get("results", []):
            for a in r.get("uniProtKBAccessions") or []:
                a = a.split(".")[0]
                if a not in out:
                    out.append(a)
                    yield a


def match(tip_seq, model_seq):
    if tip_seq == model_seq:
        return "exact"
    short, long_ = sorted((tip_seq, model_seq), key=len)
    if short in long_ and len(short) >= MIN_COVER * len(long_):
        return "tip_in_model" if short == tip_seq else "model_in_tip"
    aln = ALIGNER.align(tip_seq, model_seq)[0]
    same = sum(1 for a, b in zip(aln[0], aln[1]) if a == b and a != "-")
    if same >= MIN_IDENT * len(long_):
        return f"near_{100 * same / len(long_):.1f}"
    return None


def esmfold(seq, cache_name):
    """Fold through the ESM Atlas API, cached. Returns (pdb text, mean pLDDT)."""
    path = os.path.join(CACHE, cache_name)
    if not os.path.exists(path):
        req = urllib.request.Request(ESM_URL, data=seq.encode(), method="POST")
        for attempt in range(5):
            try:
                with urllib.request.urlopen(req, timeout=300) as r:
                    data = r.read()
                break
            except (urllib.error.URLError, TimeoutError):
                if attempt == 4:
                    raise
                time.sleep(2 ** attempt)
        with open(path, "wb") as fh:
            fh.write(data)
    pdb = open(path).read()
    b = [float(line[60:66]) for line in pdb.splitlines()
         if line.startswith("ATOM") and line[12:16].strip() == "CA"]
    plddt = sum(b) / len(b)
    return pdb, plddt * 100 if plddt <= 1.0 else plddt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-esmfold", action="store_true")
    args = ap.parse_args()
    os.makedirs(CACHE, exist_ok=True)
    os.makedirs(SDIR, exist_ok=True)
    tips = {r["tip"]: r for r in read_tsv(os.path.join(WORK, "tips.tsv"))}
    out = []
    for comp in read_tsv(os.path.join(WORK, "components.tsv")):
        c = int(comp["component"])
        seqs = read_fasta(os.path.join(cdir(c), "seqs.faa"))
        ranking = comp["medoid_ranking"].split(";")
        found = None
        for rank, tip in enumerate(ranking[:MAX_TRIES]):
            for acc in accessions(tips[tip], seqs[tip]):
                txt = fetch(f"https://alphafold.ebi.ac.uk/api/prediction/{acc}",
                            f"afdb_{acc}.json")
                if not txt:
                    continue
                for model in json.loads(txt):
                    how = match(seqs[tip], model.get("sequence", ""))
                    if how:
                        found = (rank, tip, acc, model, how)
                        break
                if found:
                    break
            if found:
                break
        medoid = ranking[0]
        if not found and not args.no_esmfold and len(seqs[medoid]) <= ESM_MAX:
            pdb, plddt = esmfold(seqs[medoid], f"esmfold_{medoid}.pdb")
            with open(os.path.join(SDIR, f"c{c:02d}.pdb"), "w") as fh:
                fh.write(pdb)
            out.append([c, comp["fold_group"], medoid,
                        tips[medoid]["enzyme_name"], "ESMFold", "ESM Atlas API",
                        "exact", 0, f"{plddt:.2f}", len(seqs[medoid])])
            continue
        if not found:
            stale = os.path.join(SDIR, f"c{c:02d}.pdb")
            if os.path.exists(stale):
                os.remove(stale)
            out.append([c, comp["fold_group"], ranking[0],
                        tips[ranking[0]]["enzyme_name"], "none",
                        "NA", "NA", "NA", "NA", len(seqs[ranking[0]])])
            print(f"c{c:02d} no usable AlphaFold DB model "
                  f"({min(len(ranking), MAX_TRIES)} tips tried)")
            continue
        rank, tip, acc, model, how = found
        pdb = fetch(model["pdbUrl"], f"{model['modelEntityId']}.pdb", binary=True)
        with open(os.path.join(SDIR, f"c{c:02d}.pdb"), "wb") as fh:
            fh.write(pdb)
        out.append([c, comp["fold_group"], tip, tips[tip]["enzyme_name"], "AFDB",
                    f"{acc} v{model['latestVersion']}", how, rank,
                    model["globalMetricValue"], len(model["sequence"])])
    write_tsv(os.path.join(WORK, "structures.tsv"),
              ["component", "fold_group", "representative", "enzyme_name",
               "source", "model", "sequence_match", "medoid_rank", "mean_plddt",
               "model_length"], out)
    n = sum(1 for r in out if r[4] != "none")
    print(f"{n} of {len(out)} components have a representative structure "
          f"({sum(1 for r in out if r[4] == 'ESMFold')} from ESMFold)")
    low = [r for r in out if r[4] != "none" and float(r[8]) < 70]
    for r in low:
        print(f"  low confidence: c{int(r[0]):02d} {r[3]} mean pLDDT {r[8]}")
    off = [r for r in out if r[4] == "AFDB" and int(r[7]) > 0]
    for r in off:
        print(f"  not the medoid: c{int(r[0]):02d} {r[3]} is rank {r[7]}")


if __name__ == "__main__":
    main()
