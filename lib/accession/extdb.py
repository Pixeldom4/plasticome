"""Sequence recovery from non-NCBI databases for Tier S.

Some accessions are valid but live outside NCBI protein: UniProt/TrEMBL and MGnify
identifiers, which Tier 0 correctly labels `external_database`. Their sequence --
needed for the Tier S supplementary match and to fill `aa_sequence` -- must come
from the source database, not NCBI.

  UniProt: `uniprotkb/{acc}.fasta` for active entries; for entries that are
           DELETED/merged (common for old TrEMBL A0A... ids) the active FASTA is
           empty, so fall back to UniParc, which archives the sequence for good.
  MGnify:  the MGYP protein database. (As of this build the MGnify Proteins
           service returns 503 -- "temporarily unavailable due to maintenance" --
           so these currently resolve to nothing; the code path is ready for when
           it returns.)

Cached under .extcache/ so reruns and negative results are free.
"""

from __future__ import annotations

import re
import subprocess
import time
from pathlib import Path

CACHE = Path(__file__).parent / ".extcache"
UA = "plasticome accession-validation (contact pixeldom04@gmail.com)"
UNIPROT = "https://rest.uniprot.org"


def _get(url: str, cache_key: str, *, timeout: int = 40) -> str:
    """Cached GET via curl. Caches the body (empty string for a miss) forever."""
    path = CACHE / cache_key
    if path.exists():
        return path.read_text()
    path.parent.mkdir(parents=True, exist_ok=True)
    time.sleep(0.34)  # be polite to EBI/UniProt
    proc = subprocess.run(
        ["curl", "-sSL", "--http1.1", "--max-time", str(timeout), "-A", UA, url],
        capture_output=True, text=True, timeout=timeout + 15,
    )
    body = proc.stdout if proc.returncode == 0 else ""
    if body.lstrip()[:1] == "<":  # HTML error/maintenance page, not FASTA
        body = ""
    path.write_text(body)
    return body


def _first_fasta_seq(text: str) -> str:
    seq: list[str] = []
    started = False
    for line in text.splitlines():
        if line.startswith(">"):
            if started:
                break          # stop at the second record
            started = True
            continue
        if started:
            seq.append(re.sub(r"[^A-Za-z]", "", line))
    return "".join(seq).upper()


def uniprot_sequence(acc: str) -> tuple[str, str]:
    """(sequence, source) for a UniProt accession; ('', '') if unrecoverable.

    source is 'uniprot' for an active entry, 'uniparc' for an archived one.
    """
    acc = re.sub(r"\.\d+$", "", (acc or "").strip())
    if not acc:
        return "", ""
    seq = _first_fasta_seq(_get(f"{UNIPROT}/uniprotkb/{acc}.fasta", f"uniprot/{acc}.fasta"))
    if seq:
        return seq, "uniprot"
    # Inactive (deleted/merged): UniParc keeps the sequence under a UPI.
    seq = _first_fasta_seq(_get(
        f"{UNIPROT}/uniparc/search?query={acc}&format=fasta&size=1", f"uniparc/{acc}.fasta"))
    if seq:
        return seq, "uniparc"
    return "", ""


def mgnify_sequence(acc: str) -> tuple[str, str]:
    """(sequence, source) for an MGnify MGYP protein accession; ('', '') if down.

    The MGnify Proteins service is under maintenance at build time (503); the FASTA
    endpoint is attempted and any HTML/maintenance body is treated as a miss.
    """
    acc = (acc or "").strip()
    if not acc:
        return "", ""
    for url, key in (
        (f"https://www.ebi.ac.uk/metagenomics/proteins/{acc}.fasta", f"mgnify/{acc}.fasta"),
        (f"https://www.ebi.ac.uk/metagenomics/proteins/{acc}", f"mgnify/{acc}.txt"),
    ):
        seq = _first_fasta_seq(_get(url, key))
        if seq:
            return seq, "mgnify"
    return "", ""


if __name__ == "__main__":
    import sys
    for a in sys.argv[1:]:
        fn = mgnify_sequence if a.upper().startswith("MGY") else uniprot_sequence
        seq, src = fn(a)
        print(f"{a}: source={src or 'none'} len={len(seq)} {seq[:40]}")
