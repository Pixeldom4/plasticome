"""Supplementary-material harvesting for Tier S.

The citation graph (Tier 2) cannot confirm the mining-paper rows: NCBI links a
protein only to the paper that *deposited* it, never to a later screening study
that merely *analysed* it. But those screening papers list the accessions and
sequences they analysed in their supplementary files. This module fetches and
parses those SI files so Tier S can match a row's accession/sequence against what
the cited paper actually reports.

Reality of the corpus:
  - The SI filenames live in each paper's main-text JATS XML (`<media xlink:href>`),
    which we already have under ../activity-annotation/papers for ~26 papers.
  - NCBI's FTP OA package and the /bin/ path are unreachable from this sandbox
    (FTP blocked, HTTPS mirror 404s), so we download SI straight from the
    publisher. ACS (pubs.acs.org) serves SI files unauthenticated and is where the
    load-bearing 2025 mining paper lives; other publishers are best-effort.

Everything is cached under .supp_cache/ (downloaded files and parsed per-DOI
indices alike), so reruns are free and failures are not re-attempted.
"""

from __future__ import annotations

import json
import re
import subprocess
import time
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

BASE = Path(__file__).parent
PAPERS_DIR = BASE.parent / "activity-annotation" / "papers"
SUPP_CACHE = BASE / ".supp_cache"
FILES_DIR = SUPP_CACHE / "files"
INDEX_DIR = SUPP_CACHE / "index"

# Browser-ish UA: publisher edges (ACS especially) reject generic agents.
UA = "Mozilla/5.0 (plasticome accession-validation; contact pixeldom04@gmail.com)"

# A protein accession as it appears in an SI table. Deliberately broad but
# anchored, covering the flavours the sheet uses (RefSeq, INSD, WGS, UniProt,
# MGnify). Versions are stripped by the caller.
ACCESSION_RE = re.compile(
    r"\b(?:"
    r"[A-Z]{2}_\d{6,9}"                       # RefSeq WP_/NP_/XP_...
    r"|[A-Z]{3}\d{5,7}"                        # INSD protein / WGS-MAG
    r"|MGYP\d{9,}"                             # MGnify
    r"|[OPQ]\d[A-Z0-9]{3}\d|[A-NR-Z]\d[A-Z][A-Z0-9]{2}\d(?:[A-Z][A-Z0-9]{2}\d)?"  # UniProt
    r")(?:\.\d+)?\b"
)

# A cell/line that is (essentially) a bare protein sequence.
SEQ_RE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWYXBZUO]{40,}$")
_XLSX_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_DOCX_NS = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


# --------------------------------------------------------------------------
# fetching
# --------------------------------------------------------------------------

# Magic bytes we accept per SI type. Publishers (ACS/Atypon especially) answer a
# throttled request with 200 + an HTML challenge page; without this check that
# page gets cached as a "supplementary file" and silently parses to nothing.
_ZIP_TYPES = {".xlsx", ".docx", ".zip", ".pptx"}


def _looks_valid(dest: Path) -> bool:
    try:
        head = dest.read_bytes()[:512]
    except OSError:
        return False
    if not head:
        return False
    suffix = dest.suffix.lower()
    if suffix in _ZIP_TYPES:
        return head[:2] == b"PK"
    if suffix == ".pdf":
        return head[:5] == b"%PDF-"
    # Text SI (csv/txt/fasta): reject an HTML challenge/error page.
    lowered = head.lstrip().lower()
    return not (lowered.startswith(b"<!doctype") or lowered.startswith(b"<html"))


def _curl(url: str, dest: Path, referer: str = "", timeout: int = 120) -> bool:
    """Download `url` to `dest`, validating the body is the expected file type.

    Paces requests: publisher edges block bursts, so we sleep before each hit and
    reject (rather than cache) any HTML challenge page returned in place of a file.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["curl", "-sSL", "--http1.1", "--max-time", str(timeout),
           "-A", UA, "-H", "Accept: */*"]
    if referer:
        cmd += ["-e", referer]
    cmd += ["-o", str(dest), "-w", "%{http_code}", url]
    time.sleep(2.0)
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout + 15)
    code = (proc.stdout or "").strip()[-3:]
    ok = code == "200" and dest.exists() and dest.stat().st_size > 0 and _looks_valid(dest)
    if not ok:
        dest.unlink(missing_ok=True)
    return ok


def _si_url(doi: str, filename: str) -> str | None:
    """Publisher download URL for an SI file named in the article XML.

    Only patterns verified to serve unauthenticated are returned; anything else
    yields None and the SI is treated as unavailable rather than guessed at.
    """
    d = doi.lower()
    if d.startswith("10.1021/"):        # ACS
        return f"https://pubs.acs.org/doi/suppl/{doi}/suppl_file/{filename}"
    if filename.startswith(("http://", "https://")):
        return filename
    return None


def _local_paper(doi: str) -> Path | None:
    stem = doi.replace("/", "_")
    for ext in (".xml", ".pdf"):
        p = PAPERS_DIR / f"{stem}{ext}"
        if p.exists():
            return p
    # Filenames occasionally keep the raw slash-as-colon form.
    for p in PAPERS_DIR.glob(f"{doi.split('/')[0]}*"):
        if p.stem.lower().replace(":", "_").replace("/", "_") == stem.lower():
            return p
    return None


def _si_filenames(xml_path: Path) -> list[str]:
    """SI filenames from a JATS `<media xlink:href>` (dedup, order-preserving)."""
    try:
        text = xml_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    seen, out = set(), []
    for m in re.finditer(r'<media[^>]*xlink:href="([^"]+)"', text):
        fn = m.group(1).strip()
        if fn and fn not in seen:
            seen.add(fn)
            out.append(fn)
    return out


# --------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------

def _xlsx_strings(path: Path) -> list[str]:
    """All string values in an .xlsx (shared strings + inline strings)."""
    out: list[str] = []
    try:
        z = zipfile.ZipFile(path)
    except (zipfile.BadZipFile, OSError):
        return out
    with z:
        if "xl/sharedStrings.xml" in z.namelist():
            try:
                root = ET.fromstring(z.read("xl/sharedStrings.xml"))
                for si in root.iter(_XLSX_NS + "si"):
                    out.append("".join(t.text or "" for t in si.iter(_XLSX_NS + "t")))
            except ET.ParseError:
                pass
        for name in z.namelist():
            if name.startswith("xl/worksheets/") and name.endswith(".xml"):
                try:
                    root = ET.fromstring(z.read(name))
                except ET.ParseError:
                    continue
                for is_el in root.iter(_XLSX_NS + "is"):
                    out.append("".join(t.text or "" for t in is_el.iter(_XLSX_NS + "t")))
    return out


def _docx_strings(path: Path) -> list[str]:
    try:
        z = zipfile.ZipFile(path)
    except (zipfile.BadZipFile, OSError):
        return []
    with z:
        if "word/document.xml" not in z.namelist():
            return []
        try:
            root = ET.fromstring(z.read("word/document.xml"))
        except ET.ParseError:
            return []
        return ["".join(t.text or "" for t in p.iter(_DOCX_NS + "t"))
                for p in root.iter(_DOCX_NS + "p")]


def _fasta_sequences(text: str) -> list[str]:
    seqs, cur = [], []
    for line in text.splitlines():
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
                cur = []
        else:
            cur.append(re.sub(r"[^A-Za-z]", "", line))
    if cur:
        seqs.append("".join(cur))
    return [s.upper() for s in seqs if len(s) >= 40]


def _extract(path: Path) -> tuple[set[str], set[str]]:
    """(accessions, sequences) from one SI file. Unparseable types yield empties."""
    accs: set[str] = set()
    seqs: set[str] = set()
    suffix = path.suffix.lower()

    if suffix == ".xlsx":
        cells = _xlsx_strings(path)
    elif suffix == ".docx":
        cells = _docx_strings(path)
    elif suffix in (".csv", ".tsv", ".txt", ".fasta", ".fa", ".faa", ".seq"):
        text = path.read_text(encoding="utf-8", errors="replace")
        cells = re.split(r"[\t,\r\n]+", text)
        seqs.update(_fasta_sequences(text))
    else:
        # .pdf / .doc / .zip and friends: no dependency-free text extraction.
        return accs, seqs

    for cell in cells:
        cell = (cell or "").strip()
        if not cell:
            continue
        for m in ACCESSION_RE.finditer(cell.upper()):
            accs.add(re.sub(r"\.\d+$", "", m.group(0)))
        compact = re.sub(r"\s+", "", cell).upper()
        if SEQ_RE.match(compact):
            seqs.add(compact)
    return accs, seqs


# --------------------------------------------------------------------------
# per-DOI index
# --------------------------------------------------------------------------

def _main_text_accessions(doi: str) -> set[str]:
    """Accessions literally printed in the local main-text XML (weak signal)."""
    p = _local_paper(doi)
    if p is None or p.suffix.lower() != ".xml":
        return set()
    text = p.read_text(encoding="utf-8", errors="replace").upper()
    # Strip PMC citation ids (PMC1234567) which otherwise pollute the match.
    text = re.sub(r"\bPMC\d+\b", " ", text)
    return {re.sub(r"\.\d+$", "", m.group(0)) for m in ACCESSION_RE.finditer(text)}


def build_index(doi: str, *, refresh: bool = False) -> dict:
    """Fetch + parse everything we can for one DOI. Cached per DOI.

    Returns {status, accessions, sequences, si_files, main_text_accessions}.
    status is one of: parsed (SI parsed), main_text_only (no SI, had local XML),
    unavailable (nothing fetchable).
    """
    INDEX_DIR.mkdir(parents=True, exist_ok=True)
    cache = INDEX_DIR / (doi.replace("/", "_") + ".json")
    if cache.exists() and not refresh:
        return json.loads(cache.read_text())

    paper = _local_paper(doi)
    si_names = _si_filenames(paper) if paper and paper.suffix.lower() == ".xml" else []

    si_accs: set[str] = set()
    si_seqs: set[str] = set()
    fetched: list[str] = []
    for fn in si_names:
        url = _si_url(doi, fn)
        if not url:
            continue
        dest = FILES_DIR / doi.replace("/", "_") / fn.split("/")[-1]
        if not dest.exists():
            referer = f"https://pubs.acs.org/doi/{doi}" if doi.lower().startswith("10.1021/") else ""
            if not _curl(url, dest, referer=referer):
                continue
        a, s = _extract(dest)
        if a or s:
            fetched.append(dest.name)
            si_accs |= a
            si_seqs |= s

    main_accs = _main_text_accessions(doi)

    if si_accs or si_seqs:
        status = "parsed"
    elif main_accs:
        status = "main_text_only"
    else:
        status = "unavailable"

    result = {
        "status": status,
        "accessions": sorted(si_accs),
        "sequences": sorted(si_seqs),
        "main_text_accessions": sorted(main_accs),
        "si_files": fetched,
    }
    cache.write_text(json.dumps(result, indent=0))
    return result


if __name__ == "__main__":
    import sys
    for d in sys.argv[1:]:
        idx = build_index(d, refresh=True)
        print(f"{d}: status={idx['status']} "
              f"accs={len(idx['accessions'])} seqs={len(idx['sequences'])} "
              f"main_text_accs={len(idx['main_text_accessions'])} files={idx['si_files']}")
