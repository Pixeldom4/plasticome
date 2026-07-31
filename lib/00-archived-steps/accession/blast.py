"""Remote NCBI BLAST (blastp) client: submit -> poll -> hit table, cached on disk.

Tier S recovers a suppressed record's live accession from its Identical Protein
Group, but only when the group contains an *identical* live member. When the whole
group is suppressed, or when a live-but-wrong-protein row needs a better-fitting
replacement, exact identity is not enough -- we need *near*-identical neighbours.
That is BLAST's job: align the stuck row's sequence against a protein database and
return the closest current records, which Tier A then ranks by organism + authors.

Transport is the NCBI BLAST URL API over curl (same reasoning as ncbi.py: curl is
the reliable transport against NCBI's edge proxy). Every search is cached under
.blastcache/ keyed by sha256(sequence)+database, because a remote BLAST is a
queued job that costs minutes; a cold run needs network + time, a rerun is free.
An empty result is cached too -- "no near-identical neighbour" is a finding, not a
reason to re-pay the queue on every rerun.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
CACHE = Path(__file__).parent / ".blastcache"
UA = "plasticome-accession-validation (contact: pixeldom04@gmail.com)"

# NCBI etiquette: submit no more than once per ~10s, and do not poll a single RID
# more than once a minute. A blastp against refseq_protein/nr can sit in NCBI's
# queue for many minutes, so we do NOT wait on one search before starting the
# next: queries are submitted in small waves that keep only _MAX_CONCURRENT RIDs
# in flight at once (spaced by _MIN_SUBMIT_GAP), polled every _BATCH_POLL, and a
# new one is submitted each time a slot frees. Submitting 40+ at once trips NCBI's
# per-user cap and *parks* the overflow WAITING indefinitely; a small in-flight
# window drains steadily instead. Override the window with BLAST_CONCURRENCY.
_MIN_SUBMIT_GAP = 11.0
_BATCH_POLL = 20.0
_MAX_CONCURRENT = int(os.environ.get("BLAST_CONCURRENCY", "5") or "5")
_BATCH_DEADLINE = 60 * 60  # overall safety net: give up after 60 minutes

_last_submit = 0.0


def _curl(args: list[str], timeout: int = 90) -> str:
    proc = subprocess.run(
        ["curl", "-sS", "--http1.1", "--max-time", str(timeout), "-A", UA, *args],
        capture_output=True, text=True, timeout=timeout + 15,
    )
    if proc.returncode != 0:
        raise OSError(f"curl exit {proc.returncode}: {proc.stderr.strip()[:200]}")
    return proc.stdout


def _qblast_field(text: str, key: str) -> str:
    """Read a value out of NCBI's `<!--QBlastInfoBegin ... QBlastInfoEnd-->` block.

    The block carries `RID = ...` / `RTOE = ...` on submit and `Status=...` /
    `ThereAreHits=...` on poll, with inconsistent spacing around the `=`.
    """
    m = re.search(rf"{key}\s*=\s*(\S+)", text)
    return m.group(1).strip() if m else ""


def _submit(seq: str, db: str, hitlist: int, expect: float) -> tuple[str, int]:
    global _last_submit
    gap = _MIN_SUBMIT_GAP - (time.monotonic() - _last_submit)
    if gap > 0:
        time.sleep(gap)
    params = [
        "--data-urlencode", "CMD=Put",
        "--data-urlencode", "PROGRAM=blastp",
        "--data-urlencode", f"DATABASE={db}",
        "--data-urlencode", f"QUERY={seq}",
        "--data-urlencode", f"HITLIST_SIZE={hitlist}",
        "--data-urlencode", f"EXPECT={expect}",
    ]
    text = _curl([*params, BLAST_URL])
    _last_submit = time.monotonic()
    rid = _qblast_field(text, "RID")
    rtoe = _qblast_field(text, "RTOE")
    if not rid:
        raise RuntimeError(f"BLAST submit returned no RID: {text[:200]!r}")
    return rid, int(rtoe) if rtoe.isdigit() else 0


def _fetch_ready(rid: str) -> str:
    """Fetch the single-file JSON2 report for a RID already known to be READY.

    JSON2_S is the one machine format the URL API serves *uncompressed* (Tabular
    comes back empty, XML2/JSON2 come back as a zip archive), and it carries the
    per-HSP identity / alignment length / query span Tier A needs.
    """
    return _curl([f"{BLAST_URL}?CMD=Get&FORMAT_TYPE=JSON2_S&RID={rid}"], timeout=120)


def _parse_json2(text: str) -> list[dict]:
    """Parse a JSON2_S BLASTP report into one dict per hit, best bit-score first.

    Each hit: {accession, pident, coverage, evalue, bitscore, title}. pident is
    the best HSP's identities/align_len; coverage is the union of that hit's HSP
    query spans over query_len -- both gate near-identical, near-full-length
    neighbours (the only re-map candidates Tier A wants).
    """
    text = (text or "").strip()
    if not text:
        return []
    try:
        search = json.loads(text)["BlastOutput2"][0]["report"]["results"]["search"]
    except (ValueError, KeyError, IndexError):
        return []
    qlen = search.get("query_len") or 0
    out: list[dict] = []
    for hit in search.get("hits", []):
        descs = hit.get("description") or [{}]
        acc = (descs[0].get("accession") or "").strip()
        hsps = hit.get("hsps") or []
        if not acc or not hsps:
            continue
        top = max(hsps, key=lambda h: h.get("bit_score", 0))
        align_len = top.get("align_len") or 0
        pident = round(100.0 * top.get("identity", 0) / align_len, 1) if align_len else 0.0
        covered = sum(abs(h.get("query_to", 0) - h.get("query_from", 0)) + 1 for h in hsps)
        coverage = round(min(100.0, 100.0 * covered / qlen), 1) if qlen else 0.0
        out.append({
            "accession": acc, "pident": pident, "coverage": coverage,
            "evalue": top.get("evalue", 0.0), "bitscore": top.get("bit_score", 0.0),
            "title": (descs[0].get("title") or "").strip(),
        })
    return sorted(out, key=lambda h: -h["bitscore"])


def _cache_path(db: str, seq: str) -> Path:
    return CACHE / f"{db}_{hashlib.sha256(f'{db}:{seq}'.encode()).hexdigest()[:24]}.json"


def blast_many(items: dict[str, str], *, db: str = "refseq_protein", hitlist: int = 50,
               expect: float = 1e-20, use_cache: bool = True) -> dict[str, list[dict]]:
    """BLAST many sequences, keeping only _MAX_CONCURRENT searches in flight.

    `items` maps a caller key -> amino-acid sequence. Returns key -> hit list
    (best bit-score first). Cached queries are served without a network round
    trip; the rest are submitted in a sliding window and a new one goes out each
    time a slot frees, so NCBI's per-user cap never parks the overflow. An
    empty/failed search yields [] for that key -- never an exception -- so one bad
    query cannot sink the batch.
    """
    CACHE.mkdir(exist_ok=True)
    # BLAST_CACHE_ONLY: serve whatever is cached and treat everything else as "no
    # hits" without submitting. Lets a run fold the results already on disk (and
    # reproduce fully offline) while the remaining searches are warmed separately.
    cache_only = bool(os.environ.get("BLAST_CACHE_ONLY"))
    results: dict[str, list[dict]] = {}
    todo: list[tuple[str, str, Path]] = []
    for key, raw in items.items():
        seq = re.sub(r"[^A-Za-z]", "", raw or "").upper()
        if not seq:
            results[key] = []
            continue
        path = _cache_path(db, seq)
        if use_cache and path.exists():
            results[key] = _parse_json2(path.read_text())
        elif cache_only:
            results[key] = []
        else:
            todo.append((key, seq, path))
    if not todo:
        if cache_only:
            print(f"  blast[{db}]: cache-only; {len(results)} served, uncached -> no hits",
                  file=sys.stderr)
        return results

    print(f"  blast[{db}]: {len(todo)} queries, <= {_MAX_CONCURRENT} in flight",
          file=sys.stderr)
    pending: dict[str, tuple[str, Path]] = {}

    def _fill() -> None:
        while todo and len(pending) < _MAX_CONCURRENT:
            key, seq, path = todo.pop(0)
            try:
                rid, _ = _submit(seq, db, hitlist, expect)
            except Exception as exc:
                print(f"  blast[{db}] {key}: submit failed: {exc}", file=sys.stderr)
                results[key] = []
                continue
            pending[key] = (rid, path)
            print(f"  blast[{db}] {key} RID={rid} submitted "
                  f"({len(pending)} in flight, {len(todo)} waiting)", file=sys.stderr)

    deadline = time.monotonic() + _BATCH_DEADLINE
    _fill()
    while pending and time.monotonic() < deadline:
        time.sleep(_BATCH_POLL)
        for key, (rid, path) in list(pending.items()):
            # Any transient curl error (502/empty reply mid-fetch) must only cost
            # this RID a retry next pass, never sink the whole batch.
            try:
                info = _curl([f"{BLAST_URL}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}"])
                status = _qblast_field(info, "Status")
                if status == "WAITING":
                    continue
                if status == "READY":
                    report = "" if _qblast_field(info, "ThereAreHits").lower() == "no" \
                        else _fetch_ready(rid)
                    path.write_text(report)
                    results[key] = _parse_json2(report)
                    del pending[key]
                    continue
            except OSError:
                continue  # transient network blip; retry this RID next pass
            # FAILED / UNKNOWN(expired) -- record as no hits and move on
            print(f"  blast[{db}] {key} RID={rid}: {status or 'no status'}", file=sys.stderr)
            results[key] = []
            del pending[key]
        _fill()  # top up the freed slots
        if pending or todo:
            print(f"  blast[{db}]: {len(pending)} searching, {len(todo)} waiting",
                  file=sys.stderr)

    for key, (rid, _) in pending.items():
        print(f"  blast[{db}] {key} RID={rid}: not READY before deadline", file=sys.stderr)
        results[key] = []
    for key, _seq, _path in todo:
        results[key] = []
    return results


def blastp(seq: str, *, db: str = "refseq_protein", **kw) -> list[dict]:
    """Single-sequence convenience wrapper over blast_many."""
    return blast_many({"_": seq}, db=db, **kw)["_"]


if __name__ == "__main__":
    # Self-test: BLAST one sequence (given on the CLI, or a default PETase) and
    # print the top hits. Verifies submit -> poll -> hit-table parse end to end.
    if len(sys.argv) > 1 and re.fullmatch(r"[A-Za-z]+", sys.argv[1]):
        query = sys.argv[1]
    else:
        # Ideonella sakaiensis PETase (A0A0K8P6T7 mature region), a safe known query.
        query = (
            "QTNPYARGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGP"
            "RLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAAN"
            "NPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNS"
            "NQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCS"
        )
    for h in blastp(query)[:15]:
        print(f"{h['accession']:<18} pid={h['pident']:>6} cov={h['coverage']:>6} "
              f"e={h['evalue']:.1e} bits={h['bitscore']}")
