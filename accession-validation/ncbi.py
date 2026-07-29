"""NCBI E-utilities client: rate-limited, retrying, and cached on disk.

Every network call is cached under .cache/ keyed by endpoint+params, so reruns
of the tier scripts are free. Delete .cache/ to force a refresh.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
import sys
import threading
import time
import urllib.parse
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CACHE_DIR = Path(__file__).parent / ".cache"

# API key is read from NCBI_ENV_VAR (NCBI_API_KEY kept as a fallback name).
# Unset simply means we run at the anonymous 3 req/s limit rather than 10 req/s.
API_KEY = (os.environ.get("NCBI_ENV_VAR") or os.environ.get("NCBI_API_KEY") or "").strip()
EMAIL = os.environ.get("NCBI_EMAIL", "").strip()
TOOL = "plasticome-accession-validation"

# NCBI allows 3 req/s anonymously, 10 req/s with an API key.
_MIN_INTERVAL = 1 / 9.0 if API_KEY else 1 / 2.7


class _RateLimiter:
    def __init__(self, min_interval: float) -> None:
        self.min_interval = min_interval
        self._lock = threading.Lock()
        self._last = 0.0

    def wait(self) -> None:
        with self._lock:
            delta = time.monotonic() - self._last
            if delta < self.min_interval:
                time.sleep(self.min_interval - delta)
            self._last = time.monotonic()


_limiter = _RateLimiter(_MIN_INTERVAL)


def _cache_path(endpoint: str, params: list[tuple[str, str]]) -> Path:
    material = json.dumps([endpoint, sorted(params)], sort_keys=True)
    digest = hashlib.sha256(material.encode()).hexdigest()[:24]
    return CACHE_DIR / f"{endpoint}_{digest}.txt"


def _post(endpoint: str, params: list[tuple[str, str]], timeout: int = 60) -> str:
    """POST to an E-utilities endpoint. POST (not GET) so long id lists fit.

    Shells out to curl rather than using urllib. NCBI's edge proxy intermittently
    answers Python's TLS stack with a 200 + `Transfer-Encoding: chunked` header
    and then zero body bytes, which surfaces as IncompleteRead. The same request,
    byte-identical down to the headers and ALPN protocol, succeeds from curl.
    Reproduced against elink; curl is the reliable transport here.
    """
    body = urllib.parse.urlencode(params)
    proc = subprocess.run(
        [
            "curl", "-sS", "--fail-with-body",
            # NCBI's edge proxy resets HTTP/2 streams mid-response
            # ("INTERNAL_ERROR"); HTTP/1.1 is stable.
            "--http1.1",
            "--max-time", str(timeout),
            "--retry", "0",
            "-A", f"{TOOL} (contact: {EMAIL or 'n/a'})",
            "-d", body,
            f"{EUTILS}/{endpoint}.fcgi",
        ],
        capture_output=True,
        text=True,
        timeout=timeout + 15,
    )
    if proc.returncode != 0:
        raise OSError(f"curl exit {proc.returncode}: {proc.stderr.strip()[:200]}")
    if not proc.stdout.strip():
        raise OSError("empty response body")
    return proc.stdout


def call(
    endpoint: str,
    params: list[tuple[str, str]],
    *,
    use_cache: bool = True,
    retries: int = 4,
    timeout: int = 30,
) -> str:
    """Rate-limited, retrying, cached E-utilities call. Returns the raw body.

    Both outcomes are cached. Caching *failures* matters as much as caching
    successes: NCBI rejects certain id sets deterministically, and without a
    negative cache every rerun re-pays `retries * timeout` seconds re-discovering
    that before it can fall back to a working strategy.
    """
    base = [("tool", TOOL)]
    if EMAIL:
        base.append(("email", EMAIL))
    full = base + params

    path = _cache_path(endpoint, full)
    fail_path = path.with_suffix(".fail")
    if use_cache and path.exists():
        return path.read_text()
    if use_cache and fail_path.exists():
        raise RuntimeError(f"{endpoint}: known-bad request (cached): {fail_path.read_text()[:120]}")

    if API_KEY:
        full = full + [("api_key", API_KEY)]

    last_err: Exception | None = None
    for attempt in range(retries):
        _limiter.wait()
        try:
            text = _post(endpoint, full, timeout=timeout)
        except (OSError, subprocess.TimeoutExpired) as exc:
            last_err = exc
        else:
            # NCBI returns 200 with an HTML error page when it throttles.
            if text.lstrip().startswith("<!DOCTYPE html") or "Too Many Requests" in text[:500]:
                last_err = RuntimeError("throttled by NCBI")
            else:
                CACHE_DIR.mkdir(exist_ok=True)
                path.write_text(text)
                return text
        if attempt < retries - 1:
            backoff = 2**attempt
            print(f"  [retry {attempt + 1}/{retries}] {endpoint}: {last_err}; sleeping {backoff}s",
                  file=sys.stderr)
            time.sleep(backoff)

    CACHE_DIR.mkdir(exist_ok=True)
    fail_path.write_text(str(last_err))
    raise RuntimeError(f"{endpoint} failed after {retries} attempts: {last_err}")


def iter_json_docs(text: str):
    """Yield every JSON document in `text`.

    efetch(retmode=json) streams *concatenated* JSON objects for large id lists
    -- one document per internal chunk, with no separator and no enclosing array.
    json.loads() chokes on that with "Extra data", so decode incrementally.
    """
    # strict=False: NCBI embeds raw newlines/tabs inside "ERROR" strings when it
    # reports a C++ exception, which strict JSON forbids.
    decoder = json.JSONDecoder(strict=False)
    idx, n = 0, len(text)
    while idx < n:
        while idx < n and text[idx] in " \t\r\n":
            idx += 1
        if idx >= n:
            return
        obj, end = decoder.raw_decode(text, idx)
        yield obj
        idx = end


def purge_cache(endpoint: str, params: list[tuple[str, str]]) -> None:
    """Drop a cached response so a poisoned body is not reused."""
    base = [("tool", TOOL)]
    if EMAIL:
        base.append(("email", EMAIL))
    _cache_path(endpoint, base + params + [("retmode", "json")]).unlink(missing_ok=True)


def call_json(endpoint: str, params: list[tuple[str, str]], **kw) -> dict:
    """Return the first JSON document (use call_json_all for streamed replies)."""
    docs = call_json_all(endpoint, params, **kw)
    return docs[0] if docs else {}


def call_json_all(
    endpoint: str, params: list[tuple[str, str]], *, reject_errors: bool = False, **kw
) -> list[dict]:
    """Fetch and decode. With reject_errors, an NCBI `ERROR` payload is a failure.

    NCBI can answer 200 with `{"linksets": [...], "ERROR": "NCBI C++ Exception..."}`
    -- sometimes alongside *partial* linksets. Parsed naively that reads as "these
    proteins have no linked papers", turning a server fault into biological
    evidence. Callers that interpret absence must set reject_errors and let the
    retry/bisect path handle it. The poisoned body is evicted from the cache too.
    """
    qparams = params + [("retmode", "json")]
    text = call(endpoint, qparams, **kw)
    try:
        docs = list(iter_json_docs(text))
    except json.JSONDecodeError:
        # NCBI's edge proxy sometimes closes the stream early, leaving a truncated
        # body that `call` cached as a success -- so every rerun re-reads the same
        # unparseable text. Re-fetch once bypassing the cache to self-heal before
        # giving up.
        text = call(endpoint, qparams, **{**kw, "use_cache": False})
        try:
            docs = list(iter_json_docs(text))
        except json.JSONDecodeError as exc:
            raise RuntimeError(f"{endpoint} returned non-JSON: {text[:200]!r}") from exc

    if reject_errors:
        for doc in docs:
            if doc.get("ERROR"):
                purge_cache(endpoint, params)
                raise RuntimeError(f"{endpoint} returned ERROR: {str(doc['ERROR'])[:120]}")
    return docs


def chunked(seq: list, size: int):
    for i in range(0, len(seq), size):
        yield seq[i : i + size]


def efetch_docsum(accessions: list[str], db: str = "protein") -> dict[str, dict]:
    """Resolve accessions -> docsum records.

    efetch(rettype=docsum) is used instead of esummary because it is the only
    endpoint that (a) accepts bare accessions of every flavour we see, including
    WGS/MAG proteins with no GI (e.g. MBX3625601.1), and (b) silently drops
    unresolvable ids instead of aborting the whole batch. Returns a dict keyed by
    the versionless caption; ids absent from the result did not resolve.
    """
    out: dict[str, dict] = {}
    if not accessions:
        return out
    docs = call_json_all(
        "efetch",
        [("db", db), ("id", ",".join(accessions)), ("rettype", "docsum")],
    )
    for data in docs:
        result = data.get("result", {})
        for uid in result.get("uids", []):
            rec = result.get(uid, {})
            caption = (rec.get("caption") or "").strip()
            if caption:
                out[caption.upper()] = rec
    return out


def _fasta_keys(header: str) -> set[str]:
    """Every versionless-uppercase accession key a FASTA header can be found by.

    NCBI emits two header shapes: a bare `>WP_0001.1 desc` and a pipe-delimited
    `>sp|P19833.1|LIP1_MORS1` / `>pdb|5ZOA|A` for Swiss-Prot, PDB and legacy INSD.
    A row's accession may be the middle field (P19833) or, for PDB, entry_chain
    (5ZOA_A), so index the sequence under all plausible keys rather than guessing.
    """
    seqid = header.split()[0] if header.split() else ""
    keys: set[str] = set()
    if "|" in seqid:
        parts = [p for p in seqid.split("|") if p]
        if len(parts) >= 2 and parts[0].isalpha():
            keys.add(parts[1])
            if parts[0].lower() == "pdb" and len(parts) >= 3:
                keys.add(f"{parts[1]}_{parts[2]}")
        else:
            keys.update(parts)
    elif seqid:
        keys.add(seqid)
    return {re.sub(r"\.\d+$", "", k).upper() for k in keys if k}


def efetch_fasta(accessions: list[str], db: str = "protein") -> dict[str, str]:
    """Resolve accessions -> amino-acid sequence via efetch(rettype=fasta).

    Returns a dict keyed by the versionless caption (upper-case), value is the
    bare sequence with all whitespace removed. Ids that do not resolve are simply
    absent from the result, exactly like efetch_docsum.
    """
    out: dict[str, str] = {}
    if not accessions:
        return out
    text = call(
        "efetch",
        [("db", db), ("id", ",".join(accessions)), ("rettype", "fasta"), ("retmode", "text")],
    )
    for block in text.split(">"):
        block = block.strip()
        if not block:
            continue
        header, _, body = block.partition("\n")
        seq = re.sub(r"[^A-Za-z]", "", body).upper()
        if not seq:
            continue
        for key in _fasta_keys(header):
            out[key] = seq
    return out


def efetch_fasta_uids(uids: list[str], db: str = "protein") -> dict[str, str]:
    """Fetch amino-acid sequences by NCBI uid, keyed by uid.

    For records efetch will not take by bare accession -- notably PDB entries
    (7NEI) and chains (5ZOA_A) -- Tier 0 already resolved a uid; fetch by that.
    """
    out: dict[str, str] = {}
    if not uids:
        return out
    text = call(
        "efetch",
        [("db", db), ("id", ",".join(uids)), ("rettype", "fasta"), ("retmode", "text")],
    )
    blocks = [b for b in text.split(">") if b.strip()]
    # efetch preserves input order; zip decoded records back onto the uid list.
    for uid, block in zip(uids, blocks):
        _, _, body = block.partition("\n")
        seq = re.sub(r"[^A-Za-z]", "", body).upper()
        if seq:
            out[str(uid)] = seq
    return out


def efetch_ipg(accession: str, db: str = "protein") -> list[dict]:
    """Identical Protein Group for one accession (efetch rettype=ipg, TSV).

    The IPG lists every accession NCBI considers to encode an identical sequence.
    This is how a suppressed record's live replacement is recovered when
    `replacedby` is empty: find a current accession in the same group.

    Returns a list of {source, protein, name, organism} dicts, one per member.
    """
    text = call("efetch", [("db", db), ("id", accession), ("rettype", "ipg"), ("retmode", "text")])
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if len(lines) < 2:
        return []
    header = [h.strip().lower() for h in lines[0].split("\t")]
    idx = {name: header.index(name) for name in header}
    out = []
    for ln in lines[1:]:
        cols = ln.split("\t")
        def col(name: str) -> str:
            i = idx.get(name)
            return cols[i].strip() if i is not None and i < len(cols) else ""
        prot = col("protein")
        if prot:
            out.append({
                "source": col("source"),
                "protein": prot,
                "name": col("protein name"),
                "organism": col("organism"),
            })
    return out


def efetch_cds_aa(nuc_accession: str) -> list[dict]:
    """CDS translations of a nucleotide record (efetch rettype=fasta_cds_aa).

    Returns one {protein_id, seq} per coding sequence. This is how a protein
    accession is verified to be the CDS product of a *cited nucleotide* accession:
    a paper may report its source as a nucleotide entry (e.g. AY026386) while the
    sheet carries the protein (AAK11724.1); the link holds iff the protein appears
    here as a `protein_id`, or its sequence matches a translation.

    Headers look like:
      >lcl|AY026386.1_prot_AAK11724.1_3 [protein=amidase] [protein_id=AAK11724.1] ...
    """
    text = call(
        "efetch",
        [("db", "nuccore"), ("id", nuc_accession),
         ("rettype", "fasta_cds_aa"), ("retmode", "text")],
    )
    out: list[dict] = []
    for block in text.split(">"):
        block = block.strip()
        if not block:
            continue
        header, _, body = block.partition("\n")
        m = re.search(r"\[protein_id=([^\]]+)\]", header)
        seq = re.sub(r"[^A-Za-z]", "", body).upper()
        if seq:
            out.append({"protein_id": (m.group(1).strip() if m else ""), "seq": seq})
    return out


def efetch_taxonomy(taxids: list[str]) -> dict[str, dict]:
    """Resolve taxids -> lineage, keyed by every taxid that maps to the record.

    Returns dict[taxid] -> {taxid, name, rank, ranks} where `ranks` maps each
    ranked ancestor (plus the record itself) to its taxid, e.g.
    {"species": "1547922", "genus": "1547921", ...}. That is what lets a
    strain-level and a species-level id be compared at a common rank.

    A record is indexed under both its primary TaxId and any AkaTaxIds (taxids
    that were merged into it), so a stale PAZy taxid still resolves to the
    current lineage instead of reading as uncomparable.
    """
    import xml.etree.ElementTree as ET

    out: dict[str, dict] = {}
    clean = sorted({t for t in (str(x).strip() for x in taxids) if t.isdigit()})
    if not clean:
        return out
    for batch in chunked(clean, 200):
        text = call("efetch", [("db", "taxonomy"), ("id", ",".join(batch)), ("retmode", "xml")])
        try:
            root = ET.fromstring(text)
        except ET.ParseError as exc:
            raise RuntimeError(f"taxonomy XML parse failed: {exc}") from exc
        for taxon in root.findall("Taxon"):
            taxid = (taxon.findtext("TaxId") or "").strip()
            if not taxid:
                continue
            rank = (taxon.findtext("Rank") or "").strip().lower()
            name = (taxon.findtext("ScientificName") or "").strip()
            ranks: dict[str, str] = {}
            for anc in taxon.findall("./LineageEx/Taxon"):
                anc_rank = (anc.findtext("Rank") or "").strip().lower()
                anc_taxid = (anc.findtext("TaxId") or "").strip()
                if anc_rank and anc_rank != "no rank" and anc_taxid:
                    ranks.setdefault(anc_rank, anc_taxid)
            if rank and rank != "no rank":
                ranks[rank] = taxid  # the record itself outranks its lineage copy
            rec = {"taxid": taxid, "name": name, "rank": rank, "ranks": ranks}
            out[taxid] = rec
            for aka in taxon.findall("./AkaTaxIds/TaxId"):
                aka_id = (aka.text or "").strip()
                if aka_id:
                    out[aka_id] = rec
    return out


def esearch_accession(acc: str, db: str = "protein") -> list[str]:
    """Resolve a bare accession to uids via the [Accession] index.

    Needed for identifiers efetch will not take directly -- notably bare PDB
    entry ids like 7NEI, which efetch rejects but which index to one uid per chain.
    """
    data = call_json(
        "esearch", [("db", db), ("term", f"{acc}[Accession]"), ("retmax", "10")]
    )
    return data.get("esearchresult", {}).get("idlist", [])


def efetch_docsum_uids(uids: list[str], db: str = "protein") -> dict[str, dict]:
    """Fetch docsums by uid, keyed by uid."""
    out: dict[str, dict] = {}
    if not uids:
        return out
    docs = call_json_all(
        "efetch", [("db", db), ("id", ",".join(uids)), ("rettype", "docsum")]
    )
    for data in docs:
        result = data.get("result", {})
        for uid in result.get("uids", []):
            out[str(uid)] = result.get(uid, {})
    return out


def elink_generic(uids: list[str], db_from: str, db_to: str) -> dict[str, list[str]]:
    """elink in any direction, keyed by input uid. Filters on `<db_from>_<db_to>`.

    Selecting the right linkname matters: a pubmed->protein request also carries
    a `pubmed_pubmed` linkset (related articles) that will happily be mistaken for
    a result if you do not filter.
    """
    out: dict[str, list[str]] = {}
    if not uids:
        return out
    params = [("dbfrom", db_from), ("db", db_to), ("cmd", "neighbor")]
    params += [("id", u) for u in uids]
    linkname = f"{db_from}_{db_to}"
    for data in call_json_all("elink", params, retries=2, timeout=30, reject_errors=True):
        for linkset in data.get("linksets", []):
            ids = linkset.get("ids") or []
            if not ids:
                continue
            links: list[str] = []
            for lsdb in linkset.get("linksetdbs", []):
                if lsdb.get("linkname") == linkname:
                    links.extend(lsdb.get("links", []))
            out[str(ids[0])] = links
    return out


def elink_pubmed_bisect(
    uids: list[str], db_from: str = "protein", failed: set[str] | None = None
) -> dict[str, list[str]]:
    """elink over a batch, halving on failure down to single uids.

    NCBI truncates elink replies for some id sets regardless of batch size (curl
    reports "transfer closed with outstanding read data remaining"). Bisecting
    isolates the offending uid instead of losing the whole batch.

    A uid that fails even alone is added to `failed`. Callers must treat those as
    "unknown", never as "no linked literature" -- an absent key in the returned
    dict would otherwise be indistinguishable from a genuine absence of PubMed
    links, silently turning a network fault into biological evidence.
    """
    if failed is None:
        failed = set()
    try:
        return elink_pubmed(uids, db_from)
    except RuntimeError:
        if len(uids) == 1:
            print(f"  elink: giving up on uid {uids[0]}", file=sys.stderr)
            failed.add(uids[0])
            return {}
    mid = len(uids) // 2
    out = elink_pubmed_bisect(uids[:mid], db_from, failed)
    out.update(elink_pubmed_bisect(uids[mid:], db_from, failed))
    return out


def elink_pubmed(uids: list[str], db_from: str = "protein") -> dict[str, list[str]]:
    """Map each protein uid -> linked PMIDs.

    Passing repeated `id=` params (rather than one comma-joined `id=`) makes NCBI
    return one linkset per input uid instead of a single merged linkset.

    Only ever queries `<db_from> -> pubmed`. Do not pass db_from="pubmed" expecting
    the reverse link: that yields pubmed->pubmed (related articles, capped at 100),
    which silently looks like a rich result. Use elink_generic for other directions.
    """
    out: dict[str, list[str]] = {}
    if not uids:
        return out
    if db_from == "pubmed":
        raise ValueError("elink_pubmed maps X->pubmed; use elink_generic for pubmed->X")
    params = [("dbfrom", db_from), ("db", "pubmed"), ("cmd", "neighbor")]
    params += [("id", u) for u in uids]
    # Fail fast and let the bisect wrapper, not the retry loop, resolve bad id
    # sets: a doomed batch should cost ~20s, not retries * 60s.
    for data in call_json_all("elink", params, retries=1, timeout=20, reject_errors=True):
        for linkset in data.get("linksets", []):
            ids = linkset.get("ids") or []
            if not ids:
                continue
            pmids: list[str] = []
            for lsdb in linkset.get("linksetdbs", []):
                if lsdb.get("linkname") == f"{db_from}_pubmed":
                    pmids.extend(lsdb.get("links", []))
            out[str(ids[0])] = pmids
    return out


def pubmed_dois(pmids: list[str]) -> dict[str, dict]:
    """Map each PMID -> {doi, title, year}."""
    out: dict[str, dict] = {}
    if not pmids:
        return out
    for data in call_json_all("esummary", [("db", "pubmed"), ("id", ",".join(pmids))]):
        result = data.get("result", {})
        for uid in result.get("uids", []):
            rec = result.get(uid, {})
            ids = {a.get("idtype"): a.get("value") for a in rec.get("articleids", [])}
            out[uid] = {
                "doi": (ids.get("doi") or "").strip().lower(),
                "title": rec.get("title", ""),
                "year": (rec.get("pubdate", "") or "")[:4],
            }
    return out


def pubmed_authors(pmids: list[str]) -> dict[str, set[str]]:
    """Map each PMID -> the set of its authors' surnames (lower-case).

    esummary carries `authors: [{name: "Smith J", authtype: "Author"}]`; the
    surname is the leading token. Surnames (not initials) are what Tier A overlaps
    between a candidate accession's linked paper and the row's cited DOI -- a
    shared author is strong evidence the two point at the same protein.
    """
    out: dict[str, set[str]] = {}
    if not pmids:
        return out
    for data in call_json_all("esummary", [("db", "pubmed"), ("id", ",".join(pmids))]):
        result = data.get("result", {})
        for uid in result.get("uids", []):
            rec = result.get(uid, {})
            names: set[str] = set()
            for a in rec.get("authors", []):
                if a.get("authtype") and a["authtype"] != "Author":
                    continue
                surname = (a.get("name") or "").strip().split(" ")[0].lower()
                if surname:
                    names.add(surname)
            out[uid] = names
    return out


def doi_to_pmid(dois: list[str]) -> dict[str, str]:
    """Resolve DOIs -> PMIDs via esearch on the [DOI] field, one call per DOI."""
    out: dict[str, str] = {}
    for doi in dois:
        data = call_json(
            "esearch", [("db", "pubmed"), ("term", f"{doi}[DOI]"), ("retmax", "5")]
        )
        idlist = data.get("esearchresult", {}).get("idlist", [])
        if idlist:
            out[doi] = idlist[0]
    return out
