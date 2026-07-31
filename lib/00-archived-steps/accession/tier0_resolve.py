"""Tier 0 -- existence and format.

Regex-classifies every accession, then batch-resolves against the NCBI protein
database. Emits resolves / is_protein / is_current / replaced_by per row.
"""

from __future__ import annotations

import sys

import dataset
import ncbi
import normalize

BATCH = 150


def main() -> int:
    rows = dataset.load_rows()
    print(f"tier0: {len(rows)} rows")

    # --- syntax pass -------------------------------------------------------
    for r in rows:
        cls, versioned = normalize.classify_accession(r.accession)
        r.flags["accession_class"] = cls
        r.flags["has_version"] = versioned

    lookupable = [r for r in rows if r.flags["accession_class"] in normalize.NCBI_PROTEIN_CLASSES]
    n_skipped = len(rows) - len(lookupable)
    print(f"tier0: {len(lookupable)} look up in NCBI protein; {n_skipped} skipped by class")

    # --- resolve pass ------------------------------------------------------
    # efetch drops unresolvable ids silently, so a caption-keyed dict is the
    # source of truth for what came back.
    records: dict[str, dict] = {}
    accessions = sorted({r.accession for r in lookupable})
    for i, batch in enumerate(ncbi.chunked(accessions, BATCH), 1):
        print(f"  efetch batch {i} ({len(batch)} ids)")
        records.update(ncbi.efetch_docsum(batch))

    def still_missing() -> list[str]:
        return [a for a in accessions if normalize.strip_version(a).upper() not in records]

    # Fallback 1: solo efetch, to separate a genuinely dead accession from a
    # batch-level hiccup.
    missing = still_missing()
    if missing:
        print(f"tier0: {len(missing)} unresolved in batch; retrying individually")
        for acc in missing:
            try:
                records.update(ncbi.efetch_docsum([acc]))
            except RuntimeError as exc:
                print(f"  {acc}: {exc}", file=sys.stderr)

    # Fallback 2: esearch on the [Accession] index. efetch rejects bare PDB entry
    # ids (7NEI) that esearch happily maps to one uid per chain.
    nuc_hits: dict[str, dict] = {}
    missing = still_missing()
    if missing:
        print(f"tier0: {len(missing)} still unresolved; trying esearch[Accession]")
        for acc in missing:
            uids = ncbi.esearch_accession(acc, db="protein")
            if not uids:
                continue
            got = ncbi.efetch_docsum_uids(uids[:1], db="protein")
            for rec in got.values():
                records[normalize.strip_version(acc).upper()] = rec

    # Fallback 3: is it a *nucleotide* accession sitting in a protein column?
    # This is a distinct, more actionable failure than "does not exist".
    missing = still_missing()
    if missing:
        print(f"tier0: {len(missing)} unresolved in protein; probing nuccore")
        for acc in missing:
            try:
                got = ncbi.efetch_docsum([acc], db="nuccore")
            except RuntimeError:
                continue
            key = normalize.strip_version(acc).upper()
            if key in got:
                nuc_hits[key] = got[key]

    # --- annotate ----------------------------------------------------------
    out = []
    for r in rows:
        key = normalize.strip_version(r.accession).upper()
        rec = records.get(key)
        cls = r.flags["accession_class"]
        note = ""

        if rec is None and key in nuc_hits:
            # Resolves, but to a nucleotide record -- wrong molecule type.
            nrec = nuc_hits[key]
            r.flags.update(
                resolves=True, is_protein=False, is_current=True, replaced_by="",
                uid=str(nrec.get("gi") or ""), title=nrec.get("title", ""),
                organism=nrec.get("organism", ""), taxid=str(nrec.get("taxid") or ""),
                slen=nrec.get("slen", ""),
                sourcedb="nuccore", status=(nrec.get("status") or "").strip(),
            )
            note = f"nucleotide record ({nrec.get('moltype', '?')}/{nrec.get('biomol', '?')}), not a protein"
        elif rec is None:
            r.flags.update(
                resolves=False, is_protein=False, is_current=False, replaced_by="",
                uid="", title="", organism="", taxid="", slen="", sourcedb="", status="",
            )
            if cls == "mgnify":
                note = "MGnify identifier; not an NCBI record (resolve via MGnify API)"
            elif cls == "uniprot":
                note = "UniProt/TrEMBL identifier absent from NCBI protein (resolve via UniProt API)"
            else:
                note = "not found in NCBI protein or nuccore"
        else:
            status = (rec.get("status") or "").strip()
            replaced_by = (rec.get("replacedby") or "").strip()
            moltype = (rec.get("moltype") or "").lower()
            biomol = (rec.get("biomol") or "").lower()
            is_protein = moltype == "aa" or biomol == "peptide"
            r.flags.update(
                resolves=True,
                is_protein=is_protein,
                # NCBI leaves `status` empty for live records; it is populated
                # with "replaced"/"suppressed"/"dead" otherwise.
                is_current=status in ("", "live") and not replaced_by,
                replaced_by=replaced_by,
                uid=str(rec.get("gi") or rec.get("uid") or ""),
                title=rec.get("title", ""),
                organism=rec.get("organism", ""),
                taxid=str(rec.get("taxid") or ""),
                slen=rec.get("slen", ""),
                sourcedb=rec.get("sourcedb", ""),
                status=status,
            )
            if status == "suppressed":
                note = "NCBI record is suppressed"
            elif status == "dead":
                note = "NCBI record is dead"

        r.flags["resolution_note"] = note
        r.flags["accession_version_resolved"] = (rec or {}).get("accessionversion", "")
        r.flags["tier0_pass"] = bool(
            r.flags["resolves"] and r.flags["is_protein"] and r.flags["is_current"]
        )
        out.append({
            "row": r.index, "pazy_id": r.pazy_id, "enzyme_name": r.enzyme_name,
            "accession": r.accession, "accession_class": cls, **r.flags,
        })

    dataset.write_json("tier0.json", out)
    dataset.write_csv(
        "tier0.csv",
        ["row", "pazy_id", "enzyme_name", "accession", "accession_class", "has_version",
         "resolves", "is_protein", "is_current", "replaced_by", "tier0_pass", "resolution_note",
         "accession_version_resolved", "uid", "title", "organism", "taxid", "slen", "sourcedb", "status"],
        out,
    )

    # --- summary -----------------------------------------------------------
    n = len(out)
    passed = sum(1 for r in out if r["tier0_pass"])
    print(f"\ntier0 summary: {passed}/{n} pass")
    for label, pred in [
        ("unresolvable", lambda r: not r["resolves"]),
        ("resolves to nucleotide record", lambda r: r["resolves"] and not r["is_protein"]),
        ("suppressed/dead", lambda r: r["resolves"] and not r["is_current"]),
        ("versionless accession (cosmetic)", lambda r: not r["has_version"]),
    ]:
        hits = [r for r in out if pred(r)]
        if hits:
            print(f"  {label}: {len(hits)}")
            for r in hits[:6]:
                print(f"     {r['accession']:<20} {r['enzyme_name']:<18} {r['resolution_note'][:46]}")
            if len(hits) > 6:
                print(f"     ... and {len(hits) - 6} more")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
