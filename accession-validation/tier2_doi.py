"""Tier 2 -- DOI cross-reference (strongest signal).

For each resolved accession, walk protein -> pubmed via elink, resolve those PMIDs
to DOIs, and intersect with the DOI set the row cites.

Evidence is gathered in two directions, because either can be missing:
  forward : accession -> linked PMIDs -> their DOIs, intersected with row DOIs
  reverse : row DOIs -> PMIDs (esearch [DOI]), intersected with linked PMIDs

A hit in either direction is a match. Absence is NOT proof of error: many
sequences were deposited without a link to the primary activity paper. Those rows
are "unlinked", and fall through to Tier 1 / manual review.
"""

from __future__ import annotations

import collections

import dataset
import ncbi
import normalize

# elink truncates its reply on large id lists; 40 is comfortably under the limit.
ELINK_BATCH = 40
ESUMMARY_BATCH = 200

FIELDS = [
    "row", "pazy_id", "enzyme_name", "accession", "uid",
    "row_dois", "linked_pmids", "linked_dois", "matched_dois", "matched_pmids",
    "n_linked_pmids", "tier2_flag",
]


def main() -> int:
    rows = {r.index: r for r in dataset.load_rows()}
    tier0 = dataset.read_json("tier0.json")

    # Only rows with a live protein record and a real GI can be elinked. WGS/MAG
    # proteins carry gi=0 and have no pubmed neighbour index.
    linkable = [t for t in tier0 if t["tier0_pass"] and t["uid"] not in ("", "0")]
    print(f"tier2: {len(linkable)}/{len(tier0)} rows are elink-able")

    # --- forward: protein uid -> pmids -------------------------------------
    uid_pmids: dict[str, list[str]] = {}
    elink_failed: set[str] = set()
    uids = sorted({t["uid"] for t in linkable})
    for i, batch in enumerate(ncbi.chunked(uids, ELINK_BATCH), 1):
        print(f"  elink batch {i} ({len(batch)} uids)")
        uid_pmids.update(ncbi.elink_pubmed_bisect(batch, failed=elink_failed))
    if elink_failed:
        print(f"tier2: WARNING -- elink failed for {len(elink_failed)} uids; "
              f"these are reported as elink_failed, not as unlinked")

    # --- pmid -> doi --------------------------------------------------------
    all_pmids = sorted({p for pmids in uid_pmids.values() for p in pmids})
    print(f"tier2: resolving {len(all_pmids)} distinct PMIDs to DOIs")
    pmid_info: dict[str, dict] = {}
    for i, batch in enumerate(ncbi.chunked(all_pmids, ESUMMARY_BATCH), 1):
        print(f"  pubmed esummary batch {i} ({len(batch)} pmids)")
        pmid_info.update(ncbi.pubmed_dois(batch))

    # --- reverse: row doi -> pmid ------------------------------------------
    # Only needed for rows that did not already match forward; each is one call.
    row_dois = {t["row"]: normalize.split_dois(rows[t["row"]].doi) for t in tier0}

    def forward_match(t: dict) -> tuple[set[str], set[str], set[str]]:
        pmids = set(uid_pmids.get(t["uid"], []))
        dois = {pmid_info.get(p, {}).get("doi", "") for p in pmids}
        dois.discard("")
        return pmids, dois, dois & set(row_dois[t["row"]])

    need_reverse: set[str] = set()
    for t in linkable:
        pmids, dois, matched = forward_match(t)
        if not matched and pmids:
            need_reverse.update(row_dois[t["row"]])
    need_reverse = {d for d in need_reverse if d}
    print(f"tier2: reverse-resolving {len(need_reverse)} DOIs to PMIDs")
    doi_pmid = ncbi.doi_to_pmid(sorted(need_reverse))

    # --- adjudicate ---------------------------------------------------------
    out = []
    for t in tier0:
        r = rows[t["row"]]
        dois = row_dois[t["row"]]
        rec = {
            "row": r.index, "pazy_id": r.pazy_id, "enzyme_name": r.enzyme_name,
            "accession": r.accession, "uid": t["uid"], "row_dois": "; ".join(dois),
        }

        if not t["tier0_pass"]:
            rec.update(linked_pmids="", linked_dois="", matched_dois="", matched_pmids="",
                       n_linked_pmids=0, tier2_flag="skipped_tier0")
        elif t["uid"] in ("", "0"):
            rec.update(linked_pmids="", linked_dois="", matched_dois="", matched_pmids="",
                       n_linked_pmids=0, tier2_flag="not_linkable")
        elif t["uid"] in elink_failed:
            # Fetch fault, not a biological finding. Never fold into `unlinked`.
            rec.update(linked_pmids="", linked_dois="", matched_dois="", matched_pmids="",
                       n_linked_pmids=0, tier2_flag="elink_failed")
        else:
            pmids, linked_dois, matched_dois = forward_match(t)
            rev_pmids = {doi_pmid[d] for d in dois if d in doi_pmid}
            matched_pmids = rev_pmids & pmids

            if matched_dois or matched_pmids:
                flag = "doi_match"
            elif not pmids:
                flag = "unlinked"          # no literature link at all
            else:
                # Linked to papers, but not the cited one -- usually a genome
                # announcement (PET5 -> 10.1038/ncomms3156), not a contradiction.
                # Deliberately not called "mismatch": this is silence, not conflict.
                flag = "no_doi_overlap"

            rec.update(
                linked_pmids="; ".join(sorted(pmids)),
                linked_dois="; ".join(sorted(linked_dois)),
                matched_dois="; ".join(sorted(matched_dois)),
                matched_pmids="; ".join(sorted(matched_pmids)),
                n_linked_pmids=len(pmids),
                tier2_flag=flag,
            )
        out.append(rec)

    dataset.write_json("tier2.json", out)
    dataset.write_csv("tier2.csv", FIELDS, out)

    counts = collections.Counter(r["tier2_flag"] for r in out)
    print(f"\ntier2: {len(out)} rows")
    for k, v in counts.most_common():
        print(f"  {k:<16} {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
