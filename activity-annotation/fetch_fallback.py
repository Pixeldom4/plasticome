#!/usr/bin/env python3
"""Fetch the DOIs that fetch_doi.py could not resolve via publisher APIs.

ACS and RSC expose no CrossRef TDM links, and the scrapy-bot prefix map routes
10.1128 (ASM) to Wiley. Each DOI below is given an explicit route, tried in
order, preferring full-text XML over PDF for downstream extraction:

  epmc   Europe PMC fullTextXML   (OA-licensed articles only)
  efetch NCBI E-utilities         (also serves author manuscripts; skipped when
                                   the response carries no <body>, which means
                                   metadata-only -- e.g. old scanned articles)
  pdf    direct publisher PDF

RSC sits behind Cloudflare: only curl_cffi's plain "chrome" profile clears it;
chrome120/124/131 all get a 403 challenge page.
"""
import os
import re
import sys
import tempfile
import time

from curl_cffi import requests as cr

EMAIL = "dennis.zhu@mail.utoronto.ca"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "papers")

EPMC_XML = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
EFETCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
          "?db=pmc&id={num}&retmode=xml")

# doi -> ordered list of (route, argument) attempts
TARGETS = {
    "10.1021/acs.est.8b04252": [
        ("epmc", "PMC12447631"),
        ("efetch", "12447631"),
    ],
    "10.1021/acscatal.2c02275": [
        ("epmc", "PMC9361285"),
    ],
    "10.1021/acscatal.5c03460": [
        ("epmc", "PMC12455559"),
        ("pdf", "https://pubs.acs.org/doi/pdf/10.1021/acscatal.5c03460?ref=article_openPDF"),
    ],
    "10.1039/D4SU00722K": [
        ("pdf", "https://pubs.rsc.org/en/content/articlepdf/2025/su/d4su00722k"),
    ],
    "10.1128/aem.65.1.189-197.1999": [
        # PMC91002 is metadata-only; the ASM PDF is the sole full text.
        ("efetch", "91002"),
        ("pdf", "https://journals.asm.org/doi/pdf/10.1128/aem.65.1.189-197.1999"),
    ],
}

# Hosts needing a homepage visit + the plain "chrome" profile to clear Cloudflare
PRIME = {"pubs.rsc.org": "https://pubs.rsc.org/"}


def doi_safe(doi):
    return re.sub(r'[<>:"|?*\\]', "_", doi.strip().lower().replace("/", "_"))


def atomic_write(path, data):
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path), suffix=".tmp")
    try:
        os.write(fd, data)
        os.close(fd)
        os.replace(tmp, path)
    except Exception:
        os.unlink(tmp)
        raise


def existing(doi):
    for ext in (".xml", ".pdf"):
        p = os.path.join(OUT, doi_safe(doi) + ext)
        if os.path.exists(p) and os.path.getsize(p) > 10_000:
            return p
    return None


def save_xml(session, doi, url, require_body):
    r = session.get(url, timeout=60)
    if r.status_code != 200 or b"<article" not in r.content[:4000]:
        return None
    if require_body and b"<body" not in r.content:
        return None  # metadata-only record
    path = os.path.join(OUT, doi_safe(doi) + ".xml")
    atomic_write(path, r.content)
    return path


def save_pdf(session, doi, url):
    host = url.split("/")[2]
    if host in PRIME:
        session = cr.Session(impersonate="chrome")
        session.headers.update({"From": EMAIL})
        session.get(PRIME[host], timeout=45, allow_redirects=True)
    r = session.get(url, timeout=90, allow_redirects=True,
                    headers={"Referer": "https://doi.org/" + doi})
    if r.status_code != 200 or r.content[:4] != b"%PDF" or len(r.content) < 10_000:
        return None
    path = os.path.join(OUT, doi_safe(doi) + ".pdf")
    atomic_write(path, r.content)
    return path


def main():
    os.makedirs(OUT, exist_ok=True)
    session = cr.Session(impersonate="chrome120")
    session.headers.update({"From": EMAIL})

    ok = fail = 0
    for doi, attempts in TARGETS.items():
        if (p := existing(doi)):
            ok += 1
            print(f"  SKIP  {doi} -> {os.path.basename(p)}")
            continue

        path = None
        for route, arg in attempts:
            try:
                if route == "epmc":
                    path = save_xml(session, doi, EPMC_XML.format(pmcid=arg), False)
                elif route == "efetch":
                    time.sleep(0.5)  # NCBI rate limit without an API key
                    path = save_xml(session, doi, EFETCH.format(num=arg), True)
                elif route == "pdf":
                    path = save_pdf(session, doi, arg)
            except Exception as exc:
                print(f"        {route} error: {type(exc).__name__}", file=sys.stderr)
            if path:
                break

        if path:
            ok += 1
            print(f"  OK    {doi} -> {os.path.basename(path)}  [{route}]")
        else:
            fail += 1
            print(f"  FAIL  {doi}")

    print(f"\nDone: {ok} ok/skip, {fail} failed of {len(TARGETS)}")


if __name__ == "__main__":
    main()
