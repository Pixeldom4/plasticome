"""Accession classification, name normalization, and concordance scoring."""

from __future__ import annotations

import difflib
import re

# ---------------------------------------------------------------------------
# Tier 0: accession syntax
# ---------------------------------------------------------------------------

ACCESSION_PATTERNS: list[tuple[str, str]] = [
    # RefSeq protein
    ("refseq", r"^(?:WP|NP|XP|YP|AP|ZP)_\d{6,9}(?:\.\d+)?$"),
    # RefSeq nucleotide prefixes -- wrong molecule type for a protein sheet
    ("refseq_nucleotide", r"^(?:NC|NM|NR|XM|XR|NZ|NT|NW)_\d+(?:\.\d+)?$"),
    # GenBank/EMBL/DDBJ protein: 3 letters + 5 digits
    ("insd_protein", r"^[A-Z]{3}\d{5}(?:\.\d+)?$"),
    # WGS / MAG / TSA derived proteins: 3 letters + 7 digits
    ("wgs_protein", r"^[A-Z]{3}\d{7}(?:\.\d+)?$"),
    # UniProt/Swiss-Prot
    ("uniprot", r"^(?:[OPQ]\d[A-Z0-9]{3}\d|[A-NR-Z]\d(?:[A-Z][A-Z0-9]{2}\d){1,2})(?:\.\d+)?$"),
    # PDB entry with chain (5ZOA_A) or bare entry id (7NEI)
    ("pdb_chain", r"^[0-9][A-Z0-9]{3}_[A-Z0-9]+$"),
    ("pdb_entry", r"^[0-9][A-Z0-9]{3}$"),
    # MGnify protein -- not an NCBI identifier at all
    ("mgnify", r"^MGYP\d+$"),
    # Legacy GenBank nucleotide-style: 1 letter + 5 digits / 2 letters + 6 digits
    ("insd_legacy", r"^(?:[A-Z]\d{5}|[A-Z]{2}\d{6})(?:\.\d+)?$"),
]

# Classes we expect to resolve in the NCBI protein database.
NCBI_PROTEIN_CLASSES = {
    "refseq", "insd_protein", "wgs_protein", "uniprot",
    "pdb_chain", "pdb_entry", "insd_legacy",
}


def classify_accession(acc: str) -> tuple[str, bool]:
    """Return (class_name, has_version)."""
    acc = (acc or "").strip()
    if not acc or acc == "-":
        return "missing", False
    has_version = bool(re.search(r"\.\d+$", acc))
    for name, pattern in ACCESSION_PATTERNS:
        if re.match(pattern, acc, re.I):
            return name, has_version
    return "unrecognized", has_version


def strip_version(acc: str) -> str:
    return re.sub(r"\.\d+$", "", (acc or "").strip())


# ---------------------------------------------------------------------------
# Tier 1: name normalization
# ---------------------------------------------------------------------------

# Titles that describe no function at all -- resolving is not confirming.
UNINFORMATIVE_TITLE = re.compile(
    r"^(?:hypothetical|uncharacterized|unnamed|predicted|putative uncharacterized)\b"
    r"|^protein of unknown function"
    r"|^(?:conserved )?(?:hypothetical )?protein\s+\S+$",
    re.I,
)

# Canonical enzyme families. Each maps surface forms -> a single token.
SYNONYMS: dict[str, list[str]] = {
    "pet_hydrolase": [
        "poly(ethylene terephthalate) hydrolase", "pet hydrolase", "petase",
        "poly(ethylene terephthalate) hydrolase family protein",
        "bta-hydrolase", "bta hydrolase", "polyester hydrolase",
    ],
    "cutinase": ["cutinase", "cutin hydrolase", "cutinase-like", "cutinase family protein"],
    "esterase": [
        "esterase", "carboxylesterase", "carboxylic ester hydrolase",
        "arylesterase", "acetyl esterase", "esterase family protein",
        "esterase/lipase", "alpha/beta fold hydrolase",
    ],
    "lipase": ["lipase", "triacylglycerol lipase", "triacylglycerol acylhydrolase", "lipolytic enzyme"],
    "protease": [
        "protease", "peptidase", "proteinase", "serine protease", "serine peptidase",
        "endopeptidase", "metalloprotease", "metallopeptidase", "subtilisin",
    ],
    "amidase": ["amidase", "amido hydrolase", "amidohydrolase"],
    "urethanase": ["urethanase", "urethane hydrolase", "polyurethanase"],
    "laccase": ["laccase", "multicopper oxidase", "multicopper polyphenol oxidase"],
    "peroxidase": ["peroxidase", "manganese peroxidase", "lignin peroxidase", "versatile peroxidase", "dye-decolorizing peroxidase"],
    "oxidase": ["oxidase", "monooxygenase", "oxygenase", "alkane monooxygenase"],
    "dehydrogenase": ["dehydrogenase", "alcohol dehydrogenase", "aldehyde dehydrogenase"],
    "dienelactone_hydrolase": ["dienelactone hydrolase", "dienelactone hydrolase family protein"],
    "depolymerase": ["depolymerase", "phb depolymerase", "polyhydroxybutyrate depolymerase", "pha depolymerase"],
    "nitrilase": ["nitrilase", "nitrile hydratase"],
    "tannase": ["tannase", "tannin acylhydrolase", "feruloyl esterase"],
    "ab_hydrolase": [
        "alpha/beta hydrolase", "alpha/beta-hydrolase", "ab hydrolase",
        "alpha/beta hydrolase fold protein", "alpha/beta hydrolase family protein",
    ],
    "mhetase": ["mhetase", "mono(2-hydroxyethyl) terephthalate hydrolase", "mhet hydrolase"],
    "hydrolase": ["hydrolase", "hydrolase family protein"],  # fallback, weakest
}

# Longest-first so "pet hydrolase" wins over bare "hydrolase".
_SYNONYM_INDEX: list[tuple[str, str]] = sorted(
    ((surface, canon) for canon, forms in SYNONYMS.items() for surface in forms),
    key=lambda kv: -len(kv[0]),
)

# Affix hints for opaque lab codes (PpCutA, EstB, PE-H, LCC). Applied to the
# raw enzyme_name only -- these are naming conventions, not English words.
CODE_HINTS: list[tuple[str, str]] = [
    (r"petase", "pet_hydrolase"),
    (r"mhetase", "mhetase"),
    (r"\bbta[-_ ]?\d*\b", "pet_hydrolase"),
    (r"\blcc\b|leaf.?compost|leaf.?branch", "cutinase"),
    (r"cut(?:in)?(?:ase)?[a-z]?\d*\b|\bcut\d+", "cutinase"),
    (r"\best(?:erase)?[a-z]?\d*\b", "esterase"),
    (r"lip(?:ase)?[a-z]?\d*\b", "lipase"),
    (r"\bpe[-_]?h\b", "pet_hydrolase"),
    (r"protease|peptidase|subtilisin", "protease"),
    (r"laccase", "laccase"),
    (r"peroxidase", "peroxidase"),
    (r"amidase", "amidase"),
    (r"urethanase", "urethanase"),
    (r"depolymerase", "depolymerase"),
    (r"tannase", "tannase"),
    (r"nitrilase", "nitrilase"),
    (r"\bpet\d*\b|\bpet[a-z]?\b", "pet_hydrolase"),
]

# Families whose annotations drift freely into one another. All are serine
# hydrolases of the same structural clan: a PETase is an alpha/beta hydrolase, and
# cutinase/esterase/lipase are routinely used interchangeably by NCBI annotators.
# Disagreement *within* the clan is weak evidence; disagreement *across* clans
# (lipase vs laccase) is the real wrong-protein signal.
SERINE_HYDROLASE_CLAN = {
    "pet_hydrolase", "cutinase", "esterase", "lipase", "mhetase",
    "dienelactone_hydrolase", "ab_hydrolase", "tannase", "urethanase",
    "depolymerase", "amidase", "hydrolase",
}

# Superfamily-level labels that are compatible with any clan member by definition.
GENERIC_FAMILIES = {"ab_hydrolase", "hydrolase"}

CLANS: list[set[str]] = [
    SERINE_HYDROLASE_CLAN,
    {"protease"},                                   # peptidases
    {"laccase", "peroxidase", "oxidase", "dehydrogenase"},  # oxidoreductases
    {"nitrilase"},
]


def _same_clan(a: str, b: str) -> bool:
    return any(a in clan and b in clan for clan in CLANS)


STOPWORDS = {
    "protein", "family", "putative", "probable", "chain", "precursor", "partial",
    "multispecies", "mag", "recname", "full", "flags", "the", "of", "and", "a",
    "isoform", "variant", "fragment", "domain", "containing", "type", "class",
    "subunit", "enzyme", "gene", "product", "unknown", "function", "uncharacterized",
}


def clean_title(title: str) -> tuple[str, str]:
    """Strip NCBI title decorations. Returns (clean_title, organism)."""
    t = (title or "").strip()
    organism = ""
    m = re.search(r"\[([^\[\]]+)\]\s*$", t)
    if m:
        organism = m.group(1).strip()
        t = t[: m.start()].strip()
    t = re.sub(r"^MULTISPECIES:\s*", "", t, flags=re.I)
    t = re.sub(r"^MAG(?:\s+TPA)?:\s*", "", t, flags=re.I)
    t = re.sub(r"^Chain\s+[A-Z0-9]+,\s*", "", t, flags=re.I)
    t = re.sub(r"^RecName:\s*Full=", "", t, flags=re.I)
    t = re.sub(r";\s*(?:AltName|Flags|Short)\b.*$", "", t, flags=re.I)
    t = re.sub(r"\s*;\s*$", "", t)
    return t.strip(), organism


def _canonical_families(text: str) -> set[str]:
    """Extract canonical family tokens from free text via synonym matching."""
    low = " " + re.sub(r"[^a-z0-9()/\- ]+", " ", (text or "").lower()) + " "
    found: set[str] = set()
    for surface, canon in _SYNONYM_INDEX:
        if surface in low:
            found.add(canon)
    # `hydrolase` is subsumed by any more specific *_hydrolase / named family.
    if len(found) > 1:
        found.discard("hydrolase")
    return found


def families_from_title(title: str) -> set[str]:
    return _canonical_families(title)


def families_from_enzyme_name(name: str) -> set[str]:
    """Families implied by an enzyme_name, including opaque-code affix hints."""
    fams = _canonical_families(name)
    low = (name or "").lower()
    for pattern, canon in CODE_HINTS:
        if re.search(pattern, low):
            fams.add(canon)
    if len(fams) > 1:
        fams.discard("hydrolase")
    return fams


def tokens(text: str) -> set[str]:
    raw = re.split(r"[^a-z0-9]+", (text or "").lower())
    return {t for t in raw if t and t not in STOPWORDS and not t.isdigit()}


def fuzzy_ratio(a: str, b: str) -> float:
    return difflib.SequenceMatcher(None, (a or "").lower(), (b or "").lower()).ratio()


def jaccard(a: set, b: set) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def score_concordance(enzyme_name: str, title: str) -> dict:
    """Compare an enzyme_name against an NCBI record title.

    Returns independent signals, not a verdict. `name_informative` is the key
    field: most PAZy enzyme_names are opaque lab codes that imply no family, in
    which case Tier 1 simply cannot speak and Tier 2 must decide.
    """
    clean, organism = clean_title(title)
    name_fams = families_from_enzyme_name(enzyme_name)
    title_fams = families_from_title(clean)

    # Two very different reasons a title yields no family, which must not be
    # conflated:
    #   unannotated -- "hypothetical protein XYZ_001": nothing to compare against.
    #   unmapped    -- "glutamyl-tRNA amidotransferase": a real function that our
    #                  lexicon does not cover. Often the actual wrong protein.
    uninformative_title = bool(UNINFORMATIVE_TITLE.search(clean))
    title_unmapped = not title_fams and not uninformative_title and bool(clean)
    name_informative = bool(name_fams)

    shared = name_fams & title_fams
    if not (name_fams and title_fams):
        family_agreement = "unknown"
    elif shared:
        family_agreement = "match"
    elif title_fams & GENERIC_FAMILIES or name_fams & GENERIC_FAMILIES:
        # e.g. name=PETase, title="alpha/beta hydrolase" -- a superfamily label
        # that is true but uninformative. Consistent, not confirmatory.
        family_agreement = "generic"
    elif any(_same_clan(a, b) for a in name_fams for b in title_fams):
        # e.g. name=cutinase, title="PET hydrolase" -- routine annotation drift.
        family_agreement = "related"
    else:
        # e.g. name=lipase, title=laccase -- cross-clan, likely wrong protein.
        family_agreement = "conflict"

    return {
        "clean_title": clean,
        "organism": organism,
        "name_families": sorted(name_fams),
        "title_families": sorted(title_fams),
        "shared_families": sorted(shared),
        "family_agreement": family_agreement,
        "name_informative": name_informative,
        "uninformative_title": uninformative_title,
        "title_unmapped": title_unmapped,
        "token_jaccard": round(jaccard(tokens(enzyme_name), tokens(clean)), 3),
        "fuzzy_ratio": round(fuzzy_ratio(enzyme_name, clean), 3),
    }


# ---------------------------------------------------------------------------
# DOI handling
# ---------------------------------------------------------------------------

def normalize_doi(doi: str) -> str:
    d = (doi or "").strip().lower()
    d = re.sub(r"^(?:https?://)?(?:dx\.)?doi\.org/", "", d)
    d = re.sub(r"^doi:\s*", "", d)
    return d.strip().rstrip(".;,")


def split_dois(cell: str) -> list[str]:
    parts = re.split(r"[;,]\s*|\s+(?=10\.)", cell or "")
    out = [normalize_doi(p) for p in parts]
    return [d for d in out if d.startswith("10.")]
