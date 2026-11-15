"""One definition of every fold: colour, short name, mechanism, structure.

The tree, the heatmap, the minimum spanning tree and the structure cards all
import this, so a fold has the same colour and the same name everywhere.

Folds are keyed by the `fold_key` of tree_tips.tsv (the catalytic CATH
superfamily), never by the slice id gNN. The slice ids are numbered by size
and differ between the all-sequences and the centroid runs.

Colour. The ten hues of palette/fold_palette.tsv (step 15, the PETadex atlas
hues) are kept as they are. The six folds that are a single sequence used to
share one grey; the reviewers asked for them to be shown, so each gets a hue
of its own, placed in the gaps the ten leave, at lower saturation so they still
read as the quieter end of the palette. Every fold has four shades:

  tint    HSL(h, S, 89)   slice backgrounds, as in the pie
  band    HSL(h, S+8, 72) rings and margin strips, where a tint is too pale
  label   HSL(h, 40, 50)  text set in the fold's colour
  strong  HSL(h, 60, 38)  the atlas base colour, for lines and outlines
"""

import colorsys
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ANNOT_DIR = os.path.join(HERE, "..", "structures", "annotations")
ANNOT_CSV = os.path.join(ANNOT_DIR, "petadex_domain_summary.csv")
PDB_DIR = os.path.join(ANNOT_DIR, "good_folds")

UNCONFIRMED = "3.40.710.10|unconfirmed"

# key: (hue, saturation, short name, mechanism class)
FOLDS = {
    "3.40.50.1820": (210, 52, "Alpha/beta hydrolase", "Ser-His-Asp triad"),
    "3.60.70.12": (6, 52, "Ntn amidase (NylC)", "N-terminal Thr nucleophile"),
    "3.90.1300.10": (0, 52, "Amidase signature", "Ser-cisSer-Lys triad"),
    "domain:oxygenase MpaB family protein": (90, 52, "Oxygenase MpaB (Lcp)", "heme b oxygenase"),
    "3.40.50.200": (258, 52, "Subtilisin S8/S53", "Asp-His-Ser triad"),
    "2.40.10.10": (270, 52, "Trypsin-like protease", "His-Asp-Ser triad"),
    "3.40.710.10": (140, 52, "Beta-lactamase / DD-peptidase", "Ser-Lys-Tyr triad"),
    "1.10.760.10": (315, 52, "Cytochrome c (RoxA)", "diheme c oxygenase"),
    "2.60.40.420": (50, 52, "Cupredoxin (laccase)", "multicopper oxidase"),
    "3.40.50.1110": (240, 52, "SGNH hydrolase", "Ser-Asp-His, shared dimer site"),
    # single sequence folds, new hues in the gaps, lower saturation. The
    # glycoside hydrolase was hue 157, CIELAB distance 12 from the aqua of
    # the modelled plastic; moved to cobalt 225 (distance 70 from the aqua,
    # 16 from its nearest fold, SGNH), at full saturation so its tint stays
    # apart from the alpha/beta hydrolase. The free gaps left were brown or
    # cyan.
    "1.10.420.10;1.10.520.10": (345, 34, "Peroxidase (MnP)", "heme + Mn(II)"),
    "1.10.630.10": (70, 34, "Cytochrome P450", "heme thiolate"),
    "2.40.70.10": (115, 34, "Acid protease", "Asp-Asp dyad"),
    "3.20.20.370": (225, 52, "Glycoside hydrolase/deacetylase", "Glu/His acid-base"),
    "3.60.21.10": (292, 34, "Purple acid phosphatase", "Fe-Zn binuclear"),
    UNCONFIRMED: (140, 14, "Beta-lactamase, unconfirmed", "no structure"),
}

# The fold types drawn in the ring tree (Figure 2a): what each enzyme does,
# in the way of the EC classes (esterases, peptidases, other C-N hydrolases,
# oxygenases, other oxidoreductases). Hydrolases first, then oxidoreductases,
# the larger type first in each; inside a type the larger fold first. The class
# of each fold is read from the name of its representative enzyme, not looked
# up per record.
TYPES = [
    ("Esterases", ["3.40.50.1820", "3.40.50.1110", "3.60.21.10"]),
    ("Amidases and other C-N hydrolases",
     ["3.60.70.12", "3.90.1300.10", "3.40.710.10", UNCONFIRMED, "3.20.20.370"]),
    ("Peptidases", ["3.40.50.200", "2.40.10.10", "2.40.70.10"]),
    ("Oxygenases", ["domain:oxygenase MpaB family protein", "1.10.760.10", "1.10.630.10"]),
    ("Other oxidoreductases", ["1.10.420.10;1.10.520.10", "2.60.40.420"]),
]

# The same folds grouped by catalytic machinery, read from FOLDS[k][3] and the
# mechanism column of fold_representatives.tsv. Every fold sits in one group;
# a group is the chemistry that cuts the bond, whatever the substrate. SGNH
# (Ser-Asp-His) is kept with the triad hydrolases, laccase with the metal
# enzymes because its copper does the work, not a heme.
MECHANISM_TYPES = [
    ("Ser-His-Asp triad", ["3.40.50.1820", "3.40.50.200", "2.40.10.10", "3.40.50.1110"]),
    ("Ser-Lys type serine hydrolases", ["3.90.1300.10", "3.40.710.10", UNCONFIRMED]),
    ("N-terminal Thr", ["3.60.70.12"]),
    ("Acid-base and aspartyl", ["2.40.70.10", "3.20.20.370"]),
    ("Metal-dependent", ["3.60.21.10", "2.60.40.420"]),
    ("Heme enzymes", ["domain:oxygenase MpaB family protein", "1.10.760.10", "1.10.630.10",
                      "1.10.420.10;1.10.520.10"]),
]
GROUPINGS = {"activity": TYPES, "mechanism": MECHANISM_TYPES}


GROUPING = "activity"


def use_grouping(which):
    """Make type_rank and type_name follow one of GROUPINGS; with "mechanism"
    the colours are the group's (tint, band, label, strong) and one fold per
    group is pictured (GROUP_PICTURED)."""
    global TYPES, GROUPING
    TYPES, GROUPING = GROUPINGS[which], which


# Mechanism groups: (hue, saturation, label lightness), one distinct colour per
# group, in MECHANISM_TYPES order. None brown, none near the aqua of the plastic,
# the pink of the catalytic residues or the amber of the secondary ones. Folds
# of one group share the hue and alternate in lightness so neighbours stay apart.
GROUP_STYLE = [(215, 52, 50), (355, 52, 50), (20, 85, 55), (130, 45, 42), (258, 45, 52), (295, 40, 46)]

# The fold pictured for each group: the largest fold in it with an annotated structure.
GROUP_PICTURED = {"3.40.50.1820", "3.90.1300.10", "3.60.70.12", "2.40.70.10", "2.60.40.420",
                  "domain:oxygenase MpaB family protein"}


def group_rank(k):
    return type_rank(k) if GROUPING == "mechanism" else None


def group_dark(rank):
    """The colour of a group's bar and name."""
    h, s, _ = GROUP_STYLE[rank]
    return hsl(h, min(s + 10, 85), 46 if 10 < h < 70 else 34)    # a dark orange reads brown


def _in_group(k):
    return GROUPING == "mechanism" and type_rank(k) < len(GROUP_STYLE)


def _alt(k):
    """0 or 1, alternating along a group's folds."""
    return TYPES[type_rank(k)][1].index(k) % 2


def tint(k):
    if not _in_group(k):
        return _fold_tint(k)
    h, s, _ = GROUP_STYLE[type_rank(k)]
    return hsl(h, s, 89 + 4 * _alt(k))


def band(k):
    if not _in_group(k):
        return _fold_band(k)
    h, s, _ = GROUP_STYLE[type_rank(k)]
    return hsl(h, min(s + 8, 80), 72 + 6 * _alt(k))


def label(k):
    if not _in_group(k):
        return _fold_label(k)
    h, s, l = GROUP_STYLE[type_rank(k)]
    return hsl(h, s, l)


def strong(k):
    if not _in_group(k):
        return _fold_strong(k)
    return group_dark(type_rank(k))


def type_rank(k):
    """Position of a fold's type in TYPES (unknown folds last)."""
    for i, (_, keys) in enumerate(TYPES):
        if k in keys:
            return i
    return len(TYPES)


def type_name(k):
    return TYPES[type_rank(k)][0] if type_rank(k) < len(TYPES) else "Other"


# What these folds do in nature, before anyone tested them on plastic. Not
# verified, so not drawn in any figure (removed 2026-10-27 at the user's
# request); kept here only as a starting list.
NATIVE = {
    "3.40.50.1820": "lipids, cutin, aromatic esters",
    "3.60.70.12": "peptides, beta-aminopeptides",
    "3.90.1300.10": "amides (amidases)",
    "domain:oxygenase MpaB family protein": "polyisoprene (rubber)",
    "3.40.50.200": "proteins",
    "2.40.10.10": "proteins",
    "3.40.710.10": "peptidoglycan, beta-lactams",
    "1.10.760.10": "polyisoprene (rubber)",
    "2.60.40.420": "lignin, phenolics",
    "3.40.50.1110": "lipids, acyl esters",
    "1.10.420.10;1.10.520.10": "lignin",
    "1.10.630.10": "hydrocarbons, fatty acids",
    "2.40.70.10": "proteins",
    "3.20.20.370": "polysaccharide acetyl groups",
    "3.60.21.10": "phosphate esters",
    UNCONFIRMED: "unknown",
}


# The representative structure of each fold: the protein each co-folded model
# in ../structures/annotations/good_folds/ was built from. Named from
# the model's accession in the annotation sheet (NCBI protein or UniProt
# record), with the experimental structure the sheet gives as its reference.
# `record` is the plasticome_id of the same protein in the union (matched by
# accession, or by sequence: 92 to 100 percent of 5-mers shared), or None when
# the model is a reference protein from outside the 609 records.
# Each model was aligned to its source sequence (representatives.faa, from
# NCBI and UniProt on 2026-10-27; the union sequence where the sheet's
# accession was wrong). Step 17 trims terminal stretches
# that are not in the source protein, such as the expression tag on the
# IsPETase model.
# key: (enzyme, organism, accession, reference structure, record)
REPRESENTATIVE = {
    "3.40.50.1820": ("IsPETase", "Piscinibacter sakaiensis", "GAP38373.1", "PDB 5XJH", "1"),
    "3.60.70.12": ("NylC, 6-aminohexanoate-oligomer endohydrolase", "Flavobacterium sp.",
                   "BAA01528.1", "PDB 5XYG", "542"),
    "3.90.1300.10": ("NylA, 6-aminohexanoate-cyclic-dimer hydrolase", "Flavobacterium sp. KI723T1",
                     "BAA05090.1", "PDB 3A2P", "548"),
    "domain:oxygenase MpaB family protein": ("Lcp, latex clearing protein", "Streptomyces sp. K30",
                                             "AAR25849.1", "PDB 5O1L", "106"),
    "3.40.50.200": ("subtilisin (precursor)", "Bacillus amyloliquefaciens", "CAA24990.1",
                    "PDB 1S01", None),
    "2.40.10.10": ("trypsin (preprotrypsin)", "Fusarium oxysporum", "AAB27568.1", "PDB 1PQ5", None),
    # the sheet gives NylB (CAA24927.1), but the model is 98% identical to
    # NylB' and 89% to NylB; its 8 differences from NylB' fit the NylB/NylB'
    # hybrids crystallised for this enzyme
    "3.40.710.10": ("NylB', 6-aminohexanoate-dimer hydrolase (hybrid variant)",
                    "Paenarthrobacter ureafaciens", "CAA26616.1", "PDB 2ZLY", "451"),
    "1.10.760.10": ("RoxA, rubber oxygenase", "Steroidobacter cummioxidans 35Y (Xanthomonas sp. 35Y)",
                    "AAP41848.1", "PDB 4B2N", "107"),
    "2.60.40.420": ("LMC021, laccase-like multicopper oxidase", "Rhodococcus opacus R7",
                    "UniProt A0A076EYX8", "AlphaFold AF-A0A076EYX8-F1", "558"),
    "3.40.50.1110": ("E4, SGNH-family esterase", "Altericroceibacterium indicum",
                     "UniProt A0A4P1LYH5", "PDB 7C82", "463"),
    "1.10.420.10;1.10.520.10": ("MnP1, manganese peroxidase isozyme 1", "Phanerodontia chrysosporium",
                                "AAA33744.1", "PDB 1MNP", "177"),
    # the sheet gives BAA00001.1 (97.7%); the model is CYP3A4 (99.6%)
    "1.10.630.10": ("CYP3A4, cytochrome P450", "Homo sapiens", "UniProt P08684", "PDB 1W0E", "221"),
    "2.40.70.10": ("pepsin (pepsinogen)", "Homo sapiens", "AAA98529.1", "PDB 1PSO", None),
    "3.20.20.370": ("conserved hypothetical protein Pfl01_1497 family", "Pseudomonas fluorescens Pf0-1",
                    "ABA73449.1", "PDB 3CL6", None),
    "3.60.21.10": ("BsMPPE, metallophosphoesterase", "Bacillus subtilis", "UniProt A0A063XDV7",
                   "AlphaFold AF-A0A063XDV7-F1", "466"),
}


REPRESENTATIVE_FAA = os.path.join(HERE, "representatives.faa")


def source_sequence(k):
    """The full-length source protein of a fold's model."""
    acc = REPRESENTATIVE[k][2].replace("UniProt ", "")
    seq, take = [], False
    with open(REPRESENTATIVE_FAA) as fh:
        for line in fh:
            if line.startswith(">"):
                take = line[1:].split()[0] == acc
            elif take:
                seq.append(line.strip())
    return "".join(seq)


def representative(k, short=False):
    """'IsPETase, Piscinibacter sakaiensis (GAP38373.1; PDB 5XJH)'."""
    if k not in REPRESENTATIVE:
        return "none yet"
    enz, org, acc, ref, rec = REPRESENTATIVE[k]
    if short:
        return f"{enz.split(',')[0]} ({ref.split()[-1]})"
    where = (f"plasticome record {rec}" if rec else "reference protein, not one of the 609 records")
    return f"{enz}, {org} ({acc}; {ref}); {where}"


def structure_code(k):
    """The structure's identifier alone: '5XJH', or 'AF-A0A076EYX8-F1'."""
    if k not in REPRESENTATIVE:
        return ""
    return REPRESENTATIVE[k][3].replace("PDB ", "").replace("AlphaFold ", "")


def hsl(h, s, l):
    r, g, b = colorsys.hls_to_rgb(h / 360, l / 100, s / 100)
    return "#%02x%02x%02x" % (round(r * 255), round(g * 255), round(b * 255))


def fold_key(tip):
    """The key of a tree_tips.tsv row, with the unconfirmed component split off."""
    k = tip["fold_key"]
    return UNCONFIRMED if tip.get("homology") == "unconfirmed" else k


def _fold_tint(k):
    h, s, *_ = FOLDS[k]
    return hsl(h, s, 89)


def _fold_band(k):
    h, s, *_ = FOLDS[k]
    return hsl(h, min(s + 8, 70), 72)


# Folds whose label and strong shades are vivid rather than muted, so they
# stay apart from a neighbour of the same hue: Ntn amidase (red-orange, next
# to the muted red of the amidase signature; muted it read as brown) and
# subtilisin (blue-violet, moved off hue 175, which was too close to the aqua
# of the modelled plastic). Picked by CIELAB distance to every other fold and
# to the structure colours (aqua, pink, amber, slate, silver).
VIVID = {"3.60.70.12": (75, 52), "3.40.50.200": (70, 48)}


def _fold_label(k):
    h, s, *_ = FOLDS[k]
    if k in VIVID:
        return hsl(h, *VIVID[k])
    return hsl(h, 40 if s > 20 else 10, 50)


def _fold_strong(k):
    h, s, *_ = FOLDS[k]
    if k in VIVID:
        return hsl(h, min(VIVID[k][0] + 5, 85), 40)
    return hsl(h, 60 if s > 20 else 12, 38)


def name(k):
    return FOLDS[k][2]


# The names drawn in the ring tree (Figure 2a): the "PETadex Domain Name" column
# of the cath_domain_annotations sheet of plasticome.v2.xlsx (first letter
# capitalised). name() stays as it is, because the card file names and the
# heatmap's ring_order.tsv are built from it.
TITLE = {
    "3.40.50.1820": "Alpha/beta hydrolase",
    "3.40.50.1110": "SGNH hydrolase",
    "3.90.1300.10": "Amidase signature domain",
    "3.40.710.10": "DD-peptidase/beta-lactamase",
    "3.60.70.12": "DmpA aminopeptidase",   # shortened 2026-11-05; sheet: L-amino peptidase D-ALA esterase/amidase
    "2.40.10.10": "Trypsin-like serine protease",
    "3.40.50.200": "Peptidase S8/S53 domain",
    "1.10.760.10": "Cytochrome c-like",
    "domain:oxygenase MpaB family protein": "Oxygenase-like",
    "3.20.20.370": "Glycoside hydrolase/deacetylase",
    "1.10.420.10;1.10.520.10": "Peroxidase",
    "1.10.630.10": "Cytochrome P450",
    "2.40.70.10": "Acid protease",
    "3.60.21.10": "Metallo-dependent phosphatases",
    "2.60.40.420": "Cupredoxins",   # shortened 2026-11-05; sheet: Cupredoxins - blue copper proteins
    UNCONFIRMED: "DD-peptidase/beta-lactamase, unconfirmed",   # not in the sheet
}


def title(k):
    return TITLE.get(k) or name(k)


def mechanism(k):
    return FOLDS[k][3]


# --- structure annotations ---------------------------------------------------

# annotation row -> fold key. Matched on the catalytic CATH, and by name for
# the oxygenase, which has none.
def _annot_key(row):
    cath = (row["Reference CATH"] or "").strip()
    if not cath or cath.lower() in ("nan", "na", "none"):
        return "domain:oxygenase MpaB family protein"
    for k in FOLDS:
        if cath in k.split(";") and k != UNCONFIRMED:
            return k
    raise KeyError(f"no fold for annotation {row['Domain Name']} ({cath})")


TOKEN = re.compile(r"(?:(?P<chain>[A-Z]):)?(?P<aa>[A-Z])(?P<num>\d+)(?P<flag>\[[^\]]*\])?")


def parse_sites(cell):
    """'Catalytic Triad (S128,D174,H205), Oxyanion Hole (Y55, M129)' ->
    [(group, [(chain, aa, num, flag), ...]), ...]."""
    out = []
    for m in re.finditer(r"([^,(][^(]*?)\s*\(([^)]*(?:\[[^\]]*\][^)]*)*)\)", cell or ""):
        group = m.group(1).strip(" ,")
        res = []
        for t in TOKEN.finditer(m.group(2)):
            res.append((t.group("chain") or "A", t.group("aa"), int(t.group("num")),
                        (t.group("flag") or "").strip("[]")))
        if res:
            out.append((group, res))
    return out


def annotations():
    """fold key -> dict of the annotation row, parsed."""
    out = {}
    with open(ANNOT_CSV, newline="") as fh:
        for row in csv.DictReader(fh):
            k = _annot_key(row)
            out[k] = {
                "domain": row["Domain Name"].strip(),
                "pdb": os.path.join(PDB_DIR, row["Fold file"].strip()),
                "reference": row["Reference PDB/AF"].strip(),
                "accession": row["Reference GenBank accession"].strip(),
                "primary": parse_sites(row["Catalytic Residues (fold numbering)"]),
                "secondary": parse_sites(row["Secondary Residues (fold numbering)"]),
                "mechanism_text": row["Catalytic Mechanism"].strip(),
                "features": row["Special structural features"].strip(),
                "ligand": row["Predicted Ligand"].strip(),
                "components": row["Components"].strip(),
                "qc": row["QC flags"].strip(),
            }
    return out


def ligand_name(a):
    """'PET dimer (C22H22O10); ligand-ipTM ...' -> ('PET dimer', 'PET')."""
    m = re.search(r"substrate chain L \(([^)]+)\)", a["components"])
    plastic = m.group(1) if m else "NR"
    return a["ligand"].split("(")[0].strip().rstrip(","), plastic
