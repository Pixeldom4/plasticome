#!/usr/bin/env python3
"""Generate a sequence logo from a HMMER3 profile HMM.

Reconstructs match-state emission probabilities from the .hmm file
(scores are -ln(p)), then renders a logo with logomaker using the
Skylign/hmmlogo default convention:

  * column total height = relative entropy (KL divergence) of the
    emission distribution against HMMER's standard amino-acid background,
    in bits.
  * each residue's height within the column is proportional to its
    emission probability.
"""
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker

# HMMER default amino-acid background frequencies (Swiss-Prot averages,
# p7_AMINO_FREQS in the HMMER source).
BG = {
    'A': 0.0787945, 'C': 0.0151600, 'D': 0.0535222, 'E': 0.0668298,
    'F': 0.0397062, 'G': 0.0695071, 'H': 0.0229198, 'I': 0.0590092,
    'K': 0.0594422, 'L': 0.0963728, 'M': 0.0237718, 'N': 0.0414386,
    'P': 0.0482904, 'Q': 0.0395639, 'R': 0.0540978, 'S': 0.0683364,
    'T': 0.0540687, 'V': 0.0673194, 'W': 0.0114135, 'Y': 0.0304133,
}


def parse_hmm(path):
    alphabet = None
    name = None
    rows = []          # per match state: dict letter -> probability
    consensus = []     # consensus residue per position
    with open(path) as fh:
        lines = fh.readlines()

    for i, line in enumerate(lines):
        toks = line.split()
        if not toks:
            continue
        if toks[0] == 'NAME':
            name = toks[1]
        elif toks[0] == 'HMM' and alphabet is None:
            alphabet = toks[1:]           # 20 amino-acid letters, in file order
        elif alphabet is not None and toks[0].isdigit():
            # match-state emission line: <pos> <20 scores> <MAP> <cons> ...
            scores = toks[1:1 + len(alphabet)]
            probs = {}
            for a, s in zip(alphabet, scores):
                probs[a] = 0.0 if s == '*' else math.exp(-float(s))
            total = sum(probs.values())
            if total > 0:
                probs = {a: p / total for a, p in probs.items()}
            rows.append(probs)
            # consensus residue (CONS yes) sits just after the MAP number
            cons_field = toks[1 + len(alphabet) + 1] if len(toks) > 1 + len(alphabet) + 1 else '-'
            consensus.append(cons_field.upper())

    return name, alphabet, rows, consensus


def relative_entropy_heights(rows, alphabet):
    """Return a DataFrame of letter heights (bits) per position."""
    data = []
    for probs in rows:
        # KL divergence vs background -> column height in bits
        col_bits = 0.0
        for a in alphabet:
            p = probs[a]
            if p > 0:
                col_bits += p * math.log2(p / BG[a])
        col_bits = max(col_bits, 0.0)
        # letters proportional to emission probability
        data.append({a: probs[a] * col_bits for a in alphabet})
    df = pd.DataFrame(data, columns=alphabet)
    df.index = np.arange(1, len(df) + 1)   # 1-based positions
    return df


def main():
    hmm_path = sys.argv[1] if len(sys.argv) > 1 else \
        '/Users/Pixel/Documents/projects/plasticome/hmm-seq-logo/pazy_hmm_c1_260507.hmm'
    out_base = sys.argv[2] if len(sys.argv) > 2 else hmm_path.rsplit('.', 1)[0]

    name, alphabet, rows, consensus = parse_hmm(hmm_path)
    print(f'model: {name}   length: {len(rows)}   alphabet: {"".join(alphabet)}')

    heights = relative_entropy_heights(rows, alphabet)
    max_bits = heights.sum(axis=1).max()
    print(f'max column info content: {max_bits:.2f} bits')

    n = len(heights)
    # single unwrapped strip; width scales with model length
    fig, ax = plt.subplots(figsize=(n * 0.18, 2.6))

    logomaker.Logo(heights, ax=ax, color_scheme='chemistry',
                   show_spines=False, vpad=0.02)
    ax.set_ylim(0, max_bits * 1.05)
    ax.set_xlim(0.5, n + 0.5)
    ax.set_ylabel('bits')
    ax.set_xticks(list(heights.index))
    ax.set_xticklabels(list(heights.index), fontsize=4, rotation=90)
    ax.tick_params(axis='y', labelsize=7)

    ax.set_title(f'Sequence logo — {name}  ({n} match states, relative-entropy heights)',
                 fontsize=11)
    fig.tight_layout()

    png = out_base + '_logo.png'
    svg = out_base + '_logo.svg'
    fig.savefig(png, dpi=200, bbox_inches='tight')
    fig.savefig(svg, bbox_inches='tight')
    print(f'wrote {png}')
    print(f'wrote {svg}')


if __name__ == '__main__':
    main()
