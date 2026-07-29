#!/usr/bin/env python3
"""Domain schematic of IsPETase (Ideonella sakaiensis PETase, 290 aa).

Draws the protein as a rectangle with the N-terminal signal peptide (1-27)
highlighted and the catalytic Ser-His-Asp triad marked as lollipops.
"""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle

SEQ = ("MNFPRASRLMQAAVLGGLMAVSAAATAQTNPYARGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTV"
       "YYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVAS"
       "LNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACEN"
       "DSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACE"
       "NPNSTRVSDFRTANCS")
N = len(SEQ)                                    # 290
SIGNAL = (1, 27)                                # signal peptide
DOMAIN = (65, 261)                              # HMM domain hit
TRIAD = [(160, 'Ser160'), (206, 'Asp206'), (237, 'His237')]

# palette
C_MATURE = '#3d6cb5'
C_SIGNAL = '#e0a53b'
C_TRIAD  = '#c0392b'
C_BACK   = '#c9d1dc'
INK      = '#1a1a1a'

fig, ax = plt.subplots(figsize=(12, 3.0))

bar_y, bar_h = 0.0, 1.0

# full-length backbone (thin), then overlaid feature boxes
back_h = bar_h * 0.42
back_y = bar_y + (bar_h - back_h) / 2
ax.add_patch(FancyBboxPatch((1, back_y), N - 1, back_h,
             boxstyle='round,pad=0,rounding_size=4', mutation_aspect=0.04,
             linewidth=1.0, edgecolor=INK, facecolor=C_BACK, zorder=2))

# domain hit (65-261)
ax.add_patch(Rectangle((DOMAIN[0], bar_y), DOMAIN[1] - DOMAIN[0], bar_h,
             linewidth=1.2, edgecolor=INK, facecolor=C_MATURE, zorder=3))
# signal peptide (1-27)
ax.add_patch(Rectangle((SIGNAL[0], bar_y), SIGNAL[1] - SIGNAL[0], bar_h,
             linewidth=1.2, edgecolor=INK, facecolor=C_SIGNAL, zorder=3))

# label inside the domain box
ax.text((DOMAIN[0] + DOMAIN[1]) / 2, bar_y + bar_h / 2,
        f'α/β-hydrolase domain ({DOMAIN[0]}–{DOMAIN[1]})',
        ha='center', va='center', color='white', fontsize=11,
        fontweight='bold', zorder=4)

# signal peptide bracket + label below
ax.annotate('', xy=(1, -0.35), xytext=(SIGNAL[1], -0.35),
            arrowprops=dict(arrowstyle='-', color=C_SIGNAL, lw=2))
ax.text((1 + SIGNAL[1]) / 2, -0.62, f'Signal peptide\n({SIGNAL[0]}–{SIGNAL[1]})',
        ha='center', va='top', fontsize=9, color=INK)

# catalytic triad lollipops above the bar
stem_top = bar_h + 0.55
for pos, lab in TRIAD:
    ax.plot([pos, pos], [bar_h, stem_top], color=C_TRIAD, lw=1.6, zorder=5)
    ax.scatter([pos], [stem_top], s=90, color=C_TRIAD, zorder=6,
               edgecolor='white', linewidth=0.8)
    ax.text(pos, stem_top + 0.12, lab, ha='center', va='bottom',
            fontsize=9.5, color=C_TRIAD, fontweight='bold')

ax.text((TRIAD[0][0] + TRIAD[-1][0]) / 2, stem_top + 0.55,
        'Catalytic triad', ha='center', va='bottom', fontsize=10,
        color=C_TRIAD, style='italic')

# residue ruler
ticks = [1, 27, 65, 100, 150, 160, 206, 237, 261, 290]
for t in ticks:
    ax.plot([t, t], [bar_y - 0.06, bar_y], color=INK, lw=0.8, zorder=2)
    ax.text(t, bar_y - 0.12, str(t), ha='center', va='top', fontsize=6.5, color=INK)

ax.set_title('IsPETase  ($Ideonella\\ sakaiensis$ PET hydrolase, 290 aa)',
             fontsize=13, pad=14)
ax.set_xlim(-8, N + 8)
ax.set_ylim(-1.1, 2.4)
ax.axis('off')
fig.tight_layout()

for ext in ('png', 'svg'):
    out = f'/Users/Pixel/Documents/projects/plasticome/hmm-seq-logo/ispetase_domain.{ext}'
    fig.savefig(out, dpi=200, bbox_inches='tight')
    print('wrote', out)
