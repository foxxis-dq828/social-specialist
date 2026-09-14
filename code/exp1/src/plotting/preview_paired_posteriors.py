"""Publication-layout example using the MATLAB plotting script's defaults."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'figures' / 'paired_posterior_example'
OUT.mkdir(parents=True, exist_ok=True)

# Positive-scale defaults in plot_custom_model_posteriors.m; ce does not affect w.
n, e = 1.0, 1.0
m = np.array([.8, .2380952, .0625, .4571429, .7999999, .8, .7999999, .8])
d, de = 2 + 2*n, 2 + 2*n*e
vt = np.array([.5, 1/d, 1/d, 1/d, (1+2*n)/d, (1+2*n)/d, .5, (1+2*n)/d])
dw = np.array([.5, 1/d, 1/de, 1/d, (1+2*n*e)/de, (1+2*n)/d, .5, (1+2*n)/d])
blue = np.array([1.870504, 1.108434, 1.108434, 1.108434])
yellow = np.array([3.085106, 1.870504, 1.108434, 1.870504])
kw = np.array([
    (1+yellow[0]*n)/(2+(blue[0]+yellow[0])*n), 1/(2+2*blue[0]*n),
    1/(2+(blue[1]+yellow[1])*n), 1/(2+2*blue[1]*n),
    (1+(blue[2]+yellow[2])*n)/(2+(blue[2]+yellow[2])*n),
    (1+2*blue[2]*n)/(2+2*blue[2]*n),
    (1+blue[3]*n)/(2+(blue[3]+yellow[3])*n),
    (1+2*blue[3]*n)/(2+2*blue[3]*n)])
weights = np.column_stack([m, vt, dw, kw])
assert weights.shape == (8, 4) and np.all((weights >= 0) & (weights <= 1))
np.testing.assert_allclose(vt, dw)
# Plot the complementary probability for Trials 1 and 2 only.
weights[:2] = 1 - weights[:2]
np.savetxt(OUT / 'posterior_weights.csv', weights, delimiter=',',
           header='M,VT,DW,KW', comments='')

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
    'axes.edgecolor': '#555555', 'axes.linewidth': .65,
    'text.color': '#242424', 'axes.labelcolor': '#242424',
    'xtick.color': '#444444', 'ytick.color': '#444444',
    'pdf.fonttype': 42, 'svg.fonttype': 'none'})
colors = ['#C39A38', '#438675', '#A27FA9', '#587FA5']
fig, axes = plt.subplots(4, 2, figsize=(9.0, 10.8),
                         gridspec_kw={'width_ratios': [1.5, 1]})
fig.subplots_adjust(left=.105, right=.975, bottom=.105, top=.81,
                    wspace=.32, hspace=.70)
fig.text(.105, .963, 'Model predictions across paired trials', fontsize=17,
         weight='bold')
fig.text(.105, .933, 'Exp. 1  |  Posterior weights at illustrative parameter settings',
         fontsize=10, color='#666666')
fig.legend(handles=[Patch(facecolor='#777777', edgecolor='#555555',
                          label='First trial in pair'),
                    Patch(facecolor='white', edgecolor='#555555',
                          label='Second trial in pair')],
           loc='upper left', bbox_to_anchor=(.098, .916), ncol=2,
           frameon=False, handlelength=1.35, columnspacing=2.2)
fig.text(.105, .864, 'a   Posterior probability', fontsize=12, weight='bold')
fig.text(.668, .864, 'b   Within-pair contrast', fontsize=12, weight='bold')
x = np.arange(4)
for row, (a, b) in enumerate(axes):
    first, second = weights[2*row:2*row+2]
    delta = first - second
    a.bar(x-.17, first, width=.30, color=colors, edgecolor=colors, linewidth=1.2)
    a.bar(x+.17, second, width=.30, color='white', edgecolor=colors, linewidth=1.2)
    a.set_ylim(0, 1)
    a.set_yticks([0, .5, 1], ['0', '0.5', '1'])
    a.set_ylabel(r'$1-P(\mathrm{Type\ II}\mid s_1,s_2)$' if row == 0
                 else r'$P(\mathrm{Type\ II}\mid s_1,s_2)$')
    a.set_title(f'Pair {row+1}   |   Trial {2*row+1} vs {2*row+2}',
                loc='left', fontsize=10, pad=9, weight='bold')
    b.axhline(0, color='#555555', lw=.9, zorder=1)
    b.vlines(x, 0, delta, color=colors, lw=2.2, zorder=2)
    b.scatter(x, delta, color=colors, s=42, zorder=3)
    for i, value in enumerate(delta):
        label = '0.00' if abs(value) < .005 else f'{value:+.2f}'
        b.annotate(label, (i, value), xytext=(0, 9 if value >= -1e-6 else -14),
                   textcoords='offset points', ha='center', fontsize=9)
    b.set_ylim(-.8, .8)
    b.set_yticks([-.5, 0, .5], ['-0.5', '0', '+0.5'])
    b.set_ylabel(r'$\Delta(1-P)$' if row == 0 else r'$\Delta P$')
    b.set_title(f'Trial {2*row+1} $-$ Trial {2*row+2}',
                loc='left', fontsize=10, pad=9, color='#666666')
    for ax in (a, b):
        ax.set_xlim(-.6, 3.6)
        ax.set_xticks(x, ['M', 'VT', 'DW', 'KW'])
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(axis='x', length=0, pad=7)
        ax.tick_params(axis='y', length=3)
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color='#EAEAEA', linewidth=.55)
fig.text(.105, .052, 'M, mentalizing; VT, vote tallying; DW and KW follow the model-file definitions.',
         fontsize=8, color='#666666')
fig.text(.105, .033, 'Defaults: n = 1, e = 1. Points and bars are model predictions, not fitted data or uncertainty intervals.',
         fontsize=8, color='#666666')
for suffix in ('png', 'pdf', 'svg'):
    fig.savefig(OUT / f'paired_posteriors.{suffix}', dpi=220, facecolor='white')
print(OUT / 'paired_posteriors.png')
