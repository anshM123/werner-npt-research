"""Figure: proven n-copy windows for qutrit Werner states (d = 3).

Left: windows alpha_n from the budget theorem (four kinds of input).
Right: the excess n(alpha_n - 1/3), with the budget, the certified correlated-factor windows (Section 7), their
all-n feed-back from n = 24, and the proved ceiling for decompositions with cores on at most m sites.
"""
from fractions import Fraction as Fr

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

d = 3


def alpha_n(n, m, a0):
    a0 = Fr(a0)
    return float((n * (1 + a0) + m * (d * a0 - 1)) / (n * d * (1 + a0) - m * (d * a0 - 1)))


BLUE, GREEN, ORANGE, GRAY, PURPLE = '#0072B2', '#009E73', '#D55E00', '#666666', '#7B3FA0'
ns = list(range(2, 31))
fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 3.3), dpi=200)

series = [
    (r'three copies, endpoint ($m=3$)', 3, Fr(1, 2), BLUE, 'o'),
    (r'three copies at $\alpha_0=37/80$ (superseded)', 3, Fr(37, 80), GREEN, 's'),
    (r'two copies, endpoint ($m=2$)', 2, Fr(1, 2), ORANGE, '^'),
    (r'five copies at $\alpha_0=47/100$ ($m=5$)', 5, Fr(47, 100), PURPLE, 'D'),
]
for label, m, a0, col, mk in series:
    xs = [n for n in ns if n >= m]
    ax.plot(xs, [alpha_n(n, m, a0) for n in xs], color=col, lw=1.8, marker=mk, ms=3, label=label)
ax.axhline(1 / 3, color=GRAY, lw=1.0, ls='--')
ax.text(2.0, 1 / 3 + 0.003, r'PPT threshold $1/d$', ha='left', va='bottom', fontsize=7, color='#333333')
ax.axhline(1 / 2, color='#999999', lw=0.9, ls=':')
ax.text(5.5, 0.5 - 0.003, r'endpoint $\alpha=1/2$', ha='left', va='top', fontsize=7, color='#333333')
ax.set_xlabel(r'number of copies $n$', fontsize=8)
ax.set_ylabel(r'window $\alpha_n$', fontsize=8)
ax.set_title('Budget windows', fontsize=8.5)
ax.set_xlim(1.5, 30.5); ax.set_ylim(0.32, 0.51)
ax.legend(fontsize=6.0, frameon=False, loc='center right', bbox_to_anchor=(1.0, 0.66))

# right panel: excess n(alpha - 1/3)
cert3 = {4: 0.45756, 5: 0.43269, 6: 0.41616, 8: 0.39543, 10: 0.38299, 12: 0.37451, 16: 0.36436, 24: 0.35403}
cert5 = {6: 0.44697, 7: 0.43062, 8: 0.41845, 10: 0.40100, 12: 0.39000, 16: 0.37580, 24: 0.36160}
cert2 = {4: 0.41603, 5: 0.39954, 6: 0.38850, 8: 0.37469, 10: 0.36622, 12: 0.36070, 24: 0.34712}


def ex(n, a):
    return n * (a - 1 / 3)


for m, a0, col, mk, cert, fb in ((3, Fr(1, 2), BLUE, 'o', cert3, (135403, 149016, 406209)),
                                 (2, Fr(1, 2), ORANGE, '^', cert2, (16839, 12408, 50517)),
                                 (5, Fr(47, 100), PURPLE, 'D', cert5, (851, 1272, 2553))):
    xs = [n for n in ns if n >= max(m, 3)]
    bx.plot(xs, [ex(n, alpha_n(n, m, a0)) for n in xs], color=col, lw=1.4, ls='-', alpha=0.55)
    bx.plot(list(cert), [ex(n, a) for n, a in cert.items()], color=col, lw=0, marker=mk, ms=4.2,
            markeredgecolor='white', markeredgewidth=0.6)
    A, B, C = fb
    xs2 = list(range(24, 31))
    bx.plot(xs2, [ex(n, (A * n + B) / (C * n - B)) for n in xs2], color=col, lw=1.6, ls=(0, (1.2, 1.2)))
    cap = [8 * m * n / (45 * n + 3 * m) for n in xs]
    bx.plot(xs, cap, color=col, lw=1.0, ls='--', alpha=0.9)
bx.text(30, 8 / 15 + 0.004, r'ceiling, $m=3$', ha='right', va='bottom', fontsize=6.5, color=BLUE)
bx.text(30, 0.4444 - 0.004, r'budget, $m=3$', ha='right', va='top', fontsize=6.5, color=BLUE)
bx.text(30, 16 / 45 + 0.004, r'ceiling, $m=2$', ha='right', va='bottom', fontsize=6.5, color=ORANGE)
bx.text(30, 0.2963 - 0.004, r'budget, $m=2$', ha='right', va='top', fontsize=6.5, color=ORANGE)
bx.text(30, 0.62 - 0.006, r'budget, $m=5$', ha='right', va='top', fontsize=6.5, color=PURPLE)
bx.text(30, 8 / 9 - 0.006, r'ceiling, $m=5$', ha='right', va='top', fontsize=6.5, color=PURPLE)
bx.text(3.0, 0.80, 'dots: exact certificates\nwith correlated factors\ndotted: all-$n$ feed-back from $n=24$',
        fontsize=6.3, color='#333333', va='center')
bx.set_xlabel(r'number of copies $n$', fontsize=8)
bx.set_ylabel(r'excess $n(\alpha_n-1/3)$', fontsize=8)
bx.set_title('Excess over the PPT threshold', fontsize=8.5)
bx.set_xlim(1.5, 30.5); bx.set_ylim(0.26, 0.92)
for a_ in (ax, bx):
    a_.grid(True, color='#e9e9e9', lw=0.6)
    a_.tick_params(labelsize=7)
    for s in ('top', 'right'):
        a_.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig('windows_d3.pdf')
fig.savefig('windows_d3.png')
print('saved')
