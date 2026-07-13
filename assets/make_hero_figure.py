"""Render the README hero figure from the Monte Carlo results.

Reads the two 500-run campaigns in MC_Results/ (weighted vs. classical B-dot) and
renders assets/detumbling_mc_results.{png,webp}. Run from the repository root:

    python3 assets/make_hero_figure.py

Requires numpy, scipy, and matplotlib. This is a documentation asset generator; it
is not part of the simulator and MATLAB does not need it.
"""
import os

import numpy as np
import scipy.io as sio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8983"
W = "#2a78d6"   # weighted B-dot  (categorical slot 1)
C = "#e34948"   # classical B-dot (categorical slot 6)

mu = 3.9860044181e14
a = 6378.137e3 + 350e3
T_ORB = 2 * np.pi / np.sqrt(mu / a**3)


def load(name):
    d = sio.loadmat(f"{REPO}/MC_Results/{name}.mat")["data"]
    d = d[:, ~np.isnan(d[0])]
    return {
        "t_det": d[0] / T_ORB,          # true detumbling time [orbits]
        "t_on": d[18:21],               # magnetorquer on-time per axis [s]
        "n": d.shape[1],
    }


wb = load("data_350_32_1_new3")   # weighted B-dot
cb = load("data_350_32_0_new3")   # classical B-dot

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 10.5,
    "axes.edgecolor": "#d8d7d2",
    "axes.labelcolor": INK2,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "text.color": INK,
})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 5.0), facecolor=SURFACE)
fig.subplots_adjust(left=0.055, right=0.985, top=0.735, bottom=0.13, wspace=0.22)

# ---------------------------------------------------------------- panel (a)
bins = np.linspace(9, 28, 30)
for data, color, label in ((cb, C, "Classical B-dot"), (wb, W, "Weighted B-dot")):
    ax1.hist(data["t_det"], bins=bins, color=color, alpha=0.55, edgecolor=SURFACE, linewidth=1.2, label=label)
for data, color in ((cb, C), (wb, W)):
    ax1.axvline(np.median(data["t_det"]), color=color, lw=2, ls=(0, (4, 2)))

ax1.set_facecolor(SURFACE)
ax1.grid(axis="y", color="#e8e7e2", lw=0.8)
ax1.set_axisbelow(True)
for s in ("top", "right"):
    ax1.spines[s].set_visible(False)
ax1.set_xlabel("Detumbling time  [orbits]")
ax1.set_ylabel(f"Monte Carlo runs  (N = {wb['n']})")
ax1.set_title("Both laws stop the tumble in ~15 orbits (~23 h)",
              loc="left", fontsize=11.5, color=INK, pad=10, fontweight="bold")

med_w, med_c = np.median(wb["t_det"]), np.median(cb["t_det"])
ax1.annotate("median\n14 orbits", xy=(med_w, ax1.get_ylim()[1] * 0.72), xytext=(9, 0),
             textcoords="offset points", color=INK2, fontsize=9.5)
ax1.legend(frameon=False, loc="upper right", labelcolor=INK2, fontsize=9.5)

# ---------------------------------------------------------------- panel (b)
axes_lbl = ["x", "y", "z"]
w_h = wb["t_on"].mean(axis=1) / 3600
c_h = cb["t_on"].mean(axis=1) / 3600
y = np.arange(3)[::-1]
h = 0.34

ax2.barh(y + h / 2 + 0.01, c_h, height=h, color=C, label="Classical B-dot")
ax2.barh(y - h / 2 - 0.01, w_h, height=h, color=W, label="Weighted B-dot")

for i, (yy, wv, cv) in enumerate(zip(y, w_h, c_h)):
    ax2.text(cv + 0.18, yy + h / 2 + 0.01, f"{cv:.1f} h", va="center", color=INK2, fontsize=9.5)
    red = 100 * (cv - wv) / cv
    ax2.text(wv + 0.18, yy - h / 2 - 0.01, f"{wv:.1f} h   −{red:.1f}%", va="center",
             color=W, fontsize=9.5, fontweight="bold")

ax2.set_facecolor(SURFACE)
ax2.set_yticks(y, [f"MTQ {c}" for c in axes_lbl])
ax2.set_xlim(0, 19.5)
ax2.grid(axis="x", color="#e8e7e2", lw=0.8)
ax2.set_axisbelow(True)
for s in ("top", "right", "left"):
    ax2.spines[s].set_visible(False)
ax2.tick_params(axis="y", length=0)
ax2.set_xlabel("Mean magnetorquer on-time until detumbled  [hours]")
tot = 100 * (c_h.sum() - w_h.sum()) / c_h.sum()
ax2.set_title(f"…but the weighted law uses {tot:.1f}% less actuation",
              loc="left", fontsize=11.5, color=INK, pad=10, fontweight="bold")
ax2.legend(frameon=False, loc="upper right", labelcolor=INK2, fontsize=9.5,
           bbox_to_anchor=(1.0, 1.02))

fig.text(0.055, 0.915, "Magnetic detumbling of a picosatellite from 180 °/s on all three axes",
         fontsize=14, fontweight="bold", color=INK)
fig.text(0.055, 0.855,
         f"{wb['n']} Monte Carlo runs per controller · 350 km SSO · 4 Hz control loop · "
         "dispersed mass, inertia, dipole, and sensor errors",
         fontsize=10, color=MUTED)

out = f"{REPO}/assets/detumbling_mc_results"
fig.savefig(out + ".png", dpi=170, facecolor=SURFACE)
fig.savefig(out + ".webp", dpi=170, facecolor=SURFACE)
print("wrote", out + ".{png,webp}")
print("weighted total %.2f h, classical total %.2f h, reduction %.2f%%" % (w_h.sum(), c_h.sum(), tot))
print("median orbits: weighted %.2f, classical %.2f" % (med_w, med_c))
print("mean orbits: weighted %.2f, classical %.2f" % (wb["t_det"].mean(), cb["t_det"].mean()))
