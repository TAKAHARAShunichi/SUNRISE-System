#!/usr/bin/env python3
"""Regenerate the figures.

    python scripts/figures.py                 # all of them
    python scripts/figures.py timecourse      # just one
    python scripts/figures.py heatmap --metric nadir --nx 30 --ny 20

The heatmaps integrate the full course at every grid point, so a 60 x 40 grid
takes a few minutes on two cores. Pass smaller --nx/--ny while iterating.

Outputs land in figures/ as PNG at 200 dpi.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm  # noqa: E402

from sunrise import metrics, parameters as prm  # noqa: E402
from sunrise.model import simulate  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "figures"

INK, INK_SOFT, STATUS = "#1A1A18", "#5B5B55", "#A33122"
GREY, WARM, COOL = "#8A8A82", "#D1794B", "#145E77"
# Diverging ramp: cool = passing, neutral midpoint, warm = failing.
RAMP = ["#08303F", "#145E77", "#4C9DB8", "#A8CBD6", "#DCDCD8",
        "#E6B48C", "#D1794B", "#A8482A", "#6B2416"]


def _style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color("#C9C9C3")
    ax.tick_params(colors=INK_SOFT, labelsize=9.5)


# ------------------------------------------------------------------ figure 1
def timecourse(p: dict, out: Path) -> None:
    """Extracellular and intracellular trajectories for three regimens."""
    arms = [("Untreated", {"a_peak": 0.0, "dose": 0.0}, GREY, (0, (1, 2))),
            ("Kanamycin alone", {"dose": 0.0}, WARM, (0, (6, 2))),
            ("Kanamycin + SUNRISE", {}, COOL, "-")]

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.2), sharey=True)
    runs = {}
    for lab, over, col, dash in arms:
        t, ex, inn, _ = simulate(p, **over)
        runs[lab] = (t, ex, inn)
        axes[0].plot(t, np.maximum(ex, 1e-4), color=col, lw=2.2, label=lab)
        axes[1].plot(t, np.maximum(inn, 1e-4), color=col, ls=dash,
                     lw=3.4 if lab == "Untreated" else 2.2, solid_capstyle="round")

    for ax, title in zip(axes, ("Extracellular  $C^\\mathrm{ex}$",
                                "Intracellular  $C^\\mathrm{in}$")):
        ax.set_yscale("log")
        ax.set_xlim(0, p["t_end"])
        ax.set_ylim(1e-3, 3e9)
        ax.set_xticks(range(0, int(p["t_end"]) + 1, 12))
        ax.set_xlabel("Time  (h)", fontsize=10.5, color=INK)
        ax.set_title(title, fontsize=11.5, color=INK, loc="left", pad=8)
        ax.axhline(p["cure_density"], color=INK, lw=1.2, ls=(0, (4, 3)))
        ax.grid(axis="y", color="#EDEDE9", lw=0.8)
        ax.set_axisbelow(True)
        _style(ax)

    axes[0].set_ylabel("Bacterial density  (CFU/mL)", fontsize=10.5, color=INK)
    axes[0].text(p["t_end"] - 2, p["cure_density"] * 1.4, "$10^2$ CFU/mL",
                 color=INK_SOFT, fontsize=9, ha="right", va="bottom")
    axes[1].axhline(p["k_in"], color=STATUS, lw=1.6, ls=(0, (2, 2)))
    axes[1].text(p["t_end"] - 2, p["k_in"] * 1.7,
                 f"$K^\\mathrm{{in}}$ = {p['k_in']:g}", color=STATUS,
                 fontsize=9, ha="right", va="bottom")
    axes[1].annotate("all three sit at capacity:\ntreatment never empties\n"
                     "the reservoir",
                     xy=(0.60 * p["t_end"], p["k_in"] * 0.78),
                     xytext=(0.56 * p["t_end"], 2.0e-2), fontsize=9.5,
                     color=INK_SOFT, ha="center", linespacing=1.5,
                     arrowprops=dict(arrowstyle="-", color=INK_SOFT, lw=1.0))

    for lab, col, dx in (("Kanamycin alone", WARM, -2), ("Kanamycin + SUNRISE", COOL, 13)):
        t, ex, inn = runs[lab]
        tot = ex + inn
        i = int(np.argmin(tot))
        axes[0].annotate(f"nadir {np.log10(tot[0] / tot.min()):.1f} log",
                         xy=(t[i], ex[i]), xytext=(t[i] + dx, ex[i] * 0.06),
                         fontsize=9, color=col, ha="center",
                         arrowprops=dict(arrowstyle="-", color=col, lw=1.0))

    h, l = axes[0].get_legend_handles_labels()
    leg = fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.5, 0.945),
                     ncol=3, frameon=False, fontsize=10)
    for txt in leg.get_texts():
        txt.set_color(INK)
    fig.suptitle("Kanamycin alone reaches the nadir but does not hold it",
                 fontsize=13.5, color=INK, x=0.008, ha="left", y=0.992)
    fig.text(0.008, 0.012, _footer(p), fontsize=8.2, color=INK_SOFT)
    fig.tight_layout(rect=(0, 0.035, 1, 0.885))
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("wrote", out)


def _footer(p: dict) -> str:
    return (f"{p['t_end']:.0f} h course · $\\tau_A$ = {p['tau_a']:.0f} h, "
            f"$\\tau_P$ = {p['tau_p']:.0f} h · $\\varphi P$ = {p['phi_p']:.1f} · "
            f"$K^\\mathrm{{in}}$ = {p['k_in']:g} CFU/mL · baseline parameters "
            "otherwise.")


# ------------------------------------------------------------------ figure 2
_CTX: dict = {}


def _init(p, metric, A, D):
    _CTX.update(p=p, metric=metric, A=A, D=D)


def _cell(ij):
    i, j = ij
    r = metrics.evaluate(_CTX["p"], a_peak=_CTX["A"][i], dose=_CTX["D"][j])
    return i, j, (r["depth_log"] if _CTX["metric"] == "nadir"
                  else r["final_total"] if _CTX["metric"] == "total"
                  else r["final_ex"])


def heatmap(p: dict, out: Path, metric: str = "nadir", nx: int = 60,
            ny: int = 40, procs: int = 2) -> None:
    """Treatment outcome over antibiotic concentration and particle dose."""
    A = np.linspace(50.0, 700.0, nx)
    D = np.logspace(10, 13, ny)
    with mp.Pool(procs, initializer=_init, initargs=(p, metric, A, D)) as pool:
        cells = pool.map(_cell, [(i, j) for i in range(nx) for j in range(ny)],
                         chunksize=24)
    Z = np.empty((ny, nx))
    for i, j, v in cells:
        Z[j, i] = v

    nadir = metric == "nadir"
    if nadir:
        L, center = np.minimum(Z, 8.0), p["nadir_target"]
        cmap = LinearSegmentedColormap.from_list("c", RAMP[::-1])
        label = ("Log-reduction at nadir  "
                 "$\\log_{10}(C^\\mathrm{tot}_0 / \\min_t C^\\mathrm{tot})$")
        passing = L >= center
    else:
        L = np.minimum(np.log10(np.maximum(Z, 1e-3)), 6.0)
        center = np.log10(p["cure_density"])
        cmap = LinearSegmentedColormap.from_list("c", RAMP)
        label = ("Final-cycle geometric mean  "
                 + ("$\\bar{C}^\\mathrm{ex}+\\bar{C}^\\mathrm{in}$"
                    if metric == "total" else "$\\bar{C}^\\mathrm{ex}$")
                 + "  (CFU/mL)")
        passing = L <= center

    norm = TwoSlopeNorm(vmin=min(L.min(), center - 0.1), vcenter=center,
                        vmax=max(L.max(), center + 0.1))

    fig, ax = plt.subplots(figsize=(10.2, 6.4))
    mesh = ax.pcolormesh(A, D, L, cmap=cmap, norm=norm, shading="gouraud",
                         rasterized=True)
    if passing.any() and not passing.all():
        ax.contour(A, D, L, levels=[center], colors=INK, linewidths=2.2)
        ax.text(0.33, 0.30,
                ("Target\n%g-log reduction" % center) if nadir
                else "Cure threshold\n$10^{%g}$ CFU/mL" % center,
                transform=ax.transAxes, fontsize=10, color=INK,
                fontweight="bold", linespacing=1.4, ha="left", va="top")
    elif not passing.any():
        ax.text(0.70, 0.52, "No point in this domain passes",
                transform=ax.transAxes, fontsize=13, color=INK,
                fontweight="bold", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.6", fc="#FFFFFF", ec=INK,
                          lw=1.4, alpha=0.92))

    ax.axvline(p["a_ceiling"], color=STATUS, lw=2.4, zorder=4)
    ax.axvspan(p["a_ceiling"], A[-1], facecolor="none", hatch="///",
               edgecolor=STATUS, alpha=0.16, lw=0, zorder=3)
    ax.annotate("Nephrotoxicity and ototoxicity\nrise with dose: kanamycin "
                f"beyond\n{p['a_ceiling']:.0f} µg/mL exceeds the clinical ceiling",
                xy=(p["a_ceiling"] + 16, 6.2e12), color=STATUS, fontsize=10,
                ha="left", va="top", linespacing=1.5, zorder=6)

    base = metrics.evaluate(p)
    shown = base["depth_log"] if nadir else base["final_total"]
    ax.plot([p["a_peak"]], [p["dose"]], marker="o", ms=8, mfc="none", mec=INK,
            mew=2.0, zorder=5)
    ax.annotate(f"Baseline  {p['a_peak']:.0f} µg/mL, {p['dose']:.1e} particles\n"
                + (f"{shown:.2f} log" if nadir else f"{shown:.0f} CFU/mL"),
                xy=(p["a_peak"], p["dose"]), xytext=(86, 1.6e12), fontsize=9.5,
                color=INK, linespacing=1.45,
                arrowprops=dict(arrowstyle="-", color=INK_SOFT, lw=1.2))

    ax.set_yscale("log")
    ax.set_xlim(A[0], A[-1])
    ax.set_ylim(D[0], D[-1])
    ax.set_xlabel("Kanamycin peak urinary concentration  $A_\\mathrm{peak}$  "
                  f"(µg/mL, q{p['tau_a']:.0f}h)", fontsize=11, color=INK)
    ax.set_ylabel(f"SUNRISE dose  $D$  (particles/dose, q{p['tau_p']:.0f}h)",
                  fontsize=11, color=INK)
    ax.set_title("Treatment outcome over kanamycin and SUNRISE dose",
                 fontsize=13.5, color=INK, pad=12, loc="left")
    _style(ax)

    cb = fig.colorbar(mesh, ax=ax, pad=0.015,
                      extend="neither" if nadir else "max")
    cb.set_label(label, fontsize=10.5, color=INK)
    cb.ax.tick_params(colors=INK_SOFT, labelsize=9.5)
    cb.outline.set_visible(False)

    fig.text(0.008, 0.030, _footer(p), fontsize=8.2, color=INK_SOFT)
    fig.text(0.008, 0.008, "The clinical ceiling is derived from systemic "
             "dosing; for intravesical instillation it needs re-deriving.",
             fontsize=8.2, color=INK_SOFT)
    fig.tight_layout(rect=(0, 0.048, 1, 1))
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("wrote", out)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("which", nargs="?", default="all",
                    choices=["all", "timecourse", "heatmap"])
    ap.add_argument("--metric", default="nadir",
                    choices=["nadir", "total", "ex"])
    ap.add_argument("--nx", type=int, default=60)
    ap.add_argument("--ny", type=int, default=40)
    ap.add_argument("--procs", type=int, default=2)
    a = ap.parse_args()

    OUT.mkdir(exist_ok=True)
    # Use the calibrated degradation rate, as validate.py establishes it.
    from validate import calibrate  # noqa: E402
    p = prm.with_overrides(k_tpd0=calibrate(verbose=False))

    if a.which in ("all", "timecourse"):
        timecourse(p, OUT / "timecourse.png")
    if a.which in ("all", "heatmap"):
        mets = ["nadir", "total"] if a.which == "all" else [a.metric]
        for m in mets:
            heatmap(p, OUT / f"heatmap_{m}.png", metric=m, nx=a.nx, ny=a.ny,
                    procs=a.procs)


if __name__ == "__main__":
    main()
