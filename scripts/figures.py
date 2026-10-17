"""Figures for the Ver.23 recalculation."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from sunrise.design import N_OUT
from sunrise.metrics import evaluate
from sunrise.model import solve_course
from sunrise.params import DOSE_MAX_G, STANDARD_PATIENT

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "figures"
RESULTS = ROOT / "results"
FIG.mkdir(exist_ok=True)

# colour-vision-safe three-arm palette
C_UNTREATED = "#7a7a7a"
C_ABX = "#d95f02"
C_SUNRISE = "#1b9e77"
C_COMBO = "#7570b3"

REF_TAU_P, REF_MOI, REF_PHASE = 6.0, 1.409, 3.0


def ref(**kw):
    p = (
        STANDARD_PATIENT.with_(tau_P=REF_TAU_P, dt_phase=REF_PHASE)
        .with_dose_g(DOSE_MAX_G)
        .with_MOI(REF_MOI)
    )
    return p.with_(**kw) if kw else p


def fig_timecourse() -> None:
    arms = [
        ("Untreated", ref(sunrise_on=False, antibiotic_on=False), C_UNTREATED),
        ("Kanamycin only", ref(sunrise_on=False), C_ABX),
        ("SUNRISE only", ref(antibiotic_on=False), C_SUNRISE),
        ("Combination", ref(), C_COMBO),
    ]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharex=True)
    for label, p, col in arms:
        s = solve_course(p, n_out=N_OUT)
        ax[0].semilogy(s.t, np.maximum(s.C_ex, 1e-2), color=col, label=label, lw=1.8)
        ax[1].semilogy(s.t, np.maximum(s.C_in, 1e-2), color=col, label=label, lw=1.8)
    for a, title in zip(ax, ("Extracellular", "Intracellular")):
        a.set_xlabel("time [h]")
        a.set_title(title)
        a.grid(alpha=0.25)
    ax[0].axhline(1e2, ls="--", c="k", lw=1, label=r"$10^2$ CFU/mL")
    ax[0].set_ylabel("CFU/mL")
    ax[0].legend(fontsize=8, loc="lower left")
    fig.suptitle(
        f"Standard patient, MOI {REF_MOI} particles/CFU, "
        f"{DOSE_MAX_G:.2f} g q12h, $\\tau_P$ = {REF_TAU_P:.0f} h", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "timecourse.png", dpi=160)
    plt.close(fig)


def fig_adhesion_tpd() -> None:
    p = ref()
    s = solve_course(p, n_out=N_OUT)
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.axhspan(p.window_lo, p.window_hi, color="#cfe8df", alpha=0.7,
               label="effective window 93-1991")
    ax.semilogy(s.t, np.maximum(s.A_cell, 1e-2), color=C_COMBO, lw=1.8,
                label=r"$A_{\rm cell}$ [molecules/cell]")
    ax2 = ax.twinx()
    ax2.semilogy(s.t, s.m, color=C_ABX, lw=1.4, ls="--",
                 label=r"$m$ [particles/cell]")
    ax.set_xlabel("time [h]")
    ax.set_ylabel("TPD per cell")
    ax2.set_ylabel("adhered particles per cell")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="lower right")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "adhesion_tpd.png", dpi=160)
    plt.close(fig)


def fig_rho_mic() -> None:
    p = ref()
    s = solve_course(p, n_out=N_OUT)
    last = s.t >= p.t_end - p.tau_A
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.plot(s.t[last], s.rho[last], color=C_COMBO, lw=1.8, label=r"$\rho$")
    ax2 = ax.twinx()
    ax2.plot(s.t[last], s.MIC[last], color=C_ABX, lw=1.8, ls="--", label="MIC")
    ax2.plot(s.t[last], s.A[last], color="k", lw=1.0, ls=":", label="A(t)")
    ax.set_xlabel("time [h]")
    ax.set_ylabel(r"resistance-factor survival $\rho$")
    ax2.set_ylabel(r"concentration [$\mu$g/mL]")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8)
    ax.grid(alpha=0.25)
    ax.set_title("Final antibiotic cycle", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "rho_mic_cycle.png", dpi=160)
    plt.close(fig)


def fig_moi_sweep() -> None:
    mois = np.geomspace(0.1, 30.0, 26)
    rows = []
    for m in mois:
        v = evaluate(solve_course(ref().with_MOI(float(m)), n_out=N_OUT))
        rows.append((m, v.delta_log, v.C_bar))
    m, dl, cb = map(np.array, zip(*rows))
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.0))
    ax[0].semilogx(m, dl, color=C_COMBO, lw=1.8)
    ax[0].axhline(3.0, ls="--", c="k", lw=1)
    ax[0].set_ylabel(r"$\Delta\log_{10}$")
    ax[1].loglog(m, cb, color=C_COMBO, lw=1.8)
    ax[1].axhline(1e2, ls="--", c="k", lw=1)
    ax[1].set_ylabel(r"$\bar C^{\rm tot}$ [CFU/mL]")
    for a in ax:
        a.set_xlabel("MOI [particles/CFU]")
        a.grid(alpha=0.25, which="both")
    fig.suptitle("Judgement against MOI at the safety-capped dose", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "moi_sweep.png", dpi=160)
    plt.close(fig)


def fig_dose_sweep() -> None:
    doses = np.linspace(0.1, 1.3, 25)
    combo, mono = [], []
    for g in doses:
        combo.append(evaluate(solve_course(
            ref().with_dose_g(float(g)), n_out=N_OUT)).delta_log)
        mono.append(evaluate(solve_course(
            ref(sunrise_on=False).with_dose_g(float(g)), n_out=N_OUT)).delta_log)
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.plot(doses, combo, color=C_COMBO, lw=1.8, label="with SUNRISE")
    ax.plot(doses, mono, color=C_ABX, lw=1.8, ls="--", label="kanamycin alone")
    ax.axhline(3.0, ls="--", c="k", lw=1, label=r"$\Delta\log_{10}=3$")
    ax.axvline(DOSE_MAX_G, color="#b2182b", lw=1.4,
               label=f"safety cap {DOSE_MAX_G:.2f} g")
    ax.set_xlabel("kanamycin per dose [g]  (q12h, 60 kg IM)")
    ax.set_ylabel(r"$\Delta\log_{10}$")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG / "dose_sweep.png", dpi=160)
    plt.close(fig)


def fig_tau_dose_heatmap() -> None:
    path = RESULTS / "task2_grid.csv"
    if not path.exists():
        return
    import csv

    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    taus = sorted({float(r["tau_P"]) for r in rows})
    doses = sorted({float(r["dose_g"]) for r in rows})
    Z = np.full((len(doses), len(taus)), np.nan)
    for r in rows:
        i = doses.index(float(r["dose_g"]))
        j = taus.index(float(r["tau_P"]))
        Z[i, j] = float(r["delta_log"])
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    im = ax.imshow(Z, origin="lower", aspect="auto", cmap="viridis",
                   extent=(-0.5, len(taus) - 0.5, -0.5, len(doses) - 0.5))
    ax.set_xticks(range(len(taus)), [f"{t:g}" for t in taus])
    ax.set_yticks(range(len(doses)), [f"{d:.2f}" for d in doses])
    ax.set_xlabel(r"$\tau_P$ [h]")
    ax.set_ylabel("kanamycin per dose [g]")
    for r in rows:
        i = doses.index(float(r["dose_g"]))
        j = taus.index(float(r["tau_P"]))
        ok = r["passes"] == "True" and r["safe"] == "True"
        ax.text(j, i, "o" if ok else "x", ha="center", va="center",
                color="w" if ok else "#ff9999", fontsize=9)
    fig.colorbar(im, ax=ax, label=r"$\Delta\log_{10}$")
    ax.set_title(f"MOI fixed at 1.297; o = passes within safety", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "tau_dose_grid.png", dpi=160)
    plt.close(fig)


def fig_zero_class() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    for moi, col in ((0.5, C_ABX), (1.409, C_COMBO), (3.2, C_SUNRISE),
                     (9.5, "#444444")):
        s = solve_course(ref().with_MOI(moi), n_out=N_OUT)
        ax.semilogy(s.t, np.maximum(s.zero_class, 1e-12), lw=1.6, color=col,
                    label=f"MOI {moi}")
    ax.axhline(0.01, ls="--", c="k", lw=1, label="1 %")
    ax.set_xlabel("time [h]")
    ax.set_ylabel(r"$e^{-\bar\varphi m(t)}$")
    ax.set_title("Fraction receiving no TPD molecule in a dose", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    fig.savefig(FIG / "zero_class.png", dpi=160)
    plt.close(fig)


def fig_tolerance() -> None:
    path = RESULTS / "task4.json"
    if not path.exists():
        return
    data = json.loads(path.read_text(encoding="utf-8"))
    axes = [a for a in data["axes"]
            if a["lowest_passing"] or a["highest_passing"]]
    names = [a["axis"] for a in axes]
    fig, ax = plt.subplots(figsize=(8.2, 0.42 * len(names) + 1.6))
    for i, a in enumerate(axes):
        lo = a["lowest_passing"] or a["search_lo"]
        hi = a["highest_passing"] or a["search_hi"]
        ax.plot([lo, hi], [i, i], lw=5, color=C_SUNRISE, solid_capstyle="butt")
        ax.plot([lo], [i], "|", color="#b2182b", ms=11)
        ax.plot([hi], [i], "|", color="#b2182b", ms=11)
    ax.set_yticks(range(len(names)), names, fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("passing range (each axis moved alone)")
    ax.grid(alpha=0.25, axis="x", which="both")
    fig.tight_layout()
    fig.savefig(FIG / "tolerance.png", dpi=160)
    plt.close(fig)


def main() -> None:
    fig_timecourse()
    fig_adhesion_tpd()
    fig_rho_mic()
    fig_moi_sweep()
    fig_dose_sweep()
    fig_tau_dose_heatmap()
    fig_zero_class()
    fig_tolerance()
    print("figures written to", FIG)
    for f in sorted(FIG.glob("*.png")):
        print("  ", f.name, f.stat().st_size // 1024, "KB")


if __name__ == "__main__":
    main()
