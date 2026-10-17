"""Competitive release scenario (Ver.23: run separately from the main line).

The standard patient is 100 % resistant, so competitive release cannot occur.
This scenario is the other claim: a susceptible-dominant infection in which a
resistant minority is present, where the point is not curing a resistant
infection but preventing one from being selected for.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from sunrise.design import N_OUT
from sunrise.metrics import evaluate
from sunrise.model import mic_of_rho, solve_course
from sunrise.params import COMPETITIVE_RELEASE, DOSE_MAX_G

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

REF_TAU_P, REF_MOI, REF_PHASE = 6.0, 3.293, 3.0


def selective_window_time(sol) -> float:
    """Hours spent with MIC_S < A(t) < MIC(t): the selective window."""
    A = sol.A
    mic = sol.MIC
    inside = (A > sol.p.MIC_S) & (A < mic)
    dt = sol.t[1] - sol.t[0]
    return float(inside.sum() * dt)


def arm(name: str, p) -> dict:
    sol = solve_course(p, n_out=N_OUT)
    v = evaluate(sol)
    CS, CR = sol.y[0], sol.y[1]
    return {
        "arm": name,
        "dose_g": p.dose_g if p.antibiotic_on else 0.0,
        "MOI": p.MOI if p.sunrise_on else 0.0,
        "delta_log": v.delta_log,
        "C_bar": v.C_bar,
        "passes": bool(v.passes),
        "CR_ex_initial": float(CR[0]),
        "CR_ex_final": float(CR[-1]),
        "CR_ex_log_change": float(np.log10(max(CR[-1], 1e-30) / max(CR[0], 1e-30))),
        "CR_ex_peak": float(CR.max()),
        "resistant_fraction_initial": float(CR[0] / (CS[0] + CR[0])),
        "resistant_fraction_final": float(
            CR[-1] / max(CS[-1] + CR[-1], 1e-30)),
        "selective_window_hours": selective_window_time(sol),
    }


def main() -> None:
    base = COMPETITIVE_RELEASE.with_(
        tau_P=REF_TAU_P, dt_phase=REF_PHASE).with_dose_g(DOSE_MAX_G)
    arms = [
        arm("untreated", base.with_(
            sunrise_on=False, antibiotic_on=False)),
        arm("kanamycin only", base.with_(sunrise_on=False)),
        arm("kanamycin + SUNRISE", base.with_MOI(REF_MOI)),
    ]
    out = {
        "scenario": {
            "S0": base.S0, "R0": base.R0,
            "resistant_fraction": base.R0 / (base.S0 + base.R0),
            "tau_P": REF_TAU_P, "MOI": REF_MOI, "dose_g": DOSE_MAX_G,
            "gamma_conj": base.gamma_conj,
        },
        "arms": arms,
        "note": (
            "The claim here is prevention, not cure: kanamycin alone clears the "
            "susceptible majority and frees carrying capacity, so the resistant "
            "minority expands. Adding SUNRISE lowers MIC(t) so the same dose "
            "also reaches the resistant subpopulation."
        ),
    }
    (RESULTS / "competitive_release.json").write_text(
        json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"{'arm':<22}{'dlog':>6}{'Cbar':>11}{'pass':>6}"
          f"{'CR_ex 0':>11}{'CR_ex 72h':>12}{'dlogCR':>8}{'R frac':>9}{'win h':>7}")
    for a in arms:
        print(f"{a['arm']:<22}{a['delta_log']:6.2f}{a['C_bar']:11.3g}"
              f"{str(a['passes']):>6}{a['CR_ex_initial']:11.2e}"
              f"{a['CR_ex_final']:12.2e}{a['CR_ex_log_change']:+8.2f}"
              f"{a['resistant_fraction_final']:9.3g}"
              f"{a['selective_window_hours']:7.1f}")


if __name__ == "__main__":
    main()
