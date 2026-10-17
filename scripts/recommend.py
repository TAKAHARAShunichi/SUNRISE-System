"""The recommended design: a margin on durability, not the bare minimum.

The minimum-MOI designs sit on the C_bar = 100 CFU/mL boundary by
construction, so a 1 h phase error or any adverse axis moves them across it.
This script reports the MOI that keeps C_bar at or below a target fraction of
the threshold at the worst phase.
"""

from __future__ import annotations

import json
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from sunrise.design import N_OUT, bisect_min_log, moi_cap
from sunrise.metrics import evaluate
from sunrise.model import m_infinity, solve_course
from sunrise.params import DOSE_MAX_G, STANDARD_PATIENT

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

N_PHASE = 8
#: C_bar target as a fraction of the 10^2 CFU/mL threshold.
#: The floor of C_bar is K_in (= 50 CFU/mL), so 0.5 is unreachable by
#: construction; 0.7 is the tightest target that leaves real headroom.
MARGIN = 0.7


def worst_C_bar(p) -> tuple[float, float, bool]:
    worst_cb, worst_dl, ok = -1.0, None, True
    for dt in np.linspace(0.0, p.tau_P, N_PHASE, endpoint=False):
        v = evaluate(solve_course(p.with_(dt_phase=float(dt)), n_out=N_OUT))
        if v.C_bar > worst_cb:
            worst_cb, worst_dl = v.C_bar, v.delta_log
        ok = ok and v.passes
    return worst_cb, worst_dl, ok


def _one(tau_P: float) -> dict:
    base = STANDARD_PATIENT.with_(tau_P=tau_P).with_dose_g(DOSE_MAX_G)
    target = STANDARD_PATIENT.cure_density * MARGIN

    def ok(moi: float) -> bool:
        cb, _, passes_all = worst_C_bar(base.with_MOI(moi))
        return passes_all and cb <= target

    moi = bisect_min_log(ok, 0.05, 3e3, tol=2e-3)
    out = {"tau_P": tau_P, "n_doses": base.n_doses_P,
           "moi_recommended": moi, "moi_cap": moi_cap(base),
           "target_C_bar": target}
    if moi:
        p = base.with_MOI(moi)
        cb, dl, _ = worst_C_bar(p)
        out.update({
            "D_per_dose": p.D,
            "particles_per_course": p.D * p.n_doses_P,
            "m_inf": m_infinity(p),
            "worst_C_bar": cb,
            "worst_delta_log": dl,
            "A_cell_max": float(solve_course(p, n_out=N_OUT).A_cell.max()),
            "within_window": bool(
                solve_course(p, n_out=N_OUT).A_cell.max() <= p.window_hi),
            "compatible_with_cap": bool(out["moi_cap"] and moi <= out["moi_cap"]),
        })
    return out


def _margin_curve_one(moi: float) -> dict:
    base = STANDARD_PATIENT.with_(tau_P=6.0).with_dose_g(DOSE_MAX_G)
    cb, dl, ok = worst_C_bar(base.with_MOI(moi))
    p = base.with_MOI(moi)
    return {
        "MOI": moi,
        "worst_C_bar": cb,
        "worst_delta_log": dl,
        "passes_all_phases": ok,
        "margin_fraction": cb / STANDARD_PATIENT.cure_density,
        "A_cell_max": float(solve_course(p, n_out=N_OUT).A_cell.max()),
        "within_window": bool(
            solve_course(p, n_out=N_OUT).A_cell.max() <= p.window_hi),
        "particles_per_course": p.D * p.n_doses_P,
    }


def main() -> None:
    with Pool(4) as pool:
        rows = pool.map(_one, (3.0, 4.0, 6.0, 8.0))
        curve = pool.map(
            _margin_curve_one,
            [float(x) for x in np.round(np.geomspace(1.4, 60.0, 9), 3)],
        )
    (RESULTS / "recommended.json").write_text(
        json.dumps({"margin": MARGIN, "rows": rows,
                    "floor_note": (
                        "C_bar cannot fall below K_in = 50 CFU/mL, so the "
                        "durability margin against the 10^2 threshold is "
                        "structurally capped at 2x."),
                    "margin_curve_tau_P_6h": curve},
                   indent=2, ensure_ascii=False), encoding="utf-8")
    print("\nworst-phase C_bar against MOI (tau_P = 6 h, dose at the cap)")
    print(f"{'MOI':>8} {'C_bar':>8} {'frac':>6} {'dlog':>6} {'A_max':>9} {'win':>5}")
    for c in curve:
        print(f"{c['MOI']:8.2f} {c['worst_C_bar']:8.2f} "
              f"{c['margin_fraction']:6.2f} {c['worst_delta_log']:6.2f} "
              f"{c['A_cell_max']:9.0f} {str(c['within_window']):>5}")
    print(f"{'tau_P':>6} {'n':>3} {'MOI':>8} {'MOI_cap':>8} {'D/dose':>11} "
          f"{'particles/course':>17} {'m_inf':>7} {'A_max':>8} {'win':>5} "
          f"{'Cbar':>7}")
    for r in rows:
        if not r.get("moi_recommended"):
            print(f"{r['tau_P']:6.1f} no solution")
            continue
        print(f"{r['tau_P']:6.1f} {r['n_doses']:3d} {r['moi_recommended']:8.3f} "
              f"{r['moi_cap']:8.3f} {r['D_per_dose']:11.3e} "
              f"{r['particles_per_course']:17.3e} {r['m_inf']:7.1f} "
              f"{r['A_cell_max']:8.1f} {str(r['within_window']):>5} "
              f"{r['worst_C_bar']:7.2f}")


if __name__ == "__main__":
    main()
