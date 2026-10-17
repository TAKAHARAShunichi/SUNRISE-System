"""Worst-case-phase design: the MOI that passes at every phase offset.

Ver.22 Table 9.4.30 took the minimum D that passes at the worst phase; the
phase is a design variable in Ver.23, but a design that only passes at one
phase is not implementable, so the worst-case figure is the one to recommend.
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
TAU_P_GRID = (2.0, 3.0, 4.0, 6.0, 8.0, 12.0)


def worst_phase_verdict(p) -> dict:
    """Judgement at the worst phase offset in [0, tau_P)."""
    worst = None
    for dt in np.linspace(0.0, p.tau_P, N_PHASE, endpoint=False):
        v = evaluate(solve_course(p.with_(dt_phase=float(dt)), n_out=N_OUT))
        key = (v.passes, v.delta_log, -v.C_bar)
        if worst is None or key < worst[0]:
            worst = (key, float(dt), v)
    (ok, dlog, _), dt, v = worst
    return {"dt_worst": dt, "delta_log": v.delta_log, "C_bar": v.C_bar,
            "passes": bool(v.passes)}


def _one(tau_P: float) -> dict:
    base = STANDARD_PATIENT.with_(tau_P=tau_P).with_dose_g(DOSE_MAX_G)
    moi = bisect_min_log(
        lambda m: worst_phase_verdict(base.with_MOI(m))["passes"], 0.05, 3e3,
        tol=2e-3,
    )
    out = {"tau_P": tau_P, "n_doses": base.n_doses_P,
           "moi_worst_phase": moi, "moi_cap": moi_cap(base)}
    if moi is not None:
        p = base.with_MOI(moi)
        out.update(
            {
                "D_per_dose": p.D,
                "particles_per_course": p.D * p.n_doses_P,
                "m_inf": m_infinity(p),
                "compatible": bool(out["moi_cap"] and moi <= out["moi_cap"]),
                **{f"worst_{k}": v for k, v in
                   worst_phase_verdict(p).items()},
            }
        )
    return out


def main() -> None:
    with Pool(4) as pool:
        rows = pool.map(_one, TAU_P_GRID)
    (RESULTS / "task1_worst_phase.json").write_text(
        json.dumps(rows, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"{'tau_P':>6} {'n':>3} {'MOI(worst)':>11} {'MOI_cap':>9} "
          f"{'D/dose':>11} {'particles/course':>17} {'m_inf':>8} {'dlog':>6} {'Cbar':>8}")
    for r in rows:
        if r.get("moi_worst_phase") is None:
            print(f"{r['tau_P']:6.1f} {r['n_doses']:3d}  no solution")
            continue
        print(
            f"{r['tau_P']:6.1f} {r['n_doses']:3d} {r['moi_worst_phase']:11.3f} "
            f"{r['moi_cap']:9.3f} {r['D_per_dose']:11.3e} "
            f"{r['particles_per_course']:17.3e} {r['m_inf']:8.1f} "
            f"{r['worst_delta_log']:6.2f} {r['worst_C_bar']:8.1f}"
        )


if __name__ == "__main__":
    main()
