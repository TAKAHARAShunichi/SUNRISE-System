"""Tasks 1-3: required MOI, the (tau_P, dose) search, and the single-agent arm.

Writes results/task1.json, results/task2.csv, results/task3.json.
"""

from __future__ import annotations

import json
import math
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from sunrise.design import (
    N_OUT,
    a_cell_max,
    bisect_min_log,
    min_moi_at_dose,
    moi_cap,
    objective,
    passes,
)
from sunrise.metrics import evaluate, summarise
from sunrise.model import m_infinity, solve_course
from sunrise.params import (
    DAILY_TOTAL_MAX_G,
    DOSE_MAX_G,
    STANDARD_PATIENT,
    Params,
    dose_to_serum,
    dose_to_urine,
)

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

TAU_P_GRID = (2.0, 3.0, 4.0, 6.0, 8.0, 12.0)
DOSE_GRID = (0.20, 0.30, 0.40, 0.45, 0.55, 0.65, DOSE_MAX_G)
N_PHASE = 8


# --- task 1 -----------------------------------------------------------------


def _task1_one(tau_P: float) -> dict:
    base = STANDARD_PATIENT.with_(tau_P=tau_P)
    cap = moi_cap(base)
    safe = min_moi_at_dose(base, DOSE_MAX_G)
    at_ref = min_moi_at_dose(base, 0.45)
    # "any dose": the safety box is lifted, so the ceiling is generous
    any_dose = min_moi_at_dose(base, 6.0)
    return {
        "tau_P": tau_P,
        "n_doses": base.n_doses_P,
        "moi_min_any_dose": any_dose,
        "moi_min_at_baseline_dose_0.45g": at_ref,
        "moi_optimal_safe": safe,
        "moi_cap": cap,
        "compatible": bool(safe is not None and cap is not None and safe <= cap),
        "A_cell_max_at_optimal": a_cell_max(base.with_MOI(safe)) if safe else None,
        "particles_per_dose_at_optimal": safe * base.total_CFU0 if safe else None,
        "particles_per_course_at_optimal": (
            safe * base.total_CFU0 * base.n_doses_P if safe else None
        ),
        "m_inf_at_optimal": m_infinity(base.with_MOI(safe)) if safe else None,
    }


def task1(pool: Pool) -> dict:
    rows = pool.map(_task1_one, TAU_P_GRID)
    compatible = [r["tau_P"] for r in rows if r["compatible"]]
    return {
        "rows": rows,
        "largest_compatible_tau_P": max(compatible) if compatible else None,
        "note": (
            "moi_min_any_dose lifts the safety box; the antibiotic alone passes "
            "at 1.05 g/dose, so this column reports the floor of the MOI search "
            "range rather than a SUNRISE requirement."
        ),
    }


# --- task 2 -----------------------------------------------------------------


def _task2_one(arg: tuple[float, float, float]) -> dict:
    moi, tau_P, dose_g = arg
    base = STANDARD_PATIENT.with_(tau_P=tau_P).with_dose_g(dose_g).with_MOI(moi)
    best = None
    for dt in np.linspace(0.0, tau_P, N_PHASE, endpoint=False):
        p = base.with_(dt_phase=float(dt))
        v = evaluate(solve_course(p, n_out=N_OUT))
        key = (v.passes, v.delta_log)
        if best is None or key > best[0]:
            best = (key, float(dt), v)
    (ok, dlog), dt_best, v = best
    p_best = base.with_(dt_phase=dt_best)
    return {
        "MOI": moi,
        "tau_P": tau_P,
        "dose_g": dose_g,
        "daily_g": p_best.daily_g,
        "serum_peak": p_best.serum_peak,
        "A_peak_urine": p_best.A_peak,
        "dt_opt": dt_best,
        "delta_log": v.delta_log,
        "C_bar": v.C_bar,
        "passes": bool(v.passes),
        "n_doses": p_best.n_doses_P,
        "D": p_best.D,
        "particles_per_course": p_best.D * p_best.n_doses_P,
        "m_inf": m_infinity(p_best),
        "safe": bool(p_best.serum_peak <= 35.0 and p_best.daily_g <= DAILY_TOTAL_MAX_G),
    }


def _phase_sensitivity(moi: float, tau_P: float, dose_g: float, dt: float) -> dict:
    base = STANDARD_PATIENT.with_(tau_P=tau_P).with_dose_g(dose_g).with_MOI(moi)
    out = {}
    for label, shift in (("minus1h", -1.0), ("plus1h", +1.0)):
        d = (dt + shift) % tau_P
        v = evaluate(solve_course(base.with_(dt_phase=d), n_out=N_OUT))
        out[label] = {"dt": d, "delta_log": v.delta_log, "passes": bool(v.passes)}
    return out


def task2(pool: Pool, moi_fixed: float, m_inf_ref: float, D_tot_ref: float) -> dict:
    args = [(moi_fixed, t, d) for t in TAU_P_GRID for d in DOSE_GRID]
    rows = pool.map(_task2_one, args)
    for r in rows:
        p = (
            STANDARD_PATIENT.with_(tau_P=r["tau_P"], dt_phase=r["dt_opt"])
            .with_dose_g(r["dose_g"])
            .with_MOI(r["MOI"])
        )
        r["J"] = objective(p, m_inf_ref, D_tot_ref)
    feasible = [r for r in rows if r["passes"] and r["safe"]]
    best = min(feasible, key=lambda r: r["J"]) if feasible else None
    out = {"grid": rows, "best": best}
    if best:
        out["phase_sensitivity"] = _phase_sensitivity(
            best["MOI"], best["tau_P"], best["dose_g"], best["dt_opt"]
        )
    return out


# --- task 2b: three-variable simultaneous search ---------------------------


def _task2b_one(arg: tuple[float, float, float]) -> dict:
    return _task2_one(arg)


def task2_joint(pool: Pool, m_inf_ref: float, D_tot_ref: float,
                moi_grid: tuple[float, ...]) -> dict:
    args = [(m, t, d) for m in moi_grid for t in TAU_P_GRID for d in DOSE_GRID]
    rows = pool.map(_task2b_one, args)
    for r in rows:
        p = (
            STANDARD_PATIENT.with_(tau_P=r["tau_P"], dt_phase=r["dt_opt"])
            .with_dose_g(r["dose_g"])
            .with_MOI(r["MOI"])
        )
        r["J"] = objective(p, m_inf_ref, D_tot_ref)
    feasible = [r for r in rows if r["passes"] and r["safe"]]
    best = min(feasible, key=lambda r: r["J"]) if feasible else None
    return {"best": best, "n_feasible": len(feasible), "n_total": len(rows)}


# --- task 3 -----------------------------------------------------------------


def task3() -> dict:
    mono = STANDARD_PATIENT.with_(sunrise_on=False)
    d = bisect_min_log(
        lambda g: evaluate(solve_course(mono.with_dose_g(g), n_out=N_OUT)).passes,
        0.02, 20.0, tol=1e-4,
    )
    arms = {}
    for name, p in (
        ("untreated", STANDARD_PATIENT.with_(sunrise_on=False, antibiotic_on=False)),
        ("antibiotic_only_baseline", mono),
        ("antibiotic_only_safety_cap", mono.with_dose_g(DOSE_MAX_G)),
    ):
        v = evaluate(solve_course(p, n_out=N_OUT))
        arms[name] = {
            "dose_g": p.dose_g if p.antibiotic_on else 0.0,
            "delta_log": v.delta_log,
            "C_bar": v.C_bar,
            "passes": bool(v.passes),
        }
    return {
        "min_dose_g_per_dose": d,
        "min_daily_g": d * 2.0,
        "min_urine_ug_per_mL": dose_to_urine(d),
        "min_serum_peak_ug_per_mL": dose_to_serum(d),
        "safety_cap_dose_g": DOSE_MAX_G,
        "safety_cap_daily_g": DAILY_TOTAL_MAX_G,
        "exceeds_per_dose_by": d / DOSE_MAX_G,
        "exceeds_daily_by": (d * 2.0) / DAILY_TOTAL_MAX_G,
        "feasible_within_safety": bool(d <= DOSE_MAX_G),
        "arms": arms,
    }


def main() -> None:
    with Pool(4) as pool:
        t1 = task1(pool)
        row6 = next(r for r in t1["rows"] if r["tau_P"] == 6.0)
        moi_fixed = row6["moi_optimal_safe"]
        ref = STANDARD_PATIENT.with_MOI(moi_fixed)
        m_inf_ref = m_infinity(ref)
        D_tot_ref = ref.D * ref.n_doses_P

        t2 = task2(pool, moi_fixed, m_inf_ref, D_tot_ref)
        moi_grid = tuple(
            float(x) for x in np.round(np.geomspace(0.3, 9.0, 7), 4)
        )
        t2b = task2_joint(pool, m_inf_ref, D_tot_ref, moi_grid)
    t3 = task3()

    (RESULTS / "task1.json").write_text(
        json.dumps(t1, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    (RESULTS / "task2.json").write_text(
        json.dumps({"fixed_MOI": moi_fixed, "normalisation": {
            "m_inf_ref": m_inf_ref, "D_tot_ref": D_tot_ref},
            "sequential": t2, "joint": t2b}, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    (RESULTS / "task3.json").write_text(
        json.dumps(t3, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    import csv

    with (RESULTS / "task2_grid.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(t2["grid"][0].keys()))
        w.writeheader()
        w.writerows(t2["grid"])

    print("=== TASK 1 ===")
    for r in t1["rows"]:
        print(
            f"tau_P={r['tau_P']:5.1f} n={r['n_doses']:3d} "
            f"opt_safe={r['moi_optimal_safe']:8.3f} cap={r['moi_cap']:8.3f} "
            f"compat={r['compatible']} "
            f"D/dose={r['particles_per_dose_at_optimal']:10.3e} "
            f"m_inf={r['m_inf_at_optimal']:8.1f}"
        )
    print("largest compatible tau_P:", t1["largest_compatible_tau_P"])

    print("\n=== TASK 2 (MOI fixed at %.3f) ===" % moi_fixed)
    b = t2["best"]
    print(
        f"best: tau_P={b['tau_P']} dose={b['dose_g']:.3f} g "
        f"daily={b['daily_g']:.3f} g serum={b['serum_peak']:.1f} "
        f"dt_opt={b['dt_opt']:.2f} h J={b['J']:.4f} dlog={b['delta_log']:.2f} "
        f"Cbar={b['C_bar']:.1f}"
    )
    print("phase +-1h:", t2["phase_sensitivity"])
    jb = t2b["best"]
    print(
        f"joint best: MOI={jb['MOI']} tau_P={jb['tau_P']} dose={jb['dose_g']:.3f} "
        f"J={jb['J']:.4f}  (sequential J={b['J']:.4f}, "
        f"gap={100*(b['J']-jb['J'])/jb['J']:+.2f} %)"
    )

    print("\n=== TASK 3 ===")
    print(
        f"single agent minimum: {t3['min_dose_g_per_dose']:.4f} g/dose, "
        f"{t3['min_daily_g']:.3f} g/day, urine {t3['min_urine_ug_per_mL']:.1f}, "
        f"serum {t3['min_serum_peak_ug_per_mL']:.2f} ug/mL -> "
        f"{t3['exceeds_daily_by']:.2f}x the daily cap"
    )


if __name__ == "__main__":
    main()
