#!/usr/bin/env python3
"""Task 1: re-seat the baseline CPI dose D against the upper window bound.

D = 3.0e11 particles/dose puts the TPD load per bacterium above the upper end
of the effective window (1991 molecules/CFU). Only the upper end binds, so we
take the largest D whose COURSE-WIDE maximum A_cell stays at or below 1991.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

import acell
from sunrise.model import cpi_attachment, dose_times
from sunrise import metrics

TARGET = acell.WINDOW_HI
VOLUME_ML = 200.0


def course_max(dose: float, p: dict) -> float:
    return acell.acell_stats(dict(p, dose=dose))["max"]


def solve_dose(p: dict, lo: float = 1.0e11, hi: float = 3.0e11) -> float:
    f = lambda D: course_max(D, p) - TARGET
    return float(brentq(f, lo, hi, xtol=1.0e8, rtol=1e-6))


def report(p: dict, dose: float, label: str) -> dict:
    p = dict(p, dose=dose)
    st = acell.acell_stats(p)
    att = cpi_attachment(p)
    ev = metrics.evaluate(p)
    n_doses = len(dose_times(p))
    # "final cycle" in the brief is the last ANTIBIOTIC cycle (tau_a = 12 h),
    # which is what reproduces its quoted 120-2257 / 18.8x / 89.7%.
    import numpy as _np
    from sunrise.model import simulate
    t, _, _, state = simulate(p, points_per_hour=200)
    ac = state[acell.A_CELL]
    m = t >= p["t_end"] - p["tau_a"]
    occ_abx = acell._occupancy(t[m], ac[m], acell.WINDOW_LO, acell.WINDOW_HI)

    out = {
        "label": label, "dose": dose, "n_doses": n_doses,
        "per_course": dose * n_doses,
        "k_enc": p["k_enc"], "f_hit0": att["f_hit"], "m0": att["m0"],
        "m_inf": att["m_inf"], "ratio": att["ratio"], "N0_0": att["tpd_load"],
        "course_max": st["max"], "t_max": st["t_max"],
        "abx_min": float(ac[m].min()), "abx_max": float(ac[m].max()),
        "abx_swing": float(ac[m].max() / ac[m].min()),
        "occ_abx": occ_abx, "occ_course": st["occupancy_course"],
        "depth": ev["depth_log"], "final_total": ev["final_total"],
        "pass_A": ev["pass_depth"], "pass_B": ev["pass_durability"],
    }
    return out


def show(r: dict) -> None:
    print(f"--- {r['label']}")
    print(f"  D                       {r['dose']:.4g} particles/dose"
          f"   ({r['n_doses']} doses, {r['per_course']:.4g} particles/course)")
    print(f"  k_enc                   {r['k_enc']:.4g} mL/h")
    print(f"  f_hit(0)                {r['f_hit0']:.4f}")
    print(f"  m(0)                    {r['m0']:.2f} particles/CFU")
    print(f"  m_inf                   {r['m_inf']:.1f} particles/CFU"
          f"   (accumulation {r['ratio']:.1f}x)")
    print(f"  N_0(0)                  {r['N0_0']:.1f} molecules/CFU")
    print(f"  A_cell course max       {r['course_max']:.1f} at t = {r['t_max']:.2f} h"
          f"   (window hi {TARGET:.0f})")
    print(f"  A_cell final abx cycle  {r['abx_min']:.1f} - {r['abx_max']:.1f}"
          f"   (swing {r['abx_swing']:.1f}x)")
    print(f"  in-window, final cycle  {100*r['occ_abx']:.1f} %")
    print(f"  in-window, whole course {100*r['occ_course']:.1f} %")
    print(f"  judgement A  dlog10     {r['depth']:.2f}   "
          f"{'PASS' if r['pass_A'] else 'FAIL'} (>= 3)")
    print(f"  judgement B  C_tot      {r['final_total']:.2f} CFU/mL "
          f"= {r['final_total']*VOLUME_ML:.4g} CFU/bladder   "
          f"{'PASS' if r['pass_B'] else 'FAIL'} (< 100 CFU/mL)")
    print()


def main() -> None:
    print("Ver.22 basis; k_tpd0 calibrated at the published conditions.\n")
    p = acell.base_params()
    print(f"  calibrated k_tpd0 = {p['k_tpd0']/60:.5f} /min")
    print(f"  phi_p = {p['phi_p']}, eta_mech = {p['eta_mech']}, "
          f"mic_r0 = {p['mic_r0']}, psi_in = {p['psi_in']}\n")

    print("Where is the A_cell maximum? (per-CPI-cycle maxima, D = 3.0e11)")
    from sunrise.model import simulate
    t, _, _, state = simulate(dict(p, dose=3.0e11), points_per_hour=200)
    ac = state[acell.A_CELL]
    for j in range(int(p["t_end"] // p["tau_p"])):
        m = (t >= j * p["tau_p"]) & (t < (j + 1) * p["tau_p"])
        if m.any():
            print(f"    cycle {j+1:2d}  [{j*6:2.0f}-{(j+1)*6:2.0f} h]  "
                  f"max {ac[m].max():7.1f}")
    print("  -> the course maximum is NOT in the final cycle.\n")

    print("Monotonicity of the course maximum in D")
    grid = [1.0e11, 1.5e11, 2.0e11, 2.5e11, 2.6e11, 3.0e11]
    vals = [course_max(D, p) for D in grid]
    for D, v in zip(grid, vals):
        print(f"    D = {D:.2e}  ->  max A_cell = {v:8.1f}")
    assert all(b > a for a, b in zip(vals, vals[1:])), "not monotone"
    print("  -> monotone increasing; bisection is well posed.\n")

    D = solve_dose(p)
    print(f"Bisection result:  D = {D:.4g} particles/dose\n")
    show(report(p, 3.0e11, "current baseline D = 3.0e11 (for comparison)"))
    show(report(p, D, f"NEW baseline D = {D:.4g}"))

    # psi_in insensitivity of the A_cell constraint
    D_old_psi = solve_dose(dict(p, psi_in=0.20))
    print(f"Sensitivity of D to psi_in: {D:.4g} (psi_in=1.3) vs "
          f"{D_old_psi:.4g} (psi_in=0.20)  -> {abs(D/D_old_psi-1)*100:.2f}% apart")
    print("  The upper-window constraint is set by TPD dynamics, which barely")
    print("  depend on the intracellular growth rate, so task 2 does not move D.")


if __name__ == "__main__":
    main()
