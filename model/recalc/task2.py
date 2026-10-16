#!/usr/bin/env python3
"""Task 2: the intracellular compartment fixed as the acute-phase IBC.

psi_in = 1.3 /h (median of the measured 1.19-1.39), swept over 0.69-1.39.
Run after task 1 has fixed D.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

import acell
from sunrise import analysis, metrics

VOLUME_ML = 200.0
D_NEW = 2.571e11
CURE = 100.0            # judgement B threshold, CFU/mL


def closed_form(p: dict) -> dict:
    """d and Lambda from the parameters, with kill_in from the simulation."""
    g = analysis.growth_rates(p)
    a_bar, d = g["a_bar"], g["d"]
    tr = a_bar + d
    det = a_bar * d - p["k_emerge"] * p["k_invade"]
    lam = 0.5 * (tr + np.sqrt(tr * tr - 4 * det))
    return {**g, "trace": tr, "det": det, "lambda": lam}


def main() -> None:
    p = acell.base_params(dose=D_NEW)
    print(f"Ver.22 basis, D = {D_NEW:.4g} particles/dose, "
          f"psi_in = {p['psi_in']} /h\n")

    print("=" * 72)
    print("1. Table 9.4.23  net rates and dominant eigenvalue")
    print("=" * 72)
    cf = closed_form(p)
    fl = analysis.floor_level(p, measure=True)
    print(f"  extracellular a_bar            {cf['a_bar']:+.4f} /h")
    print(f"  intracellular kill  k_A^in     {cf['kill_in']:+.4f} /h")
    print(f"  intracellular net   d          {cf['d']:+.4f} /h")
    print(f"  trace                          {cf['trace']:+.4f} /h")
    print(f"  determinant                    {cf['det']:+.6f} /h^2")
    print(f"  dominant eigenvalue Lambda     {cf['lambda']:+.4f} /h"
          f"   (doubling {np.log(2)/cf['lambda']:.3f} h)")
    print(f"  Lambda - d                     {cf['lambda']-cf['d']:+.2e} /h"
          f"  -> fate set by the intracellular compartment alone")
    print(f"  eradication possible           {analysis.decay_condition(p)['can_be_eradicated']}")
    print(f"  floor C_tot (closed form)      {fl['floor_over_k_in']:.4f} x k_in")
    print(f"  floor C_tot (simulated)        {fl['floor_simulated_over_k_in']:.4f} x k_in"
          f"  = {fl['floor_simulated']:.2f} CFU/mL")
    print()

    print("=" * 72)
    print("2. Table 9.4.50  single-lever value needed to bring d to zero")
    print("=" * 72)
    lv = analysis.levers(p)
    need = cf["d"]
    print(f"  psi_in     must FALL to   {max(p['psi_in']-need,0.0):>9.4f} /h"
          f"   (IBC band floor 0.69  -> unreachable)")
    print(f"  k_exf      must RISE to   {p['k_exf']+need:>9.4f} /h"
          f"   (plausible range to 0.1 -> unreachable)")
    print(f"  k_emerge   must RISE to   {p['k_emerge']+need:>9.4f} /h"
          f"   (IBC timeline gives 0.04-0.10 -> unreachable)")
    print(f"  k_A^in     must RISE to   {cf['kill_in']+need:>9.4f} /h"
          f"   (ceiling psi_max-psi_min = "
          f"{p['psi_max']-p['psi_min']:.2f} -> {100*(cf['kill_in']+need)/(p['psi_max']-p['psi_min']):.1f}% of span)")
    eta_req, a_req = lv["eta"], lv["a_peak"]
    print(f"  eta        must RISE to   "
          f"{eta_req if not np.isnan(eta_req) else float('nan'):>9.4f} -"
          f"        (baseline {p['eta']}, plausible to 0.1)")
    print(f"  a_peak     must RISE to   "
          f"{a_req if not np.isnan(a_req) else float('nan'):>9.1f} ug/mL"
          f"   (clinical ceiling 440)")
    print()

    print("=" * 72)
    print("3. Table 9.4.51  the IBC band  (QIR rows removed)")
    print("=" * 72)
    hdr = (f"  {'psi_in':>7}{'d':>10}{'Lambda':>10}{'C_tot':>11}"
           f"{'CFU/bladder':>13}{'dlog10':>9}{'A':>4}{'B':>4}")
    print(hdr); print("  " + "-" * (len(hdr) - 2))
    rows = []
    for psi in (0.69, 1.19, 1.30, 1.39):
        q = dict(p, psi_in=psi)
        c = closed_form(q)
        ev = metrics.evaluate(q)
        rows.append((psi, c, ev))
        print(f"  {psi:>7.2f}{c['d']:>+10.4f}{c['lambda']:>+10.4f}"
              f"{ev['final_total']:>11.2f}{ev['final_total']*VOLUME_ML:>13.4g}"
              f"{ev['depth_log']:>9.2f}"
              f"{'Y' if ev['pass_depth'] else 'N':>4}"
              f"{'Y' if ev['pass_durability'] else 'N':>4}")
    print()

    print("=" * 72)
    print("4. Adaptive ceiling on K^in  (judgement B at 100 CFU/mL)")
    print("=" * 72)
    ratio = fl["floor_simulated_over_k_in"]
    print(f"  C_tot / k_in (simulated)       {ratio:.4f}")
    print(f"  K^in,max = 100 / ratio         {CURE/ratio:.2f} CFU/mL"
          f"  = {CURE/ratio*VOLUME_ML:.4g} CFU/bladder   [brief's arithmetic]")
    f = lambda kin: metrics.evaluate(p, k_in=kin)["final_total"] - CURE
    kin_exact = float(brentq(f, 10.0, 300.0, xtol=0.05))
    print(f"  K^in,max  (exact solve)        {kin_exact:.2f} CFU/mL"
          f"  = {kin_exact*VOLUME_ML:.4g} CFU/bladder")
    print()

    print("=" * 72)
    print("5. Result (6)  bacteria remaining at 72 h")
    print("=" * 72)
    ev = metrics.evaluate(p)
    print(f"  C_tot final cycle              {ev['final_total']:.2f} CFU/mL")
    print(f"  remaining at 72 h              {ev['final_total']*VOLUME_ML:.4g} CFU/bladder")
    print(f"  nadir                          {ev['nadir_total']:.4g} CFU/mL"
          f"  = {ev['nadir_total']*VOLUME_ML:.4g} CFU/bladder"
          f"  at t = {ev['t_nadir']:.1f} h")


if __name__ == "__main__":
    main()
