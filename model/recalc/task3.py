#!/usr/bin/env python3
"""Tables 9.4.30 and 9.4.31 on the new baseline.

These report the EFFICACY FLOOR D_min -- the smallest dose that passes both
judgements -- which is a different quantity from task 1's D_max, the largest
dose that keeps A_cell inside the upper window bound. The usable design range
is [D_min, D_max].
"""
from __future__ import annotations

import sys
import time

import numpy as np
from scipy.optimize import brentq

import acell
from sunrise import metrics
from sunrise.model import cpi_attachment

LO, HI = 1.0e9, 2.0e13
NADIR_TARGET, CURE = 3.0, 100.0


def margin(D: float, p: dict) -> float:
    """Positive when the course fails either judgement."""
    r = metrics.evaluate(p, dose=D)
    return max(NADIR_TARGET - r["depth_log"],
               np.log10(r["final_total"]) - np.log10(CURE))


def d_min(p: dict) -> float:
    """Smallest dose passing both judgements, or NaN if unreachable."""
    if margin(LO, p) <= 0:
        return LO
    if margin(HI, p) > 0:
        return float("nan")
    return float(brentq(lambda D: margin(D, p), LO, HI, rtol=1e-4))


def m_inf_of(p: dict, D: float) -> float:
    return cpi_attachment(dict(p, dose=D))["m_inf"]


def main() -> None:
    p = acell.base_params()
    t0 = time.time()

    print("=" * 74)
    print("Table 9.4.30  minimum dose through the worst dosing phase, by tau_P")
    print("=" * 74)
    print(f"  {'tau_P':>6}{'n_dose':>8}{'worst dt':>10}{'D_min':>12}"
          f"{'per course':>12}{'m_inf':>9}{'D_max':>12}")
    print("  " + "-" * 69)
    for tau in (2.0, 3.0, 4.0, 6.0, 8.0, 12.0):
        n_phase = 6
        phases = [k * tau / n_phase for k in range(n_phase)]
        worst_D, worst_dt = -1.0, None
        for dt in phases:
            q = dict(p, tau_p=tau, delta_t=dt)
            D = d_min(q)
            if np.isnan(D):
                worst_D, worst_dt = float("nan"), dt
                break
            if D > worst_D:
                worst_D, worst_dt = D, dt
        q = dict(p, tau_p=tau, delta_t=worst_dt)
        n_dose = len(__import__("sunrise.model", fromlist=["dose_times"]).dose_times(q))
        # window ceiling at the same tau_P / worst phase
        try:
            Dmax = float(brentq(
                lambda D: acell.acell_stats(dict(q, dose=D))["max"] - acell.WINDOW_HI,
                1.0e10, 2.0e12, xtol=1e8, rtol=1e-5))
        except ValueError:
            Dmax = float("nan")
        print(f"  {tau:>6.0f}{n_dose:>8d}{worst_dt:>10.2f}{worst_D:>12.4g}"
              f"{worst_D*n_dose:>12.4g}{m_inf_of(q, worst_D):>9.1f}{Dmax:>12.4g}")
    print(f"\n  ({time.time()-t0:.0f} s)\n")

    print("=" * 74)
    print("Table 9.4.31  minimum dose and attachment against phi_p")
    print("=" * 74)
    print(f"  {'phi_p':>7}{'D_min':>12}{'m_inf at D_min':>16}"
          f"{'D_max':>12}{'m_inf at D_max':>16}")
    print("  " + "-" * 69)
    base = None
    for pp in (10.7, 20.3, 34.9):
        q = dict(p, phi_p=pp)
        Dn = d_min(q)
        Dx = float(brentq(
            lambda D: acell.acell_stats(dict(q, dose=D))["max"] - acell.WINDOW_HI,
            1.0e10, 2.0e12, xtol=1e8, rtol=1e-5))
        print(f"  {pp:>7.1f}{Dn:>12.4g}{m_inf_of(q, Dn):>16.1f}"
              f"{Dx:>12.4g}{m_inf_of(q, Dx):>16.1f}")
        if base is None:
            base = (pp, Dn, Dx)
    bp, bn, bx = base
    print(f"\n  1/phi_p scaling check against phi_p = {bp}:")
    for pp in (20.3, 34.9):
        q = dict(p, phi_p=pp)
        Dn, Dx = d_min(q), float(brentq(
            lambda D: acell.acell_stats(dict(q, dose=D))["max"] - acell.WINDOW_HI,
            1.0e10, 2.0e12, xtol=1e8, rtol=1e-5))
        print(f"    phi_p {pp:>5.1f}:  D_min x{Dn/bn:.4f} (1/ratio = {bp/pp:.4f})"
              f"   D_max x{Dx/bx:.4f} (1/ratio = {bp/pp:.4f})")
    print(f"\n  (total {time.time()-t0:.0f} s)")


if __name__ == "__main__":
    main()
