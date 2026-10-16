#!/usr/bin/env python3
"""The usable design range [D_min, D_max] at each tau_P.

D_min is the efficacy floor (Table 9.4.30, worst dosing phase); D_max is the
upper-window ceiling from task 1. If D_min > D_max the interval has no
feasible dose at all.
"""
from __future__ import annotations

from scipy.optimize import brentq

import acell
from sunrise.model import cpi_attachment

# (tau_P, worst phase, D_min) from Table 9.4.30
ROWS = [(2.0, 0.00, 4.052e10), (3.0, 1.50, 6.095e10), (4.0, 0.00, 8.296e10),
        (6.0, 3.00, 1.432e11), (8.0, 4.00, 2.091e11), (12.0, 0.00, 1.747e12)]


def d_max(p: dict) -> float:
    try:
        return float(brentq(
            lambda D: acell.acell_stats(dict(p, dose=D))["max"] - acell.WINDOW_HI,
            1.0e9, 1.0e13, xtol=1e8, rtol=1e-5))
    except ValueError:
        return float("nan")


def main() -> None:
    p = acell.base_params()
    print(f"  {'tau_P':>6}{'D_min':>12}{'D_max':>12}{'ratio':>8}"
          f"{'m_inf@Dmin':>12}{'feasible':>10}", flush=True)
    print("  " + "-" * 60)
    for tau, dt, dmin in ROWS:
        q = dict(p, tau_p=tau, delta_t=dt)
        dmax = d_max(q)
        m = cpi_attachment(dict(q, dose=dmin))["m_inf"]
        print(f"  {tau:>6.0f}{dmin:>12.4g}{dmax:>12.4g}{dmax/dmin:>8.2f}"
              f"{m:>12.1f}{'yes' if dmax > dmin else 'NO':>10}", flush=True)


if __name__ == "__main__":
    main()
