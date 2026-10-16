"""Endpoints: how a simulated course is turned into a verdict.

Two criteria, because they answer different questions and a regimen can pass
one and fail the other:

  depth       how far down the course drives the population, as a log
              reduction at the nadir. This is the form the pharmacodynamic
              literature uses (stasis / 1-log / 2-log / 3-log).

  durability  whether it stays there, as the bacterial density over the final
              dosing cycle. This is the form regulators use for microbiological
              eradication, though they assess it after therapy ends rather than
              during it.

Kanamycin alone passes depth and fails durability: it reaches a 3.3-log nadir
and then rebounds. Reporting only one of the two hides that.

A caution on the durability endpoint. Whenever the intracellular net growth
rate is positive, the system settles at a floor set by the intracellular
carrying capacity, and the final-cycle total is then ``~0.95 * k_in``
regardless of dose. Check ``analysis.floor_level()`` before reading a
durability number as a statement about the regimen.
"""
from __future__ import annotations

import numpy as np

from . import parameters as prm
from .model import simulate


def geometric_mean(t, c, t_from: float, t_to: float) -> float:
    """Geometric mean of ``c`` over a time window, by trapezoid in log space."""
    m = (t >= t_from) & (t <= t_to)
    tt, cc = t[m], np.maximum(c[m], 1e-30)
    if tt.size < 2:
        return float(cc[-1]) if cc.size else float("nan")
    return float(np.exp(np.trapezoid(np.log(cc), tt) / (tt[-1] - tt[0])))


def evaluate(p: dict | None = None, **overrides) -> dict:
    """Run one course and return every endpoint, plus the verdicts."""
    p = dict(p or prm.values())
    p.update(overrides)

    t, ex, inn, _ = simulate(p)
    tot = ex + inn
    t_end, cycle = p["t_end"], p["tau_a"]

    gm_ex = geometric_mean(t, ex, t_end - cycle, t_end)
    gm_tot = geometric_mean(t, tot, t_end - cycle, t_end)
    i = int(np.argmin(tot))
    depth = float(np.log10(tot[0] / max(tot.min(), 1e-30)))

    return {
        "depth_log": depth,
        "t_nadir": float(t[i]),
        "nadir_total": float(tot.min()),
        "final_ex": gm_ex,
        "final_total": gm_tot,
        "pass_depth": depth >= p["nadir_target"],
        "pass_durability": gm_tot < p["cure_density"],
        "pass_both": depth >= p["nadir_target"] and gm_tot < p["cure_density"],
    }


def minimum_concentration(criterion: str = "durability", p: dict | None = None,
                          lo: float = 20.0, hi: float = 3000.0,
                          tol: float = 0.5, **overrides) -> float:
    """Lowest peak urinary concentration that passes a criterion.

    ``criterion`` is one of:

        depth          log reduction at the nadir
        durability     final-cycle total bacterial density
        both           all of the above
        durability_ex  final-cycle EXTRACELLULAR density only. Kept because
                       earlier revisions judged on this, and it is what the
                       published minimum concentrations were computed with. It
                       omits the intracellular compartment, which is the one
                       that grows.

    Returns NaN if the criterion is not reachable anywhere in ``[lo, hi]``.
    """
    from scipy.optimize import brentq

    p = dict(p or prm.values())
    p.update(overrides)
    if criterion not in ("depth", "durability", "both", "durability_ex"):
        raise ValueError(f"unknown criterion: {criterion}")

    def margin(a: float) -> float:
        """Positive when failing, negative when passing."""
        r = evaluate(p, a_peak=a)
        md = p["nadir_target"] - r["depth_log"]
        mv = np.log10(r["final_total"]) - np.log10(p["cure_density"])
        mx = np.log10(r["final_ex"]) - np.log10(p["cure_density"])
        return {"depth": md, "durability": mv, "both": max(md, mv),
                "durability_ex": mx}[criterion]

    if margin(lo) <= 0:
        return lo
    if margin(hi) > 0:
        return float("nan")
    return float(brentq(margin, lo, hi, xtol=tol))


def compare_regimens(p: dict | None = None, **overrides) -> dict[str, dict]:
    """Untreated, antibiotic alone, SUNRISE alone, and the combination."""
    p = dict(p or prm.values())
    p.update(overrides)
    a, d = p["a_peak"], p["dose"]
    return {
        "untreated": evaluate(p, a_peak=0.0, dose=0.0),
        "antibiotic only": evaluate(p, a_peak=a, dose=0.0),
        "SUNRISE only": evaluate(p, a_peak=0.0, dose=d),
        "combination": evaluate(p, a_peak=a, dose=d),
    }
