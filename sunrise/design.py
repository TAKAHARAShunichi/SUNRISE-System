"""Deterministic design search for the standard patient (Ver.23 Results 7).

The probabilistic chance-constraint search is retired; what remains is
bisection and small deterministic grids under the safety box constraint.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Iterable

import numpy as np

from .metrics import evaluate, summarise
from .model import solve_course
from .params import DOSE_MAX_G, Params

N_OUT = 1441


def passes(p: Params, n_out: int = N_OUT) -> bool:
    return evaluate(solve_course(p, n_out=n_out)).passes


def a_cell_max(p: Params, n_out: int = N_OUT) -> float:
    return float(solve_course(p, n_out=n_out).A_cell.max())


def bisect_min_log(
    predicate: Callable[[float], bool],
    lo: float,
    hi: float,
    tol: float = 1e-3,
    max_iter: int = 60,
) -> float | None:
    """Smallest x in [lo, hi] (log-spaced) with predicate(x) true.

    Returns None when predicate(hi) is false.  The returned endpoint always
    satisfies the predicate: the bracket is closed on the passing side.
    """
    if not predicate(hi):
        return None
    if predicate(lo):
        return lo
    llo, lhi = math.log10(lo), math.log10(hi)
    for _ in range(max_iter):
        if lhi - llo < tol:
            break
        mid = 0.5 * (llo + lhi)
        if predicate(10.0**mid):
            lhi = mid
        else:
            llo = mid
    return 10.0**lhi


def bisect_max_log(
    predicate: Callable[[float], bool],
    lo: float,
    hi: float,
    tol: float = 1e-3,
    max_iter: int = 60,
) -> float | None:
    """Largest x in [lo, hi] (log-spaced) with predicate(x) true."""
    if not predicate(lo):
        return None
    if predicate(hi):
        return hi
    llo, lhi = math.log10(lo), math.log10(hi)
    for _ in range(max_iter):
        if lhi - llo < tol:
            break
        mid = 0.5 * (llo + lhi)
        if predicate(10.0**mid):
            llo = mid
        else:
            lhi = mid
    return 10.0**llo


# --- task 1 helpers ---------------------------------------------------------


def min_moi_at_dose(
    base: Params, dose_g: float, lo: float = 1e-2, hi: float = 1e4
) -> float | None:
    """Smallest MOI passing both judgements at a fixed per-dose mass."""
    p0 = base.with_dose_g(dose_g)
    return bisect_min_log(lambda moi: passes(p0.with_MOI(moi)), lo, hi)


def moi_cap(base: Params, lo: float = 1e-2, hi: float = 1e4) -> float | None:
    """Largest MOI keeping A_cell inside the effective window for the course."""
    top = base.window_hi
    return bisect_max_log(lambda moi: a_cell_max(base.with_MOI(moi)) <= top, lo, hi)


def min_moi_unconstrained_dose(
    base: Params, dose_hi: float = 6.0, lo: float = 1e-2, hi: float = 1e4
) -> float | None:
    """Smallest MOI passing at *any* antibiotic dose (safety ignored).

    Judgement improves monotonically with dose, so the generous ceiling
    ``dose_hi`` stands in for 'whatever dose it takes'.
    """
    return min_moi_at_dose(base, dose_hi, lo=lo, hi=hi)


@dataclass
class Task1Result:
    tau_P: float
    moi_min_any_dose: float | None
    moi_cap: float | None
    moi_optimal_safe: float | None
    dose_used_g: float
    compatible: bool = field(init=False)

    def __post_init__(self) -> None:
        self.compatible = (
            self.moi_optimal_safe is not None
            and self.moi_cap is not None
            and self.moi_optimal_safe <= self.moi_cap
        )


def task1(base: Params, tau_P: float, dose_cap_g: float = DOSE_MAX_G) -> Task1Result:
    p = base.with_(tau_P=tau_P)
    return Task1Result(
        tau_P=tau_P,
        moi_min_any_dose=min_moi_unconstrained_dose(p),
        moi_cap=moi_cap(p),
        moi_optimal_safe=min_moi_at_dose(p, dose_cap_g),
        dose_used_g=dose_cap_g,
    )


def tighten_tau_P(
    base: Params, candidates: Iterable[float] = (6.0, 4.0, 3.0, 2.0)
) -> tuple[list[Task1Result], float | None]:
    """Walk tau_P down until the optimal MOI fits under MOI_cap.

    Returns every evaluated point and the largest compatible tau_P.
    """
    results = [task1(base, t) for t in candidates]
    ok = [r.tau_P for r in results if r.compatible]
    return results, (max(ok) if ok else None)


# --- task 2: two-dimensional deterministic search --------------------------


def objective(p: Params, m_inf_ref: float, D_tot_ref: float,
              w_A: float = 1.0, w_m: float = 0.5, w_D: float = 0.25) -> float:
    """J from Ver.23 Table 9.4.16, normalised by the safety cap and a reference."""
    from .model import m_infinity

    return (
        w_A * p.daily_g / 1.43
        + w_m * m_infinity(p) / m_inf_ref
        + w_D * (p.D * p.n_doses_P) / D_tot_ref
    )


def best_phase(
    p: Params, n_phase: int = 8
) -> tuple[float, float]:
    """Phase offset maximising delta_log, and that delta_log."""
    best = (0.0, -math.inf)
    for dt in np.linspace(0.0, p.tau_P, n_phase, endpoint=False):
        v = evaluate(solve_course(p.with_(dt_phase=float(dt))))
        score = v.delta_log + (10.0 if v.passes else 0.0)
        if score > best[1]:
            best = (float(dt), score)
    dt = best[0]
    return dt, evaluate(solve_course(p.with_(dt_phase=dt))).delta_log
