"""Structural analysis: what the model can and cannot do, before simulating.

Three results that a dose sweep will not tell you, and that cost almost
nothing to compute.

**The decay condition.** In the low-density limit the resistant population is
a linear two-compartment system, so eradication is an eigenvalue question::

    a_bar = psi_max - kill(A) - k_imm - k_wash - k_invade     (extracellular)
    d     = psi_in  - kill(eta A) - k_exf - k_emerge          (intracellular)

Both eigenvalues are negative only if ``a_bar + d < 0`` and
``a_bar * d > k_emerge * k_invade``. Since ``a_bar < 0`` under any useful
regimen, the second condition forces ``d < 0``. **If the intracellular
compartment grows on its own, no amount of extracellular killing can drive the
system to extinction** -- the determinant is negative and one eigenvalue stays
positive whatever the dose.

**The floor.** When ``d > 0`` the system does not grow without bound either:
the intracellular compartment fills to its carrying capacity and stays there,
and the total settles at ``k_in * (1 + k_emerge / |a_bar|)``. The dose enters
only through the correction term. A durability endpoint evaluated in this
regime is reporting ``k_in``, not the regimen.

**Which levers can flip the sign.** Only terms that appear in ``d``. Immune
clearance, washout and the invasion rate act on the extracellular compartment
alone and cannot flip it at any value. The antibiotic concentration is the
exception among the obvious levers: it reaches the intracellular compartment
through ``eta``, so it does enter ``d``, but at the baseline access fraction it
takes about 2,100 ug/mL -- roughly five times the clinical ceiling -- to bring
``d`` to zero.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from . import parameters as prm
from .model import drug, kill_rate, mic_of, simulate


def _cycle_average(f, t_from: float, t_to: float, n: int = 400) -> float:
    ts = np.linspace(t_from, t_to, n)
    return float(np.trapezoid([f(x) for x in ts], ts) / (t_to - t_from))


def growth_rates(p: dict | None = None, **overrides) -> dict:
    """``a_bar`` and ``d``, averaged over the final dosing cycle.

    The resistance factor level is taken from a simulation, since it depends on
    the whole TPD history.
    """
    p = dict(p or prm.values())
    p.update(overrides)

    t, _, _, state = simulate(p)
    t_from, t_to = p["t_end"] - p["tau_a"], p["t_end"]
    m = (t >= t_from) & (t <= t_to)
    tt, bb = t[m], state[6][m]

    kills = [kill_rate(drug(x, p), mic_of(b, p), p) for x, b in zip(tt, bb)]
    kill_ex = float(np.trapezoid(kills, tt) / (tt[-1] - tt[0]))

    mic_in = p["mic_s"] + p["rho_in"] * (p["mic_r0"] - p["mic_s"])
    kill_in = _cycle_average(
        lambda x: kill_rate(p["eta"] * drug(x, p), mic_in, p), t_from, t_to)

    a_bar = p["psi_max"] - kill_ex - p["k_imm"] - p["k_wash"] - p["k_invade"]
    d = p["psi_in"] - kill_in - p["k_exf"] - p["k_emerge"]
    return {"a_bar": a_bar, "d": d, "kill_ex": kill_ex, "kill_in": kill_in}


def decay_condition(p: dict | None = None, **overrides) -> dict:
    """Whether the population can be driven to extinction, and why not if not."""
    p = dict(p or prm.values())
    p.update(overrides)
    g = growth_rates(p)
    a_bar, d = g["a_bar"], g["d"]

    trace = a_bar + d
    det = a_bar * d - p["k_emerge"] * p["k_invade"]
    disc = trace * trace - 4 * det
    lam = 0.5 * (trace + (np.sqrt(disc) if disc > 0 else 0.0))

    return {**g, "trace": trace, "determinant": det, "lambda_max": lam,
            "can_be_eradicated": bool(trace < 0 and det > 0),
            "blocked_by_intracellular": bool(d >= 0)}


def floor_level(p: dict | None = None, measure: bool = False, **overrides) -> dict:
    """Where the system settles when the intracellular compartment grows.

    The closed form assumes the intracellular compartment has filled and the
    extracellular one is in quasi-equilibrium with what leaks out of it. It
    overestimates the simulated final-cycle geometric mean by roughly 15%,
    because the real trajectory oscillates within the dosing cycle and the
    capacity gate never closes completely. Use it for the scaling, not the
    digits; pass ``measure=True`` for the simulated value.

    The scaling is the point: the floor is proportional to ``k_in`` and the
    dose enters only through ``a_bar`` in the correction term. Across three
    orders of magnitude in ``k_in``, the simulated final-cycle total stays at
    0.95 ``k_in`` with the published parameter set.

    Returns NaN for the floor when ``d < 0``, since there is none.
    """
    p = dict(p or prm.values())
    p.update(overrides)
    g = growth_rates(p)
    if g["d"] < 0:
        out = {**g, "floor_total": float("nan"), "floor_over_k_in": float("nan")}
    else:
        ratio = 1.0 + p["k_emerge"] / abs(g["a_bar"])
        out = {**g, "floor_total": p["k_in"] * ratio, "floor_over_k_in": ratio}
    if measure:
        from .metrics import evaluate
        r = evaluate(p)
        out["floor_simulated"] = r["final_total"]
        out["floor_simulated_over_k_in"] = r["final_total"] / p["k_in"]
    return out


def lambda_numeric(p: dict | None = None, seed_rel: float = 1e-6,
                   n_cycles: int = 4, **overrides) -> float:
    """The asymptotic growth rate, measured rather than linearised.

    Integrates from a density far below both carrying capacities and fits the
    slope over the last ``n_cycles``. The seed is scaled to ``k_in`` so the
    low-density limit genuinely holds -- a fixed seed makes the result look
    like it depends on ``k_in`` when it does not.
    """
    p = dict(p or prm.values())
    p.update(overrides)
    scale = seed_rel * p["k_in"] / p["r0"]
    t, ex, inn, _ = simulate(p, seed_scale=scale, rtol=1e-9)
    c = ex + inn
    m = t >= p["t_end"] - n_cycles * p["tau_a"]
    return float(np.polyfit(t[m], np.log(np.maximum(c[m], 1e-300)), 1)[0])


def levers(p: dict | None = None, **overrides) -> dict[str, float]:
    """For each lever, the value that would bring ``d`` to zero.

    NaN means that lever cannot do it at any value, because it does not appear
    in ``d``. ``a_peak`` does appear, through ``eta``, so it gets a number --
    usually an impractical one.
    """
    p = dict(p or prm.values())
    p.update(overrides)
    g = growth_rates(p)
    need = g["d"]  # how far d must come down

    out = {
        "psi_in": max(p["psi_in"] - need, 0.0),
        "k_exf": p["k_exf"] + need,
        "k_emerge": p["k_emerge"] + need,
        "kill_in": g["kill_in"] + need,
        "eta": required_eta(p),
        "a_peak": required_a_peak(p),
    }
    # These act on the extracellular compartment alone.
    for name in ("k_imm", "k_wash", "k_invade"):
        out[name] = float("nan")
    return out


def required_a_peak(p: dict | None = None, hi: float = 1e5, **overrides) -> float:
    """The peak concentration that would bring ``d`` to zero through ``eta``."""
    p = dict(p or prm.values())
    p.update(overrides)
    try:
        return float(brentq(lambda a: growth_rates(p, a_peak=a)["d"],
                            p["a_peak"], hi, xtol=1.0))
    except ValueError:
        return float("nan")


def required_eta(p: dict | None = None, **overrides) -> float:
    """The intracellular drug access that would make ``d`` zero, if reachable."""
    p = dict(p or prm.values())
    p.update(overrides)
    g = growth_rates(p)
    target = g["kill_in"] + g["d"]
    mic_in = p["mic_s"] + p["rho_in"] * (p["mic_r0"] - p["mic_s"])
    t_from, t_to = p["t_end"] - p["tau_a"], p["t_end"]

    def f(eta: float) -> float:
        avg = _cycle_average(
            lambda x: kill_rate(eta * drug(x, p), mic_in, p), t_from, t_to)
        return avg - target

    try:
        return float(brentq(f, 1e-5, 1.0))
    except ValueError:
        return float("nan")
