"""Judgement metrics and derived read-outs (Ver.23, 治療成否の判定基準)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import Solution, adhesion, m_infinity
from .params import Params

FLOOR = 1e-12  # CFU/mL, guard for the geometric mean


@dataclass
class Verdict:
    delta_log: float
    C_bar: float
    nadir: float
    t_nadir: float
    residual_CFU: float  # CFU/bladder at t_end
    passes_depth: bool
    passes_durability: bool

    @property
    def passes(self) -> bool:
        return self.passes_depth and self.passes_durability


def delta_log10(sol: Solution) -> float:
    C = sol.C_tot
    return float(np.log10(C[0] / max(C.min(), FLOOR)))


def c_bar_final_cycle(sol: Solution) -> float:
    """Geometric mean of C^tot over the last antibiotic cycle."""
    p = sol.p
    lo = p.t_end - p.tau_A
    sel = sol.t >= lo
    t, C = sol.t[sel], np.maximum(sol.C_tot[sel], FLOOR)
    integral = np.trapezoid(np.log(C), t)
    return float(np.exp(integral / (t[-1] - t[0])))


def evaluate(sol: Solution) -> Verdict:
    p = sol.p
    C = sol.C_tot
    i = int(np.argmin(C))
    dl = delta_log10(sol)
    cb = c_bar_final_cycle(sol)
    return Verdict(
        delta_log=dl,
        C_bar=cb,
        nadir=float(C[i]),
        t_nadir=float(sol.t[i]),
        residual_CFU=float(C[-1] * p.V),
        passes_depth=dl >= p.delta_log_target,
        passes_durability=cb < p.cure_density,
    )


def window_occupancy(sol: Solution) -> dict:
    """Fraction of time A_cell sits inside the effective window."""
    p = sol.p
    A = sol.A_cell
    inside = (A >= p.window_lo) & (A <= p.window_hi)
    last = sol.t >= p.t_end - p.tau_A
    return {
        "A_cell_max": float(A.max()),
        "A_cell_final_cycle_min": float(A[last].min()),
        "A_cell_final_cycle_max": float(A[last].max()),
        "inside_window_course": float(inside.mean()),
        "inside_window_final_cycle": float(inside[last].mean()),
        "exceeds_window_top": bool(A.max() > p.window_hi),
    }


def adhesion_summary(sol: Solution) -> dict:
    p = sol.p
    return {
        "m0": adhesion(p.C0_ex, p),
        "m_inf_closed_form": m_infinity(p),
        "m_max_simulated": float(np.nanmax(sol.m[np.isfinite(sol.m)])),
        "N0_initial": adhesion(p.C0_ex, p) * p.nu,
        "zero_class_max": float(np.nanmax(sol.zero_class)),
        "zero_class_at_t0": float(sol.zero_class[0]),
    }


def summarise(sol: Solution) -> dict:
    v = evaluate(sol)
    p = sol.p
    out = {
        "MOI": p.MOI,
        "D": p.D,
        "tau_P": p.tau_P,
        "dt_phase": p.dt_phase,
        "n_doses_P": p.n_doses_P,
        "particles_per_course": p.D * p.n_doses_P,
        "dose_g": p.dose_g,
        "daily_g": p.daily_g,
        "serum_peak": p.serum_peak,
        "A_peak_urine": p.A_peak,
        "delta_log": v.delta_log,
        "C_bar": v.C_bar,
        "nadir": v.nadir,
        "t_nadir": v.t_nadir,
        "residual_CFU_bladder": v.residual_CFU,
        "passes_depth": v.passes_depth,
        "passes_durability": v.passes_durability,
        "passes": v.passes,
        "rho_min": float(sol.rho.min()),
        "MIC_min": float(sol.MIC.min()),
    }
    out.update(adhesion_summary(sol))
    out.update(window_occupancy(sol))
    return out
