"""A_cell window diagnostics for the SUNRISE-2 recalculation (task 1).

The effective window on the TPD load per bacterium is 93-1991 molecules/CFU.
Only the UPPER end is treated as a binding constraint: below the window
k_TPD is extrapolated towards zero, which runs towards the correct limit
(no TPD, no degradation), while above it ClpXP processing capacity is
expected to saturate, so a linear extrapolation OVERSTATES the degradation
rate. See the task brief.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "sunrise-model"))

from sunrise import parameters as prm            # noqa: E402
from sunrise.model import simulate               # noqa: E402

WINDOW_LO, WINDOW_HI = 93.0, 1991.0
A_CELL = 5          # index of A_cell in the state vector


def calibrated_k_tpd0(verbose: bool = False) -> float:
    """k_tpd0 from validate.py's calibration, at the published conditions.

    Imported rather than reimplemented so there is exactly one calibration in
    the project. The calibration pins every parameter that has moved since the
    published number was produced, so changing k_in, psi_in, mic_r0 or phi_p
    afterwards does not move it.
    """
    import importlib.util
    path = Path(__file__).resolve().parents[1] / "sunrise-model" / "scripts" / "validate.py"
    spec = importlib.util.spec_from_file_location("_validate", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.calibrate(verbose=verbose)


def base_params(**overrides) -> dict:
    """Ver.22 basis: current defaults with the calibrated degradation rate."""
    p = prm.with_overrides(k_tpd0=calibrated_k_tpd0())
    p.update(overrides)
    return p


def _occupancy(t, ac, lo: float, hi: float) -> float:
    """Time-weighted fraction of the span spent inside [lo, hi]."""
    inside = ((ac >= lo) & (ac <= hi)).astype(float)
    return float(np.trapezoid(inside, t) / (t[-1] - t[0]))


def acell_stats(p: dict, points_per_hour: int = 200) -> dict:
    """A_cell range and window occupancy, over the final cycle and the course."""
    t, _, _, state = simulate(p, points_per_hour=points_per_hour)
    ac = state[A_CELL]
    t_from = p["t_end"] - p["tau_p"]
    m = t >= t_from
    i = int(np.argmax(ac))
    return {
        "max": float(ac.max()),
        "t_max": float(t[i]),
        "final_min": float(ac[m].min()),
        "final_max": float(ac[m].max()),
        "final_swing": float(ac[m].max() / max(ac[m].min(), 1e-30)),
        "occupancy_final": _occupancy(t[m], ac[m], WINDOW_LO, WINDOW_HI),
        "occupancy_course": _occupancy(t, ac, WINDOW_LO, WINDOW_HI),
        "over_window": float(ac.max()) > WINDOW_HI,
    }
