"""Explicit Poisson zero-class variant (Ver.23 task 7).

The mean-field model carries one compartment-mean ``A_cell`` for every
extracellular resistant cell, so a cell that received no TPD molecule in a
given dose is not represented.  This module adds a naive pool that keeps
``MIC_R0`` and is re-mixed at every dose:

* at each dose a fraction ``f0 = exp(-phi_bar * m)`` of the treated pool is
  moved to the naive pool (it received nothing);
* the complementary fraction ``1 - f0`` of the naive pool becomes treated, and
  its resistance factor is folded into the treated mean.

Between doses the naive pool grows, is killed at ``MIC_R0`` and is removed by
immunity and washout exactly like the treated pool.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from .metrics import FLOOR
from .model import (
    N_STATE,
    Solution,
    adhesion,
    antibiotic,
    dose_times,
    initial_state,
    k_imm_of,
    k_TPD_of,
    kill,
    mic_of_rho,
    pi_prod,
)
from .params import Params

TINY = 1e-300
#: CFU/mL well under one cell in a 200 mL bladder
CELL_FLOOR = 1e-6
IDX_NAIVE = N_STATE  # appended state
N_STATE_Z = N_STATE + 1


def rhs_z(t: float, y: np.ndarray, p: Params) -> np.ndarray:
    CSex, CRex, CSin, CRin, A_cell, B_cell, C_P, G, N, CNex = y

    CSex = max(CSex, 0.0)
    CRex = max(CRex, 0.0)
    CNex = max(CNex, 0.0)
    CSin = max(CSin, 0.0)
    CRin = max(CRin, 0.0)
    C_P = max(C_P, 0.0)

    C_ex = CSex + CRex + CNex
    C_in = CSin + CRin
    theta_in = min(C_in / p.K_in, 1.0)
    free_in = 1.0 - theta_in

    mu_div = p.psi_max * (1.0 - C_ex / p.K_ex)
    A = antibiotic(t, p)

    rho = B_cell / p.B_ref
    MIC_R = mic_of_rho(rho, p)
    MIC_R_in = mic_of_rho(p.rho_in, p)

    k_imm = k_imm_of(C_ex, N, p)
    k_out = k_imm + p.k_washout
    conj = p.gamma_conj * CRex * CSex

    dCSex = (
        mu_div * CSex - kill(A, p.MIC_S, p) * CSex - k_out * CSex
        - p.k_invade * CSex * free_in + p.k_emerge * CSin - conj
    )
    dCRex = (
        mu_div * CRex - kill(A, MIC_R, p) * CRex - k_out * CRex
        - p.k_invade * CRex * free_in + p.k_emerge * CRin + conj
    )
    # naive pool: undegraded resistance factor, no TPD
    dCNex = (
        mu_div * CNex - kill(A, p.MIC_R0, p) * CNex - k_out * CNex
        - p.k_invade * CNex * free_in
    )

    dCSin = (
        p.psi_max_in * CSin * free_in - kill(p.eta * A, p.MIC_S, p) * CSin
        - p.k_exf * CSin + p.k_invade * CSex * free_in - p.k_emerge * CSin
    )
    dCRin = (
        p.psi_max_in * CRin * free_in - kill(p.eta * A, MIC_R_in, p) * CRin
        - p.k_exf * CRin + p.k_invade * (CRex + CNex) * free_in
        - p.k_emerge * CRin
    )

    # Phi_in is a per-cell flux; below ~1 cell per bladder the compartment
    # means are meaningless and the quotient blows up, so it is switched off.
    Phi_in = ((p.k_emerge * CRin + conj) / CRex) if CRex > CELL_FLOOR else 0.0
    inflow_TPD = p.nu * p.k_enc * C_P if p.sunrise_on else 0.0
    dA = inflow_TPD - (mu_div + p.k_prot) * A_cell - Phi_in * A_cell

    k_TPD = k_TPD_of(A_cell, p)
    dB = (
        pi_prod(mu_div, p) - (mu_div + p.k_XB + k_TPD) * B_cell
        + Phi_in * (p.B_ref - B_cell)
    )

    dCP = -(p.k_washout + p.k_loss_cpi + p.k_enc * C_ex) * C_P

    L = C_ex + p.sigma_in * C_in
    sense = L / (p.K_G + L)
    dG = p.d_G * (sense - G)
    dN = p.d_N * (1.0 - N) + p.a_N * G * (1.0 - N / p.N_max)

    return np.array([dCSex, dCRex, dCSin, dCRin, dA, dB, dCP, dG, dN, dCNex])


def _apply_dose(y: np.ndarray, p: Params) -> None:
    """Add particles, then re-partition treated / naive by the zero class."""
    C_ex = y[0] + y[1] + y[9]
    m = adhesion(C_ex, p) if C_ex > 0 else 0.0
    f0 = float(np.exp(-p.phi_bar * m)) if np.isfinite(m) else 0.0

    treated, naive = y[1], y[9]
    # treated cells that received nothing fall back to naive
    to_naive = f0 * treated
    # naive cells that received something join the treated mean
    to_treated = (1.0 - f0) * naive

    new_treated = treated - to_naive + to_treated
    new_naive = naive - to_treated + to_naive

    if new_treated > 0.0:
        # joiners arrive with A = 0 and B = B_ref: dilute the treated means
        w_old = max(treated - to_naive, 0.0) / new_treated
        w_new = to_treated / new_treated
        y[4] = y[4] * w_old
        y[5] = y[5] * w_old + p.B_ref * w_new
    else:
        y[4], y[5] = 0.0, p.B_ref

    y[1], y[9] = new_treated, new_naive
    y[6] += p.D / p.V


def solve_course_z(p: Params, n_out: int = 1441, rtol: float = 1e-8,
                   atol: float = 1e-10) -> Solution:
    t_grid = np.linspace(0.0, p.t_end, n_out)
    y_out = np.empty((N_STATE_Z, n_out))

    events = sorted(set(dose_times(p)))
    bounds = sorted(set([0.0] + [e for e in events if 0.0 < e < p.t_end] + [p.t_end]))

    y = np.zeros(N_STATE_Z)
    y[:N_STATE] = initial_state(p)
    # the whole resistant inoculum starts naive
    y[9], y[1] = p.R0, 0.0
    y[5] = p.B_ref
    if events and abs(events[0]) < 1e-12:
        _apply_dose(y, p)

    filled = np.zeros(n_out, dtype=bool)
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        if hi - lo < 1e-12:
            continue
        sol = solve_ivp(rhs_z, (lo, hi), y, args=(p,), method="LSODA",
                        rtol=rtol, atol=atol, dense_output=True,
                        max_step=(hi - lo))
        if not sol.success:
            raise RuntimeError(f"zero-class integration failed: {sol.message}")
        sel = (t_grid >= lo) & (t_grid <= hi) & (~filled)
        if sel.any():
            y_out[:, sel] = sol.sol(t_grid[sel])
            filled |= sel
        y = sol.y[:, -1].copy()
        if any(abs(hi - e) < 1e-12 for e in events):
            _apply_dose(y, p)

    if not filled.all():  # pragma: no cover
        y_out[:, ~filled] = y_out[:, filled][:, -1:]

    y_out[:4] = np.maximum(y_out[:4], 0.0)
    y_out[9] = np.maximum(y_out[9], 0.0)
    y_out[6] = np.maximum(y_out[6], 0.0)

    # fold the naive pool into the extracellular resistant slot for the
    # standard Solution read-outs, keeping the split available separately
    merged = y_out[:N_STATE].copy()
    merged[1] = y_out[1] + y_out[9]
    sol_obj = Solution(t=t_grid, y=merged, p=p)
    sol_obj.naive = y_out[9]  # type: ignore[attr-defined]
    sol_obj.treated = y_out[1]  # type: ignore[attr-defined]
    return sol_obj
