"""ODE system for the SUNRISE-2 bladder model (Ver.23).

State vector
    0  C_S^ex   extracellular susceptible      [CFU/mL]
    1  C_R^ex   extracellular resistant        [CFU/mL]
    2  C_S^in   intracellular susceptible      [CFU/mL]
    3  C_R^in   intracellular resistant        [CFU/mL]
    4  A_cell   TPD per extracellular resistant cell   [molecules/cell]
    5  B_cell   resistance factor per cell             [uM]
    6  C_P      CPI particle concentration             [particles/mL]
    7  G        danger signal                          [-]
    8  N        effector level (rest = 1)              [-]
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .params import Params

TINY = 1e-300
#: CFU/mL well under one cell in a 200 mL bladder
CELL_FLOOR = 1e-6
N_STATE = 9


# --- pharmacodynamics -------------------------------------------------------


def psi(A: float, MIC: float, p: Params) -> float:
    """Hill-type net growth rate (Regoes form).  psi(MIC) = 0 exactly."""
    if A <= 0.0:
        return p.psi_max
    x = (A / MIC) ** p.kappa
    denom = x - p.psi_min / p.psi_max
    return p.psi_max - (p.psi_max - p.psi_min) * x / denom


def kill(A: float, MIC: float, p: Params) -> float:
    """Killing rate k = psi_max - psi."""
    return p.psi_max - psi(A, MIC, p)


def mic_of_rho(rho: float, p: Params) -> float:
    """MIC as a function of the resistance-factor survival ratio."""
    rho = max(rho, 0.0)
    return p.MIC_S + (rho ** p.h_mic) * (p.MIC_R0 - p.MIC_S)


def antibiotic(t: float, p: Params) -> float:
    """Closed-form superposition of q-tau_A dosing with first-order loss."""
    if not p.antibiotic_on or p.A_peak <= 0.0:
        return 0.0
    n = math.floor(t / p.tau_A)
    s = t - n * p.tau_A
    q = math.exp(-p.k_elim * p.tau_A)
    return p.A_peak * math.exp(-p.k_elim * s) * (1.0 - q ** (n + 1)) / (1.0 - q)


# --- resistance-factor production ------------------------------------------


def pi_prod(mu_div: float, p: Params) -> float:
    """Growth-dependent production, interpolated by a (Ver.23)."""
    base = (p.mu_ref + p.k_XB) * p.B_ref
    return base * ((1.0 - p.a_prod) + p.a_prod * (mu_div + p.k_XB) / (p.mu_ref + p.k_XB))


def k_TPD_of(A_cell: float, p: Params) -> float:
    """Degradation rate added by TPD, normalised at A_cell = N0_ref."""
    a = max(A_cell, 0.0)
    if p.cap_k_TPD_at_window:
        a = min(a, p.window_hi)
    return p.k_TPD0 * a / p.N0_ref


# --- immunity ---------------------------------------------------------------


def k_imm_of(C_ex: float, N: float, p: Params) -> float:
    """Removal rate acting on extracellular bacteria.

    Baseline holds N at rest, which reproduces k_imm_rest at the initial load
    (k_phag / (1 + C0/K_cap) = 0.10 /h).  With induced immunity the same
    expression is read as eps_urine * k_PMN * N_PMN(t) / (1 + C_ex/K_cap),
    where N_PMN(t) = N * N_PMN_rest and k_PMN = k_phag / N_PMN_rest.
    """
    sat = 1.0 + max(C_ex, 0.0) / p.K_cap
    if not p.induced_immunity:
        # N held at rest; eps_urine still scales the phagocytic efficiency, so
        # eps_urine is the same axis as N_PMN / N_PMN_rest in this branch.
        return p.eps_urine * p.k_phag / sat
    return p.eps_urine * p.k_PMN * (N * p.N_PMN_rest) / sat


# --- right-hand side --------------------------------------------------------


def rhs(t: float, y: np.ndarray, p: Params) -> np.ndarray:
    CSex, CRex, CSin, CRin, A_cell, B_cell, C_P, G, N = y

    CSex = max(CSex, 0.0)
    CRex = max(CRex, 0.0)
    CSin = max(CSin, 0.0)
    CRin = max(CRin, 0.0)
    C_P = max(C_P, 0.0)

    C_ex = CSex + CRex
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

    # extracellular
    dCSex = (
        mu_div * CSex
        - kill(A, p.MIC_S, p) * CSex
        - k_out * CSex
        - p.k_invade * CSex * free_in
        + p.k_emerge * CSin
        - conj
    )
    dCRex = (
        mu_div * CRex
        - kill(A, MIC_R, p) * CRex
        - k_out * CRex
        - p.k_invade * CRex * free_in
        + p.k_emerge * CRin
        + conj
    )

    # intracellular
    dCSin = (
        p.psi_max_in * CSin * free_in
        - kill(p.eta * A, p.MIC_S, p) * CSin
        - p.k_exf * CSin
        + p.k_invade * CSex * free_in
        - p.k_emerge * CSin
    )
    dCRin = (
        p.psi_max_in * CRin * free_in
        - kill(p.eta * A, MIC_R_in, p) * CRin
        - p.k_exf * CRin
        + p.k_invade * CRex * free_in
        - p.k_emerge * CRin
    )

    # per-cell molecules, carried by the extracellular resistant pool
    # Phi_in is a per-cell flux; below ~1 cell per bladder the compartment
    # means are meaningless and the quotient blows up, so it is switched off.
    Phi_in = ((p.k_emerge * CRin + conj) / CRex) if CRex > CELL_FLOOR else 0.0
    inflow_TPD = p.nu * p.k_enc * C_P if p.sunrise_on else 0.0
    dA = inflow_TPD - (mu_div + p.k_prot) * A_cell - Phi_in * A_cell

    k_TPD = k_TPD_of(A_cell, p)
    dB = (
        pi_prod(mu_div, p)
        - (mu_div + p.k_XB + k_TPD) * B_cell
        + Phi_in * (p.B_ref - B_cell)
    )

    # particles
    dCP = -(p.k_washout + p.k_loss_cpi + p.k_enc * C_ex) * C_P

    # immunity
    L = C_ex + p.sigma_in * C_in
    sense = L / (p.K_G + L)
    dG = p.d_G * (sense - G)
    dN = p.d_N * (1.0 - N) + p.a_N * G * (1.0 - N / p.N_max)

    return np.array([dCSex, dCRex, dCSin, dCRin, dA, dB, dCP, dG, dN])


# --- solution container -----------------------------------------------------


@dataclass
class Solution:
    t: np.ndarray
    y: np.ndarray
    p: Params

    @property
    def C_ex(self) -> np.ndarray:
        return self.y[0] + self.y[1]

    @property
    def C_in(self) -> np.ndarray:
        return self.y[2] + self.y[3]

    @property
    def C_tot(self) -> np.ndarray:
        return self.C_ex + self.C_in

    @property
    def A_cell(self) -> np.ndarray:
        return self.y[4]

    @property
    def rho(self) -> np.ndarray:
        return self.y[5] / self.p.B_ref

    @property
    def MIC(self) -> np.ndarray:
        return np.array([mic_of_rho(r, self.p) for r in self.rho])

    @property
    def C_P(self) -> np.ndarray:
        return self.y[6]

    @property
    def N_eff(self) -> np.ndarray:
        return self.y[8]

    @property
    def A(self) -> np.ndarray:
        return np.array([antibiotic(tt, self.p) for tt in self.t])

    @property
    def m(self) -> np.ndarray:
        """Mean adhered particles per cell at each output time."""
        return np.array([adhesion(c, self.p) for c in self.C_ex])

    @property
    def zero_class(self) -> np.ndarray:
        """Fraction receiving no TPD molecule in one dose: exp(-phi_bar m)."""
        return np.exp(-self.p.phi_bar * self.m)


# --- closed forms shared with the solver ------------------------------------


def f_hit(C_ex: float, p: Params) -> float:
    reach = p.k_enc * C_ex
    return reach / (reach + p.k_washout + p.k_loss_cpi)


def adhesion(C_ex: float, p: Params) -> float:
    """Mean adhered particles per cell, m = D f_hit / (C_ex V)."""
    if C_ex <= 0.0:
        return float("inf")
    return p.D * f_hit(C_ex, p) / (C_ex * p.V)


def m_infinity(p: Params) -> float:
    """Late-course adhesion: m_0 * (1 + k_enc C_0 / (k_wash + k_loss))."""
    return adhesion(p.C0_ex, p) * (
        1.0 + p.k_enc * p.C0_ex / (p.k_washout + p.k_loss_cpi)
    )


# --- integration ------------------------------------------------------------


def dose_times(p: Params) -> list[float]:
    if not p.sunrise_on or p.D <= 0.0:
        return []
    ts, k = [], 0
    while True:
        t = p.dt_phase + k * p.tau_P
        if t >= p.t_end - 1e-12:
            break
        ts.append(t)
        k += 1
    return ts


def initial_state(p: Params) -> np.ndarray:
    y0 = np.zeros(N_STATE)
    y0[0] = p.S0
    y0[1] = p.R0
    y0[2] = 0.0
    y0[3] = 0.0
    y0[4] = 0.0
    y0[5] = p.B_ref
    y0[6] = 0.0
    y0[7] = 0.0
    y0[8] = 1.0
    return y0


def solve_course(p: Params, n_out: int = 2881, rtol: float = 1e-8,
                 atol: float = 1e-10) -> Solution:
    """Integrate the 72 h course, restarting at every CPI dose event."""
    t_grid = np.linspace(0.0, p.t_end, n_out)
    y_out = np.empty((N_STATE, n_out))

    events = sorted(set(dose_times(p)))
    bounds = [0.0] + [e for e in events if 0.0 < e < p.t_end] + [p.t_end]
    bounds = sorted(set(bounds))

    y = initial_state(p)
    if events and abs(events[0]) < 1e-12:
        y[6] += p.D / p.V

    filled = np.zeros(n_out, dtype=bool)
    for lo, hi in zip(bounds[:-1], bounds[1:]):
        if hi - lo < 1e-12:
            continue
        sol = solve_ivp(
            rhs, (lo, hi), y, args=(p,), method="LSODA",
            rtol=rtol, atol=atol, dense_output=True, max_step=(hi - lo),
        )
        if not sol.success:
            raise RuntimeError(f"integration failed on [{lo}, {hi}]: {sol.message}")
        sel = (t_grid >= lo) & (t_grid <= hi) & (~filled)
        if sel.any():
            y_out[:, sel] = sol.sol(t_grid[sel])
            filled |= sel
        y = sol.y[:, -1].copy()
        # dose event at the right endpoint
        if any(abs(hi - e) < 1e-12 for e in events):
            y[6] += p.D / p.V

    if not filled.all():  # pragma: no cover
        y_out[:, ~filled] = y_out[:, filled][:, -1:]

    y_out[:4] = np.maximum(y_out[:4], 0.0)
    y_out[6] = np.maximum(y_out[6], 0.0)
    return Solution(t=t_grid, y=y_out, p=p)
