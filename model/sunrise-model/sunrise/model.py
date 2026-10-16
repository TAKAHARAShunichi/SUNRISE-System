"""The SUNRISE within-host model.

Five processes on one time axis:

  1. kanamycin concentration in urine, as a repeated-dose profile
  2. CPI particles encountering bacteria and injecting TPD
  3. TPD degrading the resistance factor inside each bacterium
  4. the resistance factor setting that bacterium's MIC
  5. four bacterial compartments growing, dying and exchanging

State vector (7 entries)::

    0  C_S_ex   susceptible, extracellular     CFU/mL
    1  C_R_ex   resistant, extracellular       CFU/mL
    2  C_S_in   susceptible, intracellular     CFU/mL
    3  C_R_in   resistant, intracellular       CFU/mL
    4  C_P      CPI particles                  particles/mL
    5  A_cell   TPD load per bacterium         molecules/bacterium
    6  B_cell   resistance factor              uM

Everything is read from ``parameters.py``; nothing is hard-coded here.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from . import parameters as prm

# Smoothing width for the intracellular capacity gate. Numerical only: it keeps
# max(0, 1 - x) differentiable for the stiff solver and has no biological content.
GATE_EPS = 1.2e-3


# ---------------------------------------------------------------- drug profile
def drug(t: float, p: dict) -> float:
    """Urinary kanamycin concentration at time ``t``.

    Doses of height ``a_peak`` every ``tau_a`` hours, each decaying at
    ``k_elim``, summed. The superposition is a geometric series, so it has a
    closed form and needs no extra state. Dosing stops after ``n_dose_a``.
    """
    if t < 0:
        return 0.0
    n = min(int(t // p["tau_a"]), int(p["n_dose_a"]) - 1)
    s = t - n * p["tau_a"]
    q = np.exp(-p["k_elim"] * p["tau_a"])
    return p["a_peak"] * np.exp(-p["k_elim"] * s) * (1 - q ** (n + 1)) / (1 - q)


# ------------------------------------------------------------- pharmacodynamics
def kill_rate(a: float, mic: float, p: dict) -> float:
    """Kill rate at concentration ``a`` for a strain of the given MIC.

    The four-parameter form in A/MIC: growth falls from ``psi_max`` towards
    ``psi_min``, and this returns the difference, i.e. the kill contribution.
    """
    if a <= 0.0:
        return 0.0
    x = (a / mic) ** p["kappa"]
    return (p["psi_max"] - p["psi_min"]) * x / (x - p["psi_min"] / p["psi_max"])


def mic_of(b_cell: float, p: dict) -> float:
    """MIC of a bacterium carrying ``b_cell`` of the resistance factor.

    Linear between the susceptible MIC at zero resistance factor and the fully
    resistant MIC at the untreated level.
    """
    frac = b_cell / p["b_ref"]
    return max(p["mic_s"] + frac * (p["mic_r0"] - p["mic_s"]), p["mic_s"])


# ----------------------------------------------------------------- helper terms
def _gate(theta: float) -> float:
    """Smoothed free capacity max(0, 1 - theta)."""
    x = (1.0 - theta) / GATE_EPS
    if x > 40.0:
        return 1.0 - theta
    if x < -40.0:
        return 0.0
    return GATE_EPS * np.log1p(np.exp(x))


def production(mu_div: float, p: dict) -> float:
    """Resistance factor production, interpolated between two regimes.

    ``a_interp`` = 0 gives production independent of growth rate; = 1 gives a
    concentration independent of growth rate. Calibrated so that an untreated
    bacterium growing at ``mu_ref`` sits at ``b_ref``.
    """
    f = ((1.0 - p["a_interp"])
         + p["a_interp"] * (mu_div + p["k_xb"]) / (p["mu_ref"] + p["k_xb"]))
    return (p["mu_ref"] + p["k_xb"]) * p["b_ref"] * f


# ------------------------------------------------------------------------- rhs
def rhs(t: float, y, p: dict):
    """Right-hand side of the system. ``p`` is a plain dict of floats."""
    Sx, Rx, Si, Ri, CP, Ac, Bc = (max(v, 0.0) for v in y)

    Cex, Cin = Sx + Rx, Si + Ri
    gate = _gate(Cin / p["k_in"])
    A = drug(t, p)

    mic_r = mic_of(Bc, p)
    mic_r_in = p["mic_s"] + p["rho_in"] * (p["mic_r0"] - p["mic_s"])
    mu_div = p["psi_max"] * (1.0 - Cex / p["k_ex"])
    conj = p["gamma"] * Rx * Sx

    out_ex = p["k_imm"] + p["k_wash"]
    dSx = (mu_div * Sx - kill_rate(A, p["mic_s"], p) * Sx - out_ex * Sx
           - p["k_invade"] * Sx * gate + p["k_emerge"] * Si - conj)
    dRx = (mu_div * Rx - kill_rate(A, mic_r, p) * Rx - out_ex * Rx
           - p["k_invade"] * Rx * gate + p["k_emerge"] * Ri + conj)

    a_in = p["eta"] * A
    dSi = (p["psi_in"] * Si * gate - kill_rate(a_in, p["mic_s"], p) * Si
           - p["k_exf"] * Si + p["k_invade"] * Sx * gate - p["k_emerge"] * Si)
    dRi = (p["psi_in"] * Ri * gate - kill_rate(a_in, mic_r_in, p) * Ri
           - p["k_exf"] * Ri + p["k_invade"] * Rx * gate - p["k_emerge"] * Ri)

    # Particles are lost to washout, inactivation, and being consumed by encounters.
    dCP = -(p["k_wash"] + p["k_loss"] + p["k_enc"] * Cex) * CP

    # TPD accumulates by injection and is diluted by division and by untreated
    # bacteria entering the resistant compartment.
    lam_p = p["k_enc"] * CP
    phi_in = (p["k_emerge"] * Ri + conj) / Rx if Rx > 1e-30 else 0.0
    dAc = (lam_p * p["eta_mech"] * p["phi_p"]
           - (mu_div + p["k_prot"] + phi_in) * Ac)

    k_tpd = p["k_tpd0"] * Ac / p["n0_ref"]
    dBc = (production(mu_div, p) - (mu_div + p["k_xb"] + k_tpd) * Bc
           + phi_in * (p["b_ref"] - Bc))

    return [dSx, dRx, dSi, dRi, dCP, dAc, dBc]


# -------------------------------------------------------------------- simulate
def dose_times(p: dict) -> list[float]:
    """When CPI doses are given. Empty if the dose is zero."""
    if p["dose"] <= 0:
        return []
    d0 = p.get("delta_t", 0.0) % p["tau_p"]
    n = int(np.ceil(p["t_end"] / p["tau_p"])) + 1
    return [d0 + j * p["tau_p"] for j in range(n)
            if 0.0 <= d0 + j * p["tau_p"] < p["t_end"]]


def simulate(p: dict | None = None, seed_scale: float = 1.0, points_per_hour: int = 55,
             rtol: float = 1e-8, atol: float = 1e-10, **overrides):
    """Integrate the course and return ``(t, C_ex, C_in, state)``.

    CPI doses are discrete events, so the course is integrated in segments and
    the particle compartment is topped up between them.

    ``seed_scale`` multiplies the initial densities, for probing the
    low-density limit without changing anything else.
    """
    p = dict(p or prm.values())
    p.update(overrides)

    doses = dose_times(p)
    breaks = sorted(set([0.0] + doses + [p["t_end"]]))

    y = [p["s0"] * seed_scale, p["r0"] * seed_scale, 0.0, 0.0, 0.0, 0.0, p["b_ref"]]
    ts, ex, inn, full = [], [], [], []
    for t0, t1 in zip(breaks[:-1], breaks[1:]):
        if any(abs(t0 - d) < 1e-9 for d in doses):
            y[4] += p["dose"] / p["volume"]
        grid = np.linspace(t0, t1, max(int((t1 - t0) * points_per_hour) + 2, 20))
        sol = solve_ivp(rhs, (t0, t1), y, method="LSODA", t_eval=grid, args=(p,),
                        rtol=rtol, atol=atol * max(seed_scale, 1e-6), max_step=0.5)
        if not sol.success:
            raise RuntimeError(f"integration failed on [{t0}, {t1}]: {sol.message}")
        ts.append(sol.t)
        ex.append(sol.y[0] + sol.y[1])
        inn.append(sol.y[2] + sol.y[3])
        full.append(sol.y)
        y = list(sol.y[:, -1])

    return (np.concatenate(ts), np.concatenate(ex), np.concatenate(inn),
            np.concatenate(full, axis=1))


# ------------------------------------------------- closed forms worth having
def cpi_attachment(p: dict | None = None, **overrides) -> dict:
    """Particle attachment numbers, which follow without integrating anything.

    ``f_hit`` is the fraction of particles that meet a bacterium before being
    lost; ``m0`` the attachments per bacterium from the first dose; ``m_inf``
    the accumulated attachments at the end of a long course.
    """
    p = dict(p or prm.values())
    p.update(overrides)
    c0 = p["s0"] + p["r0"]
    f_hit = p["k_enc"] * c0 / (p["k_enc"] * c0 + p["k_wash"] + p["k_loss"])
    m0 = p["dose"] * f_hit / (c0 * p["volume"])
    ratio = 1.0 + p["k_enc"] * c0 / (p["k_wash"] + p["k_loss"])
    return {"f_hit": f_hit, "m0": m0, "m_inf": m0 * ratio, "ratio": ratio,
            "tpd_load": p["eta_mech"] * p["phi_p"] * m0}
