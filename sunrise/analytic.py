"""Closed forms that need no integration (Ver.23, Results (1))."""

from __future__ import annotations

import math

from .model import adhesion, f_hit, k_imm_of, kill, m_infinity, mic_of_rho
from .params import Params


def k_bar_crit(p: Params) -> float:
    """Mean killing rate above which the extracellular pool turns over."""
    return p.psi_max - p.k_imm_rest - p.k_washout


def adhesion_span(p: Params) -> float:
    return 1.0 + p.k_enc * p.C0_ex / (p.k_washout + p.k_loss_cpi)


def tau_response(p: Params, mu_div: float, k_TPD: float, Phi_in: float = 0.0) -> float:
    """Relaxation time of B_cell, in minutes."""
    return 60.0 / (mu_div + p.k_XB + k_TPD + Phi_in)


def tau_TPD_eff(p: Params, C_ex: float, mu_div: float) -> float:
    """Effective TPD time constant (slower of the two routes), in minutes."""
    intra = mu_div + p.k_prot
    particle = p.k_washout + p.k_loss_cpi + p.k_enc * C_ex
    return 60.0 / min(intra, particle)


def crossover_density(p: Params) -> float:
    """Bacterial density at which the rate-limiting route swaps."""
    return (p.psi_max - p.k_washout - p.k_loss_cpi) / p.k_enc


def rho_ss(p: Params, mu_div: float, k_TPD: float, Phi_in: float = 0.0) -> float:
    """Quasi-steady resistance-factor survival ratio."""
    from .model import pi_prod

    num = pi_prod(mu_div, p) + Phi_in * p.B_ref
    den = (mu_div + p.k_XB + k_TPD + Phi_in) * p.B_ref
    return num / den


def net_growth_rates(p: Params, k_A_ex: float, k_A_in: float) -> dict:
    """Extracellular / intracellular net growth and the dominant eigenvalue."""
    a_ex = p.psi_max - k_A_ex - p.k_imm_rest - p.k_washout - p.k_invade
    d_in = p.psi_max_in - k_A_in - p.k_exf - p.k_emerge
    # two-compartment linear system [[a_ex, k_emerge], [k_invade, d_in]]
    tr = a_ex + d_in
    det = a_ex * d_in - p.k_emerge * p.k_invade
    disc = math.sqrt(max(tr * tr - 4.0 * det, 0.0))
    return {
        "a_ex": a_ex,
        "d_in": d_in,
        "Lambda": 0.5 * (tr + disc),
    }


def theta_sat_requirement(p: Params, theta_sat: float) -> float:
    """Required k_max given the saturated fraction of a dosing cycle."""
    return k_bar_crit(p) / theta_sat


def rho_in_crit(p: Params) -> float:
    """Intracellular survival ratio at which the drug becomes bactericidal."""
    return p.eta * 1.0


def analytic_table(p: Params) -> dict:
    mu0 = p.psi_max * (1.0 - p.C0_ex / p.K_ex)
    k_TPD_ref = p.k_TPD0  # at A_cell = N0_ref
    return {
        "k_enc_mL_per_h": p.k_enc,
        "k_bar_crit_per_h": k_bar_crit(p),
        "adhesion_span": adhesion_span(p),
        "m0": adhesion(p.C0_ex, p),
        "m_inf": m_infinity(p),
        "N0_at_t0": adhesion(p.C0_ex, p) * p.nu,
        "f_hit_at_t0": f_hit(p.C0_ex, p),
        "tau_response_min": tau_response(p, mu0, k_TPD_ref),
        "tau_TPD_eff_initial_min": tau_TPD_eff(p, p.C0_ex, mu0),
        "tau_TPD_eff_late_min": tau_TPD_eff(p, 1e2, p.psi_max),
        "crossover_density_CFU_per_mL": crossover_density(p),
        "MIC_R0": p.MIC_R0,
        "rho_in_crit": rho_in_crit(p),
        "k_max_requirement_at_theta_0.61": theta_sat_requirement(p, 0.61),
        "mu_ref_per_h": p.mu_ref,
        "k_XB_per_h": p.k_XB,
        "k_TPD0_per_h": p.k_TPD0,
        "k_PMN_implied": p.k_PMN,
        "zero_class_at_m0": math.exp(-p.phi_bar * adhesion(p.C0_ex, p)),
    }
