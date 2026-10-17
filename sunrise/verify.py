"""Invariant checks (Ver.23, 実装と基準条件).

Every item listed in the manuscript is one function here.  Run as
``python -m sunrise.verify``; a non-zero exit status means an invariant broke.
"""

from __future__ import annotations

import math
import sys

import numpy as np

from . import analytic
from .metrics import c_bar_final_cycle, delta_log10, evaluate
from .model import (
    adhesion,
    initial_state,
    k_TPD_of,
    k_imm_of,
    kill,
    mic_of_rho,
    m_infinity,
    pi_prod,
    psi,
    rhs,
    solve_course,
)
from .params import STANDARD_PATIENT, Params, dose_to_serum, urine_to_dose

CHECKS: list[tuple[str, callable]] = []


def check(name: str):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn

    return deco


def close(a: float, b: float, rtol: float = 1e-6, atol: float = 0.0) -> bool:
    return abs(a - b) <= atol + rtol * abs(b)


# --- pharmacodynamics -------------------------------------------------------


@check("psi(A=MIC) = 0")
def _c1() -> bool:
    p = STANDARD_PATIENT
    return all(close(psi(mic, mic, p), 0.0, atol=1e-12)
               for mic in (p.MIC_S, p.MIC_R0, 48.2))


@check("psi(A=0) = psi_max")
def _c2() -> bool:
    p = STANDARD_PATIENT
    return close(psi(0.0, p.MIC_S, p), p.psi_max)


@check("psi(A->inf) -> psi_min")
def _c3() -> bool:
    p = STANDARD_PATIENT
    return close(psi(1e12, p.MIC_S, p), p.psi_min, rtol=1e-9)


@check("MIC exactly linear in B_cell at h = 1")
def _c4() -> bool:
    p = STANDARD_PATIENT
    lo, hi = mic_of_rho(0.0, p), mic_of_rho(1.0, p)
    mid = mic_of_rho(0.5, p)
    return (close(lo, p.MIC_S) and close(hi, p.MIC_R0)
            and close(mid, 0.5 * (lo + hi)))


@check("k_max = psi_max - psi_min = 1.70 /h")
def _c5() -> bool:
    p = STANDARD_PATIENT
    return close(p.psi_max - p.psi_min, 1.70, rtol=1e-9)


# --- resistance-factor production ------------------------------------------


@check("rho = 1 at mu_ref with no TPD, for any a")
def _c6() -> bool:
    for a in (0.0, 0.25, 0.5, 0.75, 1.0):
        p = STANDARD_PATIENT.with_(a_prod=a)
        r = analytic.rho_ss(p, p.mu_ref, 0.0)
        if not close(r, 1.0, rtol=1e-12):
            return False
    return True


@check("rho = 1 at any mu only when a = 1")
def _c7() -> bool:
    p1 = STANDARD_PATIENT.with_(a_prod=1.0)
    ok_a1 = all(close(analytic.rho_ss(p1, mu, 0.0), 1.0, rtol=1e-12)
                for mu in (0.1, 0.5, 1.16))
    p0 = STANDARD_PATIENT.with_(a_prod=0.0)
    differs = not close(analytic.rho_ss(p0, 0.1, 0.0), 1.0, rtol=1e-6)
    return ok_a1 and differs


@check("k_TPD normalisation: k_TPD(N0_ref) = k_TPD0")
def _c8() -> bool:
    p = STANDARD_PATIENT
    return close(k_TPD_of(p.N0_ref, p), p.k_TPD0)


# --- CPI kinetics -----------------------------------------------------------


@check("adhesion span matches the closed form")
def _c9() -> bool:
    p = STANDARD_PATIENT
    return close(m_infinity(p) / adhesion(p.C0_ex, p), analytic.adhesion_span(p))


@check("N0 = nu * m, phi_bar not applied twice")
def _c10() -> bool:
    p = STANDARD_PATIENT
    m0 = adhesion(p.C0_ex, p)
    return close(analytic.analytic_table(p)["N0_at_t0"], p.nu * m0)


@check("TPD time constant is set by the slower route")
def _c11() -> bool:
    p = STANDARD_PATIENT
    mu0 = p.psi_max * (1.0 - p.C0_ex / p.K_ex)
    early = analytic.tau_TPD_eff(p, p.C0_ex, mu0)
    late = analytic.tau_TPD_eff(p, 1e2, p.psi_max)
    intra_early = 60.0 / (mu0 + p.k_prot)
    particle_late = 60.0 / (p.k_washout + p.k_loss_cpi + p.k_enc * 1e2)
    return close(early, intra_early) and close(late, particle_late)


@check("k_washout from residual volume reproduces Table 9.4.3")
def _c12() -> bool:
    from .params import k_washout_from_residual

    return (close(k_washout_from_residual(125.0), 0.30, rtol=0.02)
            and close(k_washout_from_residual(200.0), 0.24, rtol=0.03)
            and close(k_washout_from_residual(10.0), 0.70, rtol=0.03))


# --- immunity ---------------------------------------------------------------


@check("rest-state immunity reproduces 0.10 /h at the initial load")
def _c13() -> bool:
    p = STANDARD_PATIENT
    return close(k_imm_of(p.C0_ex, 1.0, p), p.k_imm_rest, rtol=1e-9)


@check("induced form equals the fixed form at rest")
def _c14() -> bool:
    p = STANDARD_PATIENT
    q = p.with_(induced_immunity=True, eps_urine=1.0)
    return close(k_imm_of(p.C0_ex, 1.0, q), k_imm_of(p.C0_ex, 1.0, p))


@check("induced N equilibrium at G=1 is 7.75")
def _c15() -> bool:
    p = STANDARD_PATIENT
    eq = (p.d_N + p.a_N) / (p.d_N + p.a_N / p.N_max)
    return close(eq, 7.75, rtol=2e-3)


# --- antibiotic -------------------------------------------------------------


@check("dose <-> urine <-> serum conversions are consistent")
def _c16() -> bool:
    return (close(urine_to_dose(300.0), 0.45)
            and close(dose_to_serum(0.45), 22.0)
            and close(dose_to_serum(0.7159090909), 35.0, rtol=1e-6))


@check("antibiotic closed form matches the iterated sum")
def _c17() -> bool:
    from .model import antibiotic

    p = STANDARD_PATIENT
    for t in (0.0, 3.7, 11.9, 12.0, 35.5, 71.9):
        n = math.floor(t / p.tau_A)
        s = t - n * p.tau_A
        ref = sum(p.A_peak * math.exp(-p.k_elim * (s + j * p.tau_A))
                  for j in range(n + 1))
        if not close(antibiotic(t, p), ref, rtol=1e-13):
            return False
    return True


# --- integration ------------------------------------------------------------


@check("states stay non-negative and bounded")
def _c18() -> bool:
    sol = solve_course(STANDARD_PATIENT.with_MOI(8.0), n_out=1441)
    p = sol.p
    return (bool((sol.y[:4] >= -1e-9).all())
            and bool((sol.C_ex <= p.K_ex * 1.001).all())
            and bool((sol.C_in <= p.K_in * 1.001).all())
            and bool((sol.y[5] >= -1e-12).all())
            and bool((sol.y[8] >= 1.0 - 1e-9).all()))


@check("judgement is independent of the output grid")
def _c19() -> bool:
    p = STANDARD_PATIENT.with_MOI(8.0)
    a = evaluate(solve_course(p, n_out=1441))
    b = evaluate(solve_course(p, n_out=5761))
    return (close(a.delta_log, b.delta_log, rtol=2e-3)
            and close(a.C_bar, b.C_bar, rtol=5e-3))


@check("C_bar scales with K_in (passing design, floor pinned at K_in)")
def _c20() -> bool:
    base = STANDARD_PATIENT.with_MOI(30.0)
    ratios = []
    for K in (15.0, 50.0, 150.0, 500.0):
        v = evaluate(solve_course(base.with_(K_in=K), n_out=1441))
        ratios.append(v.C_bar / K)
    return max(ratios) / min(ratios) < 1.08


@check("deeper judgement is monotone in MOI")
def _c21() -> bool:
    base = STANDARD_PATIENT
    vals = [evaluate(solve_course(base.with_MOI(m), n_out=1441)).delta_log
            for m in (1.0, 3.0, 10.0, 30.0)]
    return all(b >= a - 1e-6 for a, b in zip(vals, vals[1:]))


@check("deeper judgement is monotone in antibiotic dose")
def _c22() -> bool:
    base = STANDARD_PATIENT.with_MOI(8.0)
    vals = [evaluate(solve_course(base.with_dose_g(d), n_out=1441)).delta_log
            for d in (0.1, 0.3, 0.5, 0.72)]
    return all(b >= a - 1e-6 for a, b in zip(vals, vals[1:]))


@check("judgement is monotone in N0 (nu), not scale-invariant")
def _c23() -> bool:
    base = STANDARD_PATIENT.with_MOI(8.0)
    vals = [evaluate(solve_course(base.with_(nu=n), n_out=1441)).delta_log
            for n in (10.7, 20.3, 32.3)]
    return all(b >= a - 1e-6 for a, b in zip(vals, vals[1:])) and vals[-1] > vals[0] + 1e-3


@check("no SUNRISE leaves rho at 1 (a = 1, growth-independent of mu)")
def _c24() -> bool:
    # Ver.23: rho == 1 holds for any a only at mu_div = mu_ref; over a course
    # mu_div moves, so the course-wide invariant is the a = 1 case.
    sol = solve_course(
        STANDARD_PATIENT.with_(sunrise_on=False, a_prod=1.0), n_out=721
    )
    at_ref = analytic.rho_ss(STANDARD_PATIENT, STANDARD_PATIENT.mu_ref, 0.0)
    return bool(np.all(np.abs(sol.rho - 1.0) < 1e-6)) and close(at_ref, 1.0)


@check("quasi-steady rho matches the integrated rho when TPD is flat")
def _c25() -> bool:
    p = STANDARD_PATIENT.with_MOI(8.0)
    sol = solve_course(p, n_out=2881)
    i = int(np.argmin(np.abs(sol.t - 65.0)))
    mu = p.psi_max * (1.0 - sol.C_ex[i] / p.K_ex)
    Phi = p.k_emerge * sol.y[3][i] / max(sol.y[1][i], 1e-300)
    pred = analytic.rho_ss(p, mu, k_TPD_of(sol.A_cell[i], p), Phi)
    return close(sol.rho[i], pred, rtol=0.25)


def main() -> int:
    failed = []
    for name, fn in CHECKS:
        try:
            ok = bool(fn())
        except Exception as exc:  # pragma: no cover
            ok, name = False, f"{name}  [raised {exc!r}]"
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
        if not ok:
            failed.append(name)
    print(f"\n{len(CHECKS) - len(failed)}/{len(CHECKS)} invariants passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
