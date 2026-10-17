"""Parameters for the SUNRISE-2 bladder pharmacodynamic model.

All values follow SUNRISE-2 Ver.23 (Table 9.4.5 / 9.4.7 / 9.4.10 / 9.4.14 /
9.4.22).  Rates are per hour unless the name says otherwise.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

# --- unit helpers -----------------------------------------------------------

PER_MIN = 60.0  # multiply a /min rate by this to get /h

# --- antibiotic dose conversion (Ver.23, 抗生物質の扱い) --------------------
# 1 dose of 0.45 g (7.5 mg/kg IM, 60 kg) = 22 ug/mL serum peak
#                                        = 300 ug/mL model-internal urine
DOSE_REF_G = 0.45
URINE_REF = 300.0  # ug/mL
SERUM_REF = 22.0  # ug/mL
SERUM_PEAK_MAX = 35.0  # ug/mL, safety ceiling
DAILY_TOTAL_MAX_G = 1.5  # g/day, all routes combined


def dose_to_urine(dose_g: float) -> float:
    """Per-dose mass (g) -> model-internal urinary peak concentration."""
    return dose_g / DOSE_REF_G * URINE_REF


def urine_to_dose(urine: float) -> float:
    """Model-internal urinary peak concentration -> per-dose mass (g)."""
    return urine / URINE_REF * DOSE_REF_G


def dose_to_serum(dose_g: float) -> float:
    """Per-dose mass (g) -> serum peak concentration."""
    return dose_g / DOSE_REF_G * SERUM_REF


def serum_to_dose(serum: float) -> float:
    return serum / SERUM_REF * DOSE_REF_G


#: per-dose ceiling implied by the serum peak limit (0.7159 g, quoted as 0.72)
DOSE_MAX_G = serum_to_dose(SERUM_PEAK_MAX)
#: urinary-equivalent ceiling (480 ug/mL)
URINE_MAX = dose_to_urine(DOSE_MAX_G)


@dataclass(frozen=True)
class Params:
    """Full parameter set.  Defaults are the Ver.23 standard patient."""

    # ---- course / judgement -------------------------------------------
    t_end: float = 72.0  # h, course length
    tau_A: float = 12.0  # h, antibiotic dosing interval
    cure_density: float = 1e2  # CFU/mL, Cbar^tot target
    delta_log_target: float = 3.0

    # ---- CPI pharmacokinetics (Table 9.4.5) ---------------------------
    V: float = 200.0  # mL, bladder volume
    k_enc: float = 1.78e-7  # mL/h, Smoluchowski encounter rate
    k_loss_cpi: float = 0.35  # /h, loss of injection competence in urine
    nu: float = 32.3  # molecules per ADHERED particle (expected delivery)
    phi_bar: float = 0.770  # marginal penetration probability
    D: float = 1.0e11  # particles/dose (design variable)
    tau_P: float = 6.0  # h, SUNRISE dosing interval
    dt_phase: float = 0.0  # h, phase offset of SUNRISE vs antibiotic
    q_urine: float = 60.0  # mL/h, urine production
    V_max: float = 300.0  # mL, volume at onset of voiding
    V_r: float = 125.0  # mL, residual volume
    k_washout: float = 0.30  # /h, dilution (derived from V_r)

    # ---- pharmacodynamics (Table 9.4.7) -------------------------------
    psi_max: float = 1.16  # /h
    psi_min: float = -0.54  # /h
    kappa: float = 2.95
    MIC_S: float = 1.20  # ug/mL
    MIC_R0: float = 163.2  # ug/mL
    B_ref: float = 0.3233  # uM
    h_mic: float = 1.0  # curvature of MIC vs resistance-factor load
    rho_in: float = 1.0  # intracellular resistance-factor survival
    eta: float = 0.05  # intracellular drug discount
    k_elim: float = 0.30  # /h
    A_peak: float = URINE_REF  # ug/mL urinary peak (= 0.45 g/dose)

    # ---- innate immunity (Table 9.4.10) -------------------------------
    induced_immunity: bool = False  # baseline fixes N at rest
    k_phag: float = 0.20  # /h, = k_PMN * N_PMN at rest
    K_cap: float = 1.0e8  # CFU/mL
    K_G: float = 1.0e6  # CFU/mL
    sigma_in: float = 0.30
    d_G: float = 1.00  # /h
    a_N: float = 0.90  # /h
    d_N: float = 0.030  # /h
    N_max: float = 10.0
    eps_urine: float = 1.0
    N_PMN_rest: float = 1.0e5  # neutrophils/mL at rest
    N_PMN: float = 9.0e5  # neutrophils/mL, standard patient (900 /uL)
    k_imm_rest: float = 0.10  # /h, rest-state removal rate

    # ---- population dynamics (Table 9.4.14) ---------------------------
    k_TPD0_min: float = 0.0204  # /min at A_cell = N0_ref
    N0_ref: float = 231.0  # molecules/cell, fixed upstream reference
    window_lo: float = 93.0  # molecules/cell, effective window
    window_hi: float = 1991.0  # molecules/cell
    k_XB_min: float = 9.8185e-4  # /min, endogenous ClpXP baseline
    k_prot_min: float = 0.0  # /min, intrabacterial TPD degradation
    a_prod: float = 0.5  # growth dependence of resistance-factor output
    mu_ref_doubling_min: float = 46.0  # min, calibration doubling time
    K_ex: float = 1.0e9  # CFU/mL
    K_in: float = 50.0  # CFU/mL
    k_exf: float = 0.0  # /h
    k_invade: float = 0.01  # /h
    k_emerge: float = 0.065  # /h
    psi_max_in: float = 1.3  # /h, early IBC
    gamma_conj: float = 1.0e-11  # mL/cell/h
    S0: float = 0.0  # CFU/mL, standard patient is 100 % resistant
    R0: float = 1.0e8  # CFU/mL

    # ---- switches ------------------------------------------------------
    sunrise_on: bool = True
    antibiotic_on: bool = True
    cap_k_TPD_at_window: bool = False  # clamp k_TPD above the window top

    # ---- derived -------------------------------------------------------
    @property
    def k_TPD0(self) -> float:
        return self.k_TPD0_min * PER_MIN

    @property
    def k_XB(self) -> float:
        return self.k_XB_min * PER_MIN

    @property
    def k_prot(self) -> float:
        return self.k_prot_min * PER_MIN

    @property
    def mu_ref(self) -> float:
        """Division rate at the TPD-module calibration condition, /h."""
        return math.log(2.0) / self.mu_ref_doubling_min * PER_MIN

    @property
    def k_PMN(self) -> float:
        """Per-neutrophil killing rate implied by the rest-state calibration."""
        return self.k_phag / self.N_PMN_rest

    @property
    def C0_ex(self) -> float:
        return self.S0 + self.R0

    @property
    def total_CFU0(self) -> float:
        """Initial bacterial burden per bladder."""
        return self.C0_ex * self.V

    @property
    def dose_g(self) -> float:
        return urine_to_dose(self.A_peak)

    @property
    def daily_g(self) -> float:
        return self.dose_g * (24.0 / self.tau_A)

    @property
    def serum_peak(self) -> float:
        return dose_to_serum(self.dose_g)

    @property
    def n_doses_P(self) -> int:
        """Number of SUNRISE doses inside the course."""
        if not self.sunrise_on:
            return 0
        return int(math.ceil((self.t_end - self.dt_phase) / self.tau_P))

    @property
    def MOI(self) -> float:
        """Particles per CFU at t = 0."""
        return self.D / self.total_CFU0

    def with_(self, **kw) -> "Params":
        return replace(self, **kw)

    def with_MOI(self, moi: float) -> "Params":
        return replace(self, D=moi * self.total_CFU0)

    def with_dose_g(self, dose_g: float) -> "Params":
        return replace(self, A_peak=dose_to_urine(dose_g))


def k_washout_from_residual(V_r: float, q: float = 60.0, V_max: float = 300.0) -> float:
    """Table 9.4.3: cycle-averaged dilution rate from the residual volume."""
    T = (V_max - V_r) / q
    return math.log(V_max / V_r) / T


STANDARD_PATIENT = Params()

#: competitive-release scenario kept separate from the main line (Ver.23)
COMPETITIVE_RELEASE = Params(S0=1.0e8, R0=1.0e5)
