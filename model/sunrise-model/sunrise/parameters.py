"""Parameters for the SUNRISE within-host model, each with its provenance.

Every parameter carries a ``source`` field with one of four values:

    measured    a value measured by this project's wet lab
    literature  a value taken from a published measurement
    derived     computed from other entries here, with the derivation stated
    assumed     no measurement exists; treated as a sweep axis

Nothing in this file is a bare number. If you are adapting this model to a
different organism, drug or resistance factor, this is the only file you
should need to rewrite -- ``model.py`` reads everything from here.

Units are stated per entry and are consistent throughout: time in hours,
extracellular and intracellular bacterial density in CFU/mL, drug in ug/mL,
resistance factor in uM, CPI particles in particles, TPD in molecules
per bacterium.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from math import log


@dataclass(frozen=True)
class P:
    """One parameter: value, unit, provenance, and the range to sweep it over."""
    value: float
    unit: str
    source: str              # measured / literature / derived / assumed
    note: str = ""
    low: float | None = None   # lower end of the plausible range
    high: float | None = None  # upper end

    def __float__(self) -> float:
        return float(self.value)


# --------------------------------------------------------------------------
# Pharmacodynamics -- how kanamycin kills, and how the resistance factor
# shifts that curve.
# --------------------------------------------------------------------------
PD = {
    "psi_max": P(1.16, "1/h", "literature",
                 "Maximum growth rate of UPEC in urine. "
                 "mBio 2018, 10.1128/mBio.00186-18 reports a doubling time of "
                 "36.0 +- 2.6 min for CFT073 in pooled human urine ex vivo."),
    "psi_min": P(-0.54, "1/h", "literature",
                 "Maximum net kill rate at saturating kanamycin. "
                 "Lee & Lee 2024 time-kill curves."),
    "kappa": P(2.95, "-", "literature",
               "Hill exponent of the dose-response in A/MIC. Inherited from the "
               "pre-recalibration fit (n_1 = 2.95). See docs/parameters.md: the "
               "recalibrated fit gives 1.73 for the susceptible strain and 3.05 "
               "for the resistant strain, which a single exponent cannot represent.",
               low=1.7, high=3.1),
    "mic_s": P(1.20, "ug/mL", "measured",
               "MIC of the susceptible strain. The recalibrated pharmacodynamic "
               "parameters of the standalone model imply 1.11; the two differ by 8%.",
               low=1.11, high=1.60),
    "mic_r0": P(163.2, "ug/mL", "measured",
                "MIC of the resistant strain carrying a full complement of the "
                "resistance factor. Ver.22 adopts the value implied by the "
                "recalibrated standalone parameters (163.2); the earlier direct "
                "measurement was 125, and refitting this functional form to the "
                "recalibrated curve gives 191.",
                low=125.0, high=191.0),
    "b_ref": P(0.3233, "uM", "derived",
               "Resistance factor concentration in the untreated resistant strain, "
               "converted from the detox capacity of the standalone model."),
}

# --------------------------------------------------------------------------
# Antibiotic pharmacokinetics in the bladder (intravesical or systemic;
# the model only needs the urinary concentration profile).
# --------------------------------------------------------------------------
PK_DRUG = {
    "a_peak": P(300.0, "ug/mL", "assumed",
                "Peak urinary kanamycin concentration. Design variable.",
                low=50.0, high=700.0),
    "k_elim": P(0.30, "1/h", "literature",
                "First-order decay of the urinary concentration between doses."),
    "tau_a": P(12.0, "h", "assumed", "Antibiotic dosing interval.",
               low=6.0, high=24.0),
    "n_dose_a": P(6, "-", "assumed", "Number of antibiotic doses in the course."),
    "a_ceiling": P(440.0, "ug/mL", "literature",
                   "Clinical ceiling on the urinary concentration, derived from the "
                   "labelled serum peak times the urine/serum ratio. This is a "
                   "SYSTEMIC-dosing ceiling: for intravesical instillation it must be "
                   "re-derived from the fraction absorbed."),
}

# --------------------------------------------------------------------------
# CPI particle pharmacokinetics -- encounter, loss, and payload delivery.
# --------------------------------------------------------------------------
PK_CPI = {
    "volume": P(200.0, "mL", "assumed", "Bladder volume."),
    "k_enc": P(1.78e-7, "mL/h", "derived",
               "Smoluchowski encounter rate coefficient between one particle and "
               "one bacterium."),
    "k_loss": P(0.35, "1/h", "assumed", "Particle inactivation rate.",
                low=0.1, high=1.0),
    "phi_p": P(20.3, "molecules/particle", "derived",
               "EXPECTED DELIVERY per encounter: the length-weighted product of "
               "the geometric penetration probability phi(L) and the loading "
               "capacity N_max(L), summed over the measured length distribution "
               "(32,800 particles, WLC, 4 nm bins). Geometrical condition only -- "
               "mechanical puncture efficiency and actual loading/release "
               "efficiency are NOT included, so this is an UPPER BOUND, and the "
               "range spans structural scenarios only (entry angle 0-30 deg, "
               "125/133 residues, 0.4/0.5 nm discretisation). "
               "This must stay ONE parameter: phi and N_max covary positively "
               "with L, so the product of the marginals understates the true "
               "value by 4-18% (16% at psi=0 deg: 0.4366 x 40.07 = 17.5 against "
               "a true 20.29). The marginal phi is 0.4366 (0 deg) / 0.4681 (10) / "
               "0.5763 (20) / 0.7697 (30); it is NOT 0.98.",
               low=10.7, high=34.9),
    "eta_mech": P(1.0, "-", "assumed",
                  "Mechanical puncture x loading/release efficiency, the factor "
                  "phi_p excludes. N_0 = eta_mech * phi_p * m, so the required "
                  "attachment m_inf, the dose and the particles per course all "
                  "scale as 1/eta_mech. 1.0 is the upper-bound limit, not an "
                  "estimate.",
                  low=0.2, high=1.0),
    "dose": P(3.0e11, "particles/dose", "assumed",
              "SUNRISE dose. Design variable.", low=1e10, high=1e13),
    "tau_p": P(6.0, "h", "assumed", "SUNRISE dosing interval. Design variable.",
               low=3.0, high=12.0),
    "delta_t": P(0.0, "h", "assumed",
                 "Phase offset of the CPI schedule relative to the antibiotic "
                 "schedule, taken modulo tau_p. 0 means the two are given "
                 "together. Design variable; Table 9.4.30 reports the "
                 "worst-case phase at each tau_p.",
                 low=0.0, high=12.0),
}

# --------------------------------------------------------------------------
# Targeted protein degradation inside the bacterium.
# --------------------------------------------------------------------------
TPD = {
    "k_tpd0": P(0.14 * 60, "1/h", "assumed",
                "Degradation rate at the reference TPD load. Calibrated in "
                "validate.py against the published baseline; the nominal value "
                "here is a placeholder. The standalone degradation model's closed "
                "form gives 0.041-0.091 /min over the relevant load range, i.e. "
                "the placeholder is optimistic by about 1.5x.",
                low=0.0338 * 60, high=0.2516 * 60),
    "n0_ref": P(231.0, "molecules/bacterium", "derived",
                "Reference TPD load, = phi_p * m(0) at the published "
                "calibration conditions (16.0 x 14.458 = 231). FIXED: it is the "
                "TPD module's normalisation constant, not a per-trial value. "
                "Substituting the trial's own load here would cancel the dose "
                "dependence of k_tpd out of the model entirely."),
    "k_xb": P(9.8185e-4 * 60, "1/h", "measured",
              "Basal degradation rate conferred by the degron tag alone "
              "(pBAD chase experiment)."),
    "mu_ref": P(log(2) / 46 * 60, "1/h", "measured",
                "Growth rate at which the resistance factor production rate was "
                "calibrated (doubling time 46 min in culture)."),
    "a_interp": P(0.5, "-", "assumed",
                  "How strongly resistance factor production tracks growth rate. "
                  "0 = constant production, 1 = constant concentration.",
                  low=0.0, high=1.0),
    "k_prot": P(0.0, "1/h", "assumed",
                "Proteolytic loss of TPD inside the bacterium.", low=0.0, high=0.1),
}

# --------------------------------------------------------------------------
# Population dynamics -- four compartments (susceptible/resistant x
# extracellular/intracellular).
# --------------------------------------------------------------------------
POP = {
    "k_ex": P(1.0e9, "CFU/mL", "measured",
              "Extracellular carrying capacity. The standalone model measured "
              "6.29e8 at the drug-free endpoint.",
              low=6.29e8, high=1.0e9),
    "k_in": P(50.0, "CFU/mL", "literature",
              "Intracellular carrying capacity = total intracellular burden per "
              "bladder / volume. Direct measurements give 1e4 CFU/bladder at the "
              "acute stage (Infection and Immunity 2011, 10.1128/IAI.05339-11) and "
              "~1e3 at the chronic stage, i.e. 50 and 5 CFU/mL in 200 mL. "
              "The acute band spans 3e3-1e5 per bladder.",
              low=15.0, high=500.0),
    "k_imm": P(0.10, "1/h", "assumed",
               "Immune clearance of extracellular bacteria.", low=0.0, high=0.3),
    "k_wash": P(0.30, "1/h", "derived",
                "Washout by urine production and voiding."),
    "k_invade": P(0.01, "1/h", "literature",
                  "Extracellular -> intracellular. Converted from gentamicin "
                  "protection assays: 4.4% of adherent bacteria invaded by 1 h in "
                  "mouse gives 0.045 /h; 0.5-10% of cell-associated over 2 h in "
                  "human 5637 cells gives 0.0025-0.053 /h. All denominators are "
                  "adherent or cell-associated, never the luminal pool, so these "
                  "are upper bounds.",
                  low=0.001, high=0.05),
    "k_emerge": P(0.065, "1/h", "literature",
                  "Intracellular -> extracellular, inferred as 1/residence time "
                  "from the IBC cycle timeline (invasion ~1 h, flux ~16 h in mouse; "
                  "first release ~15.6 h on a human bladder chip, eLife 2021, "
                  "10.7554/eLife.66481). Never measured as a rate.",
                  low=0.04, high=0.10),
    "k_exf": P(0.0, "1/h", "assumed",
               "Loss of intracellular bacteria to urothelial exfoliation and "
               "autophagy. No study reports this as a rate. Zero is both the "
               "conservative choice and arguably the best point estimate for "
               "wild-type UPEC, which neutralises the lysosome and is expelled "
               "alive. The only derivable anchor is 0.043 /h.",
               low=0.0, high=0.1),
    "psi_in": P(1.3, "1/h", "literature",
                "Intracellular growth rate. Ver.22 FIXES the intracellular "
                "compartment as the acute-phase IBC (intracellular bacterial "
                "community), NOT the quiescent reservoir: the 72 h uncomplicated "
                "cystitis course is an acute-phase setting, and k_in and "
                "k_emerge are already IBC-derived. 1.3 /h is the median of the "
                "measured 1.19-1.39 /h (doubling 30-35 min, PNAS 2004, "
                "10.1073/pnas.0308125100). The sweep band 0.69-1.39 runs from the "
                "mid-IBC upper bound to the early-IBC measurement. Quiescent "
                "reservoirs (<0.002 /h) are OUT of scope for this course. "
                "See analysis.decay_condition().",
                low=0.69, high=1.39),
    "eta": P(0.05, "-", "assumed",
             "Effective intracellular drug concentration as a fraction of the "
             "extracellular one. A lumped engineering approximation, not a "
             "measurement: intracellularly the EC50 is roughly preserved while "
             "Emax collapses and the Hill slope flattens ~10-fold, so a single "
             "multiplier mis-extrapolates across a wide concentration sweep.",
             low=0.005, high=0.1),
    "rho_in": P(1.0, "-", "assumed",
                "Fraction of the resistance phenotype retained intracellularly. "
                "1.0 means TPD does not reach intracellular bacteria."),
    "gamma": P(1.0e-11, "mL/(CFU h)", "assumed",
               "Conjugative transfer of the resistance determinant.",
               low=1e-13, high=1e-9),
    "s0": P(1.0e8, "CFU/mL", "assumed", "Initial susceptible density."),
    "r0": P(1.0e5, "CFU/mL", "assumed", "Initial resistant density."),
}

# --------------------------------------------------------------------------
# Endpoints.
# --------------------------------------------------------------------------
ENDPOINT = {
    "t_end": P(72.0, "h", "assumed", "Length of the treatment course."),
    "cure_density": P(1.0e2, "CFU/mL", "literature",
                      "Density below which the course is called a success. "
                      "10^2 CFU/mL is a DIAGNOSTIC threshold for symptomatic "
                      "cystitis with pyuria (N Engl J Med 2013, "
                      "10.1056/NEJMoa1302186), not a cure threshold. The "
                      "regulatory analogue is <10^3 CFU/mL at a post-therapy "
                      "test-of-cure visit.",
                      low=1.0e2, high=1.0e3),
    "nadir_target": P(3.0, "log10", "literature",
                      "Required log reduction at the nadir. Stasis / 1-log / "
                      "2-log / 3-log is the standard pharmacodynamic ladder."),
}

ALL = {"PD": PD, "PK_DRUG": PK_DRUG, "PK_CPI": PK_CPI, "TPD": TPD,
       "POP": POP, "ENDPOINT": ENDPOINT}


def flat() -> dict[str, P]:
    """Every parameter in one dict, keyed by name."""
    out: dict[str, P] = {}
    for group in ALL.values():
        out.update(group)
    return out


def values() -> dict[str, float]:
    """Every parameter as a plain float, for passing into the model."""
    return {k: float(v) for k, v in flat().items()}


def with_overrides(**kw) -> dict[str, float]:
    """Parameter values with some replaced, e.g. ``with_overrides(k_in=500.0)``."""
    v = values()
    unknown = set(kw) - set(v)
    if unknown:
        raise KeyError(f"unknown parameter(s): {sorted(unknown)}")
    v.update(kw)
    return v


def by_source(source: str) -> list[str]:
    """Names of every parameter with the given provenance."""
    return sorted(k for k, p in flat().items() if p.source == source)


def sweep_ranges() -> dict[str, tuple[float, float]]:
    """The plausible range of every parameter that has one."""
    return {k: (p.low, p.high) for k, p in flat().items()
            if p.low is not None and p.high is not None}


def table(fmt: str = "markdown") -> str:
    """The provenance table, for pasting into a wiki or a report."""
    rows = [("Parameter", "Value", "Unit", "Source", "Range")]
    for group_name, group in ALL.items():
        for name, p in group.items():
            rng = (f"{p.low:g} - {p.high:g}" if p.low is not None else "")
            rows.append((name, f"{p.value:g}", p.unit, p.source, rng))
    if fmt != "markdown":
        return "\n".join("\t".join(r) for r in rows)
    w = [max(len(r[i]) for r in rows) for i in range(5)]
    out = ["| " + " | ".join(h.ljust(w[i]) for i, h in enumerate(rows[0])) + " |",
           "|" + "|".join("-" * (w[i] + 2) for i in range(5)) + "|"]
    out += ["| " + " | ".join(c.ljust(w[i]) for i, c in enumerate(r)) + " |"
            for r in rows[1:]]
    return "\n".join(out)


if __name__ == "__main__":
    print(table())
    print()
    for s in ("measured", "literature", "derived", "assumed"):
        names = by_source(s)
        print(f"{s:<11} {len(names):2d}  {', '.join(names)}")
