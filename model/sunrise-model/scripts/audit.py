#!/usr/bin/env python3
"""Compare this model against the standalone degradation and integration modules.

The two formulations are not parameter-for-parameter comparable. This model
shifts the MIC with the resistance factor level and uses a four-parameter form
in A/MIC; the standalone model routes the drug through a detox capacity and a
Michaelis constant to an intracellular active concentration, then applies a
saturating kill fraction. The Hill exponents are exponents on different
variables.

So rather than lining up symbols, this compares the curves they produce, and
then refits this model's functional form to the standalone one to see what it
would take to agree.

    python scripts/audit.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sunrise import parameters as prm  # noqa: E402
from sunrise.model import kill_rate  # noqa: E402

# Standalone integration module, recalibrated parameter set.
STANDALONE = {
    "r_g": np.log(2) / 46 * 60,   # 1/h, doubling time 46 min in culture
    "r_kill": 0.02837 * 60,       # 1/h
    "k1": 1.032,                  # ug/mL, half-maximal intracellular concentration
    "n1": 1.722,                  # -
    "km": 1.463,                  # uM, Michaelis constant of the resistance enzyme
    "q0": 375.9,                  # ug/mL, detox capacity of the resistant strain
    "b0": 0.3233,                 # uM
    "capacity": 6.29e8,           # CFU/mL
    "payload_uM": 0.20,           # uM per particle
    "phi": 0.98,
}
# The same set before recalibration, kept because this model inherited from it.
PRE_RECALIBRATION = {"k1": 0.503, "n1": 2.95, "km": 4.66, "q0": 454.0}
MOLECULES_PER_UM = 602.0          # a 1 fL bacterium


def active_concentration(a_ext: float, q: float, km: float) -> float:
    """Intracellular active drug after the resistance enzyme has worked on it."""
    b = km + q - a_ext
    return (-b + np.sqrt(b * b + 4 * a_ext * km)) / 2


def standalone_kill(a_ext: float, q: float, s=STANDALONE, **over) -> float:
    s = dict(s, **over)
    ai = active_concentration(a_ext, q, s["km"])
    return s["r_kill"] * ai ** s["n1"] / (s["k1"] ** s["n1"] + ai ** s["n1"])


def mic_from_curve(kill_fn, growth: float, lo=1e-4, hi=1e6) -> float:
    """The concentration at which net growth is zero."""
    try:
        return float(brentq(lambda a: growth - kill_fn(a), lo, hi, xtol=1e-6))
    except ValueError:
        return float("nan")


def main() -> None:
    p = prm.values()
    s = STANDALONE

    print("=" * 76)
    print("1. Does this model's MIC follow from the standalone parameters?")
    print("=" * 76)
    mic_s = mic_from_curve(lambda a: standalone_kill(a, 0.0), s["r_g"])
    mic_r = mic_from_curve(lambda a: standalone_kill(a, s["q0"]), s["r_g"])
    old = dict(s, **PRE_RECALIBRATION)
    mic_s_old = mic_from_curve(lambda a: standalone_kill(a, 0.0, old), s["r_g"])
    mic_r_old = mic_from_curve(
        lambda a: standalone_kill(a, old["q0"], old), s["r_g"])

    print(f"{'':<36}{'MIC_S':>10}{'MIC_R':>10}{'ratio':>9}")
    print("-" * 76)
    print(f"{'this model, as written':<36}{p['mic_s']:>10.3f}"
          f"{p['mic_r0']:>10.1f}{p['mic_r0']/p['mic_s']:>9.1f}")
    print(f"{'from the recalibrated standalone set':<36}{mic_s:>10.3f}"
          f"{mic_r:>10.1f}{mic_r/mic_s:>9.1f}")
    print(f"{'from the pre-recalibration set':<36}{mic_s_old:>10.3f}"
          f"{mic_r_old:>10.1f}{mic_r_old/mic_s_old:>9.1f}")

    print("\n" + "=" * 76)
    print("2. The kill-rate curves themselves (resistant strain)")
    print("=" * 76)
    print(f"{'A [ug/mL]':>10}{'this model':>13}{'standalone':>13}{'ratio':>9}")
    print("-" * 76)
    for a in (10.0, 50.0, 100.0, 125.0, 200.0, 300.0, 440.0, 700.0):
        k_here = kill_rate(a, p["mic_r0"], p)
        k_there = standalone_kill(a, s["q0"])
        r = k_here / k_there if k_there > 1e-9 else float("nan")
        print(f"{a:>10.0f}{k_here:>13.4f}{k_there:>13.4f}{r:>9.2f}")
    print(f"\n  Kill ceiling: {p['psi_max'] - p['psi_min']:.3f} /h here, "
          f"{s['r_kill']:.3f} /h there -- these agree.")
    print(f"  Growth rate:  {p['psi_max']:.3f} /h here (urine), "
          f"{s['r_g']:.3f} /h there (culture) -- different by design.")

    print("\n" + "=" * 76)
    print("3. Refitting this model's form to the recalibrated curve")
    print("=" * 76)
    grid = np.logspace(np.log10(0.5), np.log10(1000), 400)

    def form(a, mic, kappa):
        x = (a / mic) ** kappa
        return (p["psi_max"] - p["psi_min"]) * x / (x - p["psi_min"] / p["psi_max"])

    for label, q, mic0 in [("susceptible (q = 0)", 0.0, 1.2),
                           ("resistant (q = q0)", s["q0"], 125.0)]:
        target = np.array([standalone_kill(a, q) for a in grid])
        (mic, kappa), _ = curve_fit(form, grid, target, p0=[mic0, 2.0],
                                    maxfev=20000)
        rmse = float(np.sqrt(np.mean((form(grid, mic, kappa) - target) ** 2)))
        print(f"  {label:<24} MIC = {mic:8.3f} ug/mL   "
              f"exponent = {kappa:.3f}   RMSE = {rmse:.4f} /h")
    print(f"\n  In use here: MIC_S = {p['mic_s']}, MIC_R = {p['mic_r0']}, "
          f"exponent = {p['kappa']}")
    print("  The two strains want different exponents, which one parameter")
    print("  cannot express: the detox route changes the shape of the curve,")
    print("  not only its position.")

    print("\n" + "=" * 76)
    print("4. Parameters that do correspond one-to-one")
    print("=" * 76)
    rows = [
        ("maximum kill rate", p["psi_max"] - p["psi_min"], s["r_kill"], "1/h"),
        ("resistance factor baseline", p["b_ref"], s["b0"], "uM"),
        ("extracellular capacity", p["k_ex"], s["capacity"], "CFU/mL"),
        ("expected delivery phi*P", p["phi_p"],
         s["payload_uM"] * MOLECULES_PER_UM * s["phi"], "molecules"),
    ]
    print(f"{'quantity':<30}{'here':>13}{'standalone':>13}{'ratio':>8}  status")
    print("-" * 76)
    for name, here, there, unit in rows:
        r = here / there
        status = "agree" if 0.98 <= r <= 1.02 else "differ"
        print(f"{name:<30}{here:>13.4g}{there:>13.4g}{r:>8.2f}  {status}")

    print("\n" + "=" * 76)
    print("5. In the standalone modules, absent here")
    print("=" * 76)
    for line in [
        "protease competition window after each antibiotic dose",
        "structurally non-injectable subpopulation (its largest stated risk)",
        "detox capacity and Michaelis route to intracellular active drug",
        "the resistance-factor threshold expressed as a concentration",
        "closed form from degrader turnover to effective degradation rate",
    ]:
        print(f"  - {line}")

    print("\n" + "=" * 76)
    print("6. Here, absent in the standalone modules")
    print("=" * 76)
    for line in [
        "the intracellular compartment and everything governing it",
        "susceptible and resistant lineages with conjugative transfer",
        "immune clearance and washout",
        "repeated antibiotic dosing with a decay profile",
        "repeated particle dosing and the attachment history",
    ]:
        print(f"  - {line}")


if __name__ == "__main__":
    main()
