"""Tests for the SUNRISE model.

Run with ``python -m pytest tests/`` from the repository root, or
``python tests/test_model.py`` without pytest installed.

These check internal consistency and the structural results, not agreement
with the published figures -- that is what ``scripts/validate.py`` does.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sunrise import analysis, metrics, parameters as prm  # noqa: E402
from sunrise.model import cpi_attachment, drug, kill_rate, simulate  # noqa: E402


def test_every_parameter_has_a_source():
    allowed = {"measured", "literature", "derived", "assumed"}
    for name, p in prm.flat().items():
        assert p.source in allowed, f"{name} has source {p.source!r}"
        assert p.unit, f"{name} has no unit"


def test_sweep_ranges_contain_their_value():
    for name, (lo, hi) in prm.sweep_ranges().items():
        v = prm.flat()[name].value
        assert lo <= v <= hi, f"{name}: {v} outside [{lo}, {hi}]"


def test_overrides_reject_typos():
    try:
        prm.with_overrides(k_inn=1.0)
    except KeyError:
        return
    raise AssertionError("a misspelled parameter was accepted")


def test_drug_profile_accumulates_then_stops():
    p = prm.values()
    first = drug(0.0, p)
    later = drug(5 * p["tau_a"], p)
    assert later > first, "repeated doses should accumulate"
    # Dosing stops after n_dose_a, so the profile decays thereafter.
    end = p["n_dose_a"] * p["tau_a"]
    assert drug(end + 24.0, p) < drug(end, p)


def test_kill_rate_is_monotone_and_saturates():
    p = prm.values()
    xs = [1.0, 10.0, 100.0, 1000.0, 1e5]
    ks = [kill_rate(x, p["mic_s"], p) for x in xs]
    assert all(b >= a for a, b in zip(ks, ks[1:])), "kill rate must be monotone"
    ceiling = p["psi_max"] - p["psi_min"]
    assert ks[-1] <= ceiling + 1e-9
    assert ks[-1] > 0.99 * ceiling, "should approach the ceiling"


def test_susceptible_mic_is_lower_than_resistant():
    from sunrise.model import mic_of
    p = prm.values()
    assert mic_of(0.0, p) == p["mic_s"]
    assert mic_of(p["b_ref"], p) > mic_of(0.0, p)


def test_untreated_population_reaches_carrying_capacity():
    p = prm.values()
    t, ex, inn, _ = simulate(p, a_peak=0.0, dose=0.0)
    assert ex[-1] > 0.5 * p["k_ex"], "untreated should approach capacity"


def test_state_stays_non_negative():
    t, ex, inn, state = simulate()
    assert (state >= -1e-9).all(), "no compartment may go negative"


def test_intracellular_never_exceeds_capacity():
    p = prm.values()
    t, ex, inn, _ = simulate(p)
    assert inn.max() <= p["k_in"] * 1.02, "capacity gate should hold"


def test_combination_beats_antibiotic_alone():
    alone = metrics.evaluate(dose=0.0)
    both = metrics.evaluate()
    assert both["final_total"] < alone["final_total"]
    assert both["depth_log"] > alone["depth_log"]


def test_sunrise_alone_does_nothing():
    """Degrading the resistance factor is useless without the antibiotic."""
    none = metrics.evaluate(a_peak=0.0, dose=0.0)
    sun = metrics.evaluate(a_peak=0.0)
    assert abs(np.log10(sun["final_total"] / none["final_total"])) < 0.05


def test_attachment_closed_form_matches_simulation():
    """The TPD load at the first dose should match the closed form."""
    p = prm.values()
    a = cpi_attachment(p)
    t, _, _, state = simulate(p)
    peak = state[5][(t > 0) & (t < p["tau_p"])].max()
    assert 0.5 * a["tpd_load"] < peak < 2.0 * a["tpd_load"]


def test_eradication_is_blocked_by_the_intracellular_compartment():
    """With the default intracellular growth rate, d > 0 and nothing flips it."""
    dc = analysis.decay_condition()
    assert dc["d"] > 0
    assert not dc["can_be_eradicated"]
    assert dc["blocked_by_intracellular"]


def test_purely_extracellular_levers_cannot_flip_the_sign():
    """Immune clearance, washout and invasion do not appear in d, at any value."""
    base = analysis.decay_condition()["lambda_max"]
    for name, huge in (("k_imm", 50.0), ("k_wash", 50.0), ("k_invade", 0.05)):
        lam = analysis.decay_condition(**{name: huge})["lambda_max"]
        assert lam > 0, f"{name} should not make eradication possible"
        assert abs(lam - base) < 0.05 * abs(base) + 0.02


def test_the_antibiotic_reaches_d_but_only_at_an_impractical_dose():
    """a_peak enters d through eta, so it is not a purely extracellular lever."""
    p = prm.values()
    need = analysis.required_a_peak(p)
    assert np.isfinite(need), "the concentration axis should reach d eventually"
    assert need > 3 * p["a_ceiling"], "and should need far more than is safe"
    assert analysis.decay_condition(a_peak=need * 1.2)["d"] < 0


def test_quiescent_reservoir_interpretation_flips_the_sign():
    """At the measured quiescent growth rate, d < 0 and the system decays."""
    dc = analysis.decay_condition(psi_in=0.002, k_in=5.0)
    assert dc["d"] < 0
    assert dc["can_be_eradicated"]


def test_lambda_is_independent_of_the_intracellular_capacity():
    """k_in is a capacity, not a rate: the exponent must not depend on it."""
    lams = [analysis.lambda_numeric(k_in=k) for k in (500.0, 50.0, 5.0)]
    assert max(lams) - min(lams) < 0.01, f"lambda drifted with k_in: {lams}"


def test_floor_scales_with_the_intracellular_capacity():
    """The durability endpoint is proportional to k_in, which is the warning."""
    ratios = []
    for k in (500.0, 50.0, 5.0):
        r = metrics.evaluate(k_in=k)
        ratios.append(r["final_total"] / k)
    assert max(ratios) / min(ratios) < 1.1, f"not proportional: {ratios}"


def test_depth_and_durability_can_disagree():
    """A regimen can go deep and rebound; that is why both are reported.

    Ver.22 raised MIC_R(0) from 125 to 163.2 ug/mL, and at that MIC kanamycin
    alone no longer clears the 3-log depth bar at the baseline 300 ug/mL -- it
    reaches 2.74 log and fails BOTH endpoints. The disagreement the two
    endpoints exist to expose is still there, but it now needs a concentration
    at the clinical ceiling to show it.
    """
    ceiling = metrics.evaluate(dose=0.0, a_peak=440.0)
    assert ceiling["pass_depth"], "at the ceiling, depth should be reached"
    assert not ceiling["pass_durability"], "and durability should still fail"

    baseline = metrics.evaluate(dose=0.0)
    assert not baseline["pass_depth"], (
        "Ver.22 regression check: at MIC_R(0) = 163.2 kanamycin alone must FAIL "
        "depth at 300 ug/mL. If this passes, MIC_R(0) has moved back.")


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for fn in fns:
        try:
            fn()
            print(f"  ok    {fn.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL  {fn.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            failed += 1
            print(f"  ERROR {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    raise SystemExit(1 if failed else 0)
