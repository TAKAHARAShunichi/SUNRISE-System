"""Tasks 4-8: tolerance limits, the IBC band, burden dependence, the zero
class, and the cross-checks."""

from __future__ import annotations

import csv
import json
import math
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from sunrise import analytic
from sunrise.design import N_OUT, bisect_max_log, bisect_min_log, moi_cap
from sunrise.metrics import evaluate
from sunrise.model import kill, m_infinity, psi, solve_course
from sunrise.params import DOSE_MAX_G, STANDARD_PATIENT, Params
from sunrise.zeroclass import solve_course_z

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

#: reference design used for the tolerance table: tau_P = 6 h, dose at the cap,
#: MOI from the worst-phase search, with the phase at its optimum
REF_TAU_P = 6.0
REF_MOI = 1.409
REF_PHASE = 3.0


def ref_params(**kw) -> Params:
    p = (
        STANDARD_PATIENT.with_(tau_P=REF_TAU_P, dt_phase=REF_PHASE)
        .with_dose_g(DOSE_MAX_G)
        .with_MOI(REF_MOI)
    )
    return p.with_(**kw) if kw else p


def verdict(p: Params) -> bool:
    return evaluate(solve_course(p, n_out=N_OUT)).passes


# --- task 4: tolerance limits ----------------------------------------------
# k_TPD0 collects the three upstream enzymology axes.  Under a
# Michaelis-Menten reading k_TPD0 ~ k_cat * [ClpXP]_eff / (K_D + B_ref), so a
# threshold in k_TPD0 maps to each of them; the absolute mapping needs the TPD
# module, the fold-change does not.
K_CAT_REF = 120.0  # /h
CLPXP_REF = 0.166  # uM
K_D_REF = 0.154  # uM


def _axis_min(name: str, lo: float, hi: float, build) -> dict:
    """Largest-tolerated / smallest-tolerated value on one axis."""
    base_ok = verdict(build(None))
    up = bisect_max_log(lambda x: verdict(build(x)), lo, hi, tol=2e-3)
    down = bisect_min_log(lambda x: verdict(build(x)), lo, hi, tol=2e-3)
    return {
        "axis": name,
        "baseline_passes": base_ok,
        "lowest_passing": down,
        "highest_passing": up,
        "search_lo": lo,
        "search_hi": hi,
    }


AXES = [
    ("nu [molecules/adhered particle]", 1.0, 200.0,
     lambda x: ref_params() if x is None else ref_params(nu=x)),
    ("k_TPD0 [/min]", 1e-4, 10.0, lambda x: ref_params()
     if x is None else ref_params(k_TPD0_min=x)),
    ("psi_max_in [/h]", 0.0, 5.0, lambda x: ref_params()
     if x is None else ref_params(psi_max_in=x)),
    ("K_in (n_IBC) [CFU/mL]", 1.0, 5000.0, lambda x: ref_params()
     if x is None else ref_params(K_in=x)),
    ("eps_urine = N_PMN/1e5 [-]", 1e-3, 20.0, lambda x: ref_params()
     if x is None else ref_params(eps_urine=x)),
    ("k_washout [/h]", 0.02, 2.0, lambda x: ref_params()
     if x is None else ref_params(k_washout=x)),
    ("MIC_R0 [ug/mL]", 10.0, 3000.0, lambda x: ref_params()
     if x is None else ref_params(MIC_R0=x)),
    ("h (MIC curvature) [-]", 0.1, 3.0, lambda x: ref_params()
     if x is None else ref_params(h_mic=x)),
    ("k_prot [/min]", 1e-5, 1.0, lambda x: ref_params()
     if x is None else ref_params(k_prot_min=x)),
    ("eta [-]", 1e-4, 1.0, lambda x: ref_params()
     if x is None else ref_params(eta=x)),
    ("rho_in [-]", 1e-3, 1.0, lambda x: ref_params()
     if x is None else ref_params(rho_in=x)),
    ("k_loss_cpi [/h]", 0.01, 20.0, lambda x: ref_params()
     if x is None else ref_params(k_loss_cpi=x)),
    ("k_exf [/h]", 1e-4, 1.0, lambda x: ref_params()
     if x is None else ref_params(k_exf=x)),
    ("k_emerge [/h]", 1e-3, 2.0, lambda x: ref_params()
     if x is None else ref_params(k_emerge=x)),
    ("k_invade [/h]", 1e-4, 1.0, lambda x: ref_params()
     if x is None else ref_params(k_invade=x)),
]


def _axis_worker(i: int) -> dict:
    name, lo, hi, build = AXES[i]
    return _axis_min(name, lo, hi, build)


def task4(pool: Pool) -> dict:
    rows = pool.map(_axis_worker, range(len(AXES)))
    k_row = next(r for r in rows if r["axis"].startswith("k_TPD0"))
    k_crit = k_row["lowest_passing"]
    enzym = None
    if k_crit:
        fold = k_crit / STANDARD_PATIENT.k_TPD0_min
        enzym = {
            "k_TPD0_crit_per_min": k_crit,
            "fold_of_baseline": fold,
            "k_cat_crit_per_h": K_CAT_REF * fold,
            "ClpXP_available_fraction_crit": fold,
            "ClpXP_eff_crit_uM": CLPXP_REF * fold,
            "K_D_crit_uM": (K_D_REF + STANDARD_PATIENT.B_ref) / fold
            - STANDARD_PATIENT.B_ref,
            "mapping": (
                "k_TPD0 ~ k_cat * [ClpXP]_eff / (K_D + B_ref); only the "
                "fold-change is model-internal, the absolute values need the "
                "TPD module."
            ),
        }
    return {"reference_design": {
        "tau_P": REF_TAU_P, "MOI": REF_MOI, "dose_g": DOSE_MAX_G,
        "dt_phase": REF_PHASE}, "axes": rows, "enzymology": enzym}


# --- task 5: IBC band ------------------------------------------------------


def _task5_one(psi_in: float) -> dict:
    p = ref_params(psi_max_in=psi_in)
    sol = solve_course(p, n_out=N_OUT)
    v = evaluate(sol)
    k_in_max = bisect_max_log(
        lambda K: verdict(ref_params(psi_max_in=psi_in, K_in=K)), 1.0, 5000.0,
        tol=2e-3,
    )
    return {
        "psi_max_in": psi_in,
        "C_bar": v.C_bar,
        "delta_log": v.delta_log,
        "passes_depth": v.passes_depth,
        "passes_durability": v.passes_durability,
        "passes": bool(v.passes),
        "residual_CFU_bladder_72h": v.residual_CFU,
        "C_bar_over_K_in": v.C_bar / p.K_in,
        "K_in_max_passing": k_in_max,
    }


def task5(pool: Pool) -> dict:
    rows = pool.map(_task5_one, [0.69, 1.19, 1.3, 1.39])
    return {"rows": rows}


# --- task 6: burden dependence --------------------------------------------


def _task6_one(R0: float) -> dict:
    base = (
        STANDARD_PATIENT.with_(tau_P=REF_TAU_P, dt_phase=REF_PHASE, R0=R0)
        .with_dose_g(DOSE_MAX_G)
    )
    moi = bisect_min_log(lambda m: verdict(base.with_MOI(m)), 1e-3, 1e4, tol=2e-3)
    out = {
        "C0_ex": R0,
        "total_CFU0": base.total_CFU0,
        "moi_optimal": moi,
    }
    if moi:
        p = base.with_MOI(moi)
        out.update({
            "D_per_dose": p.D,
            "particles_per_course": p.D * p.n_doses_P,
            "m_inf": m_infinity(p),
        })
    return out


def task6(pool: Pool) -> dict:
    rows = pool.map(_task6_one, [1e6, 3e6, 1e7, 3e7, 1e8])
    for r in rows:
        from sunrise.model import f_hit

        r["f_hit_at_t0"] = f_hit(r["C0_ex"], STANDARD_PATIENT)
    monotone = all(
        (a["moi_optimal"] or 0) <= (b["moi_optimal"] or 0) + 1e-9
        for a, b in zip(rows, rows[1:])
    )
    p = STANDARD_PATIENT
    crossover = (p.k_washout + p.k_loss_cpi) / p.k_enc
    Ds = [r["D_per_dose"] for r in rows if r.get("D_per_dose")]
    return {
        "rows": rows,
        "monotone_increasing_with_burden": monotone,
        "f_hit_crossover_CFU_per_mL": crossover,
        "D_per_dose_span": (min(Ds), max(Ds)) if Ds else None,
        "D_per_dose_fold_over_100x_burden": (max(Ds) / min(Ds)) if Ds else None,
        "interpretation": (
            "Required MOI falls as the burden rises, i.e. the opposite of the "
            "brief's expectation, and it is not an implementation fault. "
            "f_hit = k_enc C / (k_enc C + k_washout + k_loss) crosses over at "
            f"{crossover:.2e} CFU/mL. Below it most particles are washed out "
            "before meeting a cell, so m = D k_enc / (V (k_wash + k_loss)) is "
            "independent of the burden and the invariant quantity is the "
            "absolute D, not the MOI. Handing Implementation an MOI to "
            "multiply by the patient's burden is therefore only valid above "
            "that crossover."
        ),
    }


# --- task 7: Poisson zero class -------------------------------------------


def task7() -> dict:
    p = ref_params()
    sol = solve_course(p, n_out=N_OUT)
    z = sol.zero_class
    frac_above_1pct = float((z > 0.01).mean())

    # the explicit zero-class variant, and the MOI it needs
    moi_mf = bisect_min_log(
        lambda m: evaluate(solve_course(
            p.with_MOI(m), n_out=N_OUT)).passes, 1e-2, 1e3, tol=2e-3)
    moi_z = bisect_min_log(
        lambda m: evaluate(solve_course_z(
            p.with_MOI(m), n_out=N_OUT)).passes, 1e-2, 1e3, tol=2e-3)

    samples = []
    for moi in (0.5, 1.409, 3.2, 9.5):
        q = p.with_MOI(moi)
        s_mf = solve_course(q, n_out=N_OUT)
        s_z = solve_course_z(q, n_out=N_OUT)
        v_mf, v_z = evaluate(s_mf), evaluate(s_z)
        samples.append({
            "MOI": moi,
            "m0": float(s_mf.m[0]),
            "zero_class_at_t0": float(s_mf.zero_class[0]),
            "zero_class_max": float(np.nanmax(s_mf.zero_class)),
            "meanfield_delta_log": v_mf.delta_log,
            "meanfield_C_bar": v_mf.C_bar,
            "meanfield_passes": bool(v_mf.passes),
            "zeroclass_delta_log": v_z.delta_log,
            "zeroclass_C_bar": v_z.C_bar,
            "zeroclass_passes": bool(v_z.passes),
        })

    return {
        "zero_class_at_t0": float(z[0]),
        "zero_class_max": float(np.nanmax(z)),
        "fraction_of_course_above_1pct": frac_above_1pct,
        "mean_field_has_discrete_zero_class": False,
        "note": (
            "The mean-field model carries one compartment mean A_cell and has "
            "neither a discrete per-dose zero class nor a continuous lambda' "
            "conversion between treated and untreated pools; sunrise.zeroclass "
            "adds the discrete split."
        ),
        "moi_meanfield": moi_mf,
        "moi_with_zero_class": moi_z,
        "moi_inflation": (moi_z / moi_mf) if (moi_mf and moi_z) else None,
        "samples": samples,
    }


# --- task 8: cross-checks --------------------------------------------------


def eta_threshold() -> dict:
    """Reconcile Table 9.4.50's eta >= 0.7575 with the implemented Hill form.

    Ver.23 derived it by inverting psi = psi_max - (psi_max - psi_min)
    x^kappa / (1 + x^kappa) with x = eta*A/MIC_R.  The implemented form has
    (x^kappa - psi_min/psi_max) in the denominator, so the inverse differs.
    """
    p = STANDARD_PATIENT
    A = p.A_peak
    MIC = p.MIC_R0

    def psi_quoted(eta: float) -> float:
        x = (eta * A / MIC) ** p.kappa
        return p.psi_max - (p.psi_max - p.psi_min) * x / (1.0 + x)

    def psi_impl(eta: float) -> float:
        return psi(eta * A, MIC, p)

    def solve_zero(f) -> float:
        lo, hi = 1e-6, 10.0
        for _ in range(200):
            mid = math.sqrt(lo * hi)
            if f(mid) > 0.0:
                lo = mid
            else:
                hi = mid
        return math.sqrt(lo * hi)

    e_quoted = solve_zero(psi_quoted)
    e_impl = solve_zero(psi_impl)

    # the quoted value may have been derived before the MIC_R0 recalibration
    p_old = STANDARD_PATIENT.with_(MIC_R0=125.0)

    def psi_impl_old(eta: float) -> float:
        return psi(eta * A, 125.0, p_old)

    def psi_quoted_old(eta: float) -> float:
        x = (eta * A / 125.0) ** p.kappa
        return p.psi_max - (p.psi_max - p.psi_min) * x / (1.0 + x)

    return {
        "quoted_in_manuscript": 0.7575,
        "A_used_ug_per_mL": A,
        "MIC_R0_used": MIC,
        "eta_at_psi_zero_quoted_form": e_quoted,
        "eta_at_psi_zero_implemented_form": e_impl,
        "eta_quoted_form_at_MIC_R0_125": solve_zero(psi_quoted_old),
        "eta_implemented_form_at_MIC_R0_125": solve_zero(psi_impl_old),
        "eta_scales_as": "eta is proportional to MIC_R0 / A",
        "forms_agree": abs(e_quoted - e_impl) < 1e-3,
        "reconciles_with_0.7575": abs(e_quoted - 0.7575) < 5e-3,
        "note": (
            "The implemented denominator is x^kappa - psi_min/psi_max, which "
            "makes psi(A = MIC) = 0 exactly; the quoted derivation used "
            "1 + x^kappa, which does not. The two give different eta."
        ),
    }


def immune_replacement_invariance() -> dict:
    """Results (1)-(9) must not move when the immune term is swapped in."""
    p = ref_params()
    a = evaluate(solve_course(p, n_out=N_OUT))
    b = evaluate(solve_course(p.with_(
        induced_immunity=True, a_N=0.0, d_N=0.0), n_out=N_OUT))
    return {
        "fixed_delta_log": a.delta_log,
        "swapped_at_rest_delta_log": b.delta_log,
        "fixed_C_bar": a.C_bar,
        "swapped_at_rest_C_bar": b.C_bar,
        "unchanged": bool(abs(a.delta_log - b.delta_log) < 1e-6
                          and abs(a.C_bar - b.C_bar) / a.C_bar < 1e-6),
    }


def induced_immunity_contribution() -> dict:
    """Results (10): what the induced term buys at the standard patient."""
    p = ref_params()
    out = {}
    for name, q in (
        ("fixed_at_rest", p),
        ("induced", p.with_(induced_immunity=True)),
        ("induced_untreated", p.with_(
            induced_immunity=True, sunrise_on=False, antibiotic_on=False)),
    ):
        sol = solve_course(q, n_out=N_OUT)
        v = evaluate(sol)
        out[name] = {
            "delta_log": v.delta_log,
            "C_bar": v.C_bar,
            "passes": bool(v.passes),
            "N_final": float(sol.N_eff[-1]),
            "N_PMN_final_per_mL": float(sol.N_eff[-1] * q.N_PMN_rest),
            "N_max_over_course": float(sol.N_eff.max()),
            "N_PMN_max_per_mL": float(sol.N_eff.max() * q.N_PMN_rest),
        }
    out["observed_N_PMN_standard_patient"] = STANDARD_PATIENT.N_PMN
    out["reproduces_observed_pyuria"] = bool(
        out["induced_untreated"]["N_PMN_max_per_mL"]
        >= 0.5 * STANDARD_PATIENT.N_PMN
    )
    out["note"] = (
        "The standard patient's 900 WBC/uL (9e5 /mL) is the value observed at "
        "presentation. The induced model is driven by G = L/(K_G + L), which "
        "saturates at the initial burden, so N_PMN_max is the figure to "
        "compare against the clinical observation."
    )
    return out


def calibration_structure_check() -> dict:
    """The grid.py-class bug: a calibrated value applied in workers only.

    Params is a frozen dataclass passed explicitly into every worker, and no
    module-level mutable holds a calibrated value, so the parent/worker split
    cannot diverge.  The check below proves it for the quantity that bug hit.
    """
    p = ref_params()
    direct = evaluate(solve_course(p, n_out=N_OUT)).delta_log
    with Pool(2) as pool:
        viaworkers = pool.map(_dlog_worker, [p, p])
    return {
        "delta_log_parent": direct,
        "delta_log_workers": viaworkers,
        "identical": bool(all(abs(v - direct) < 1e-12 for v in viaworkers)),
        "module_level_mutable_state": False,
    }


def _dlog_worker(p: Params) -> float:
    return evaluate(solve_course(p, n_out=N_OUT)).delta_log


def task8() -> dict:
    return {
        "eta_reconciliation": eta_threshold(),
        "immune_replacement_invariance": immune_replacement_invariance(),
        "induced_immunity_contribution": induced_immunity_contribution(),
        "calibration_structure": calibration_structure_check(),
        "analytic_table": analytic.analytic_table(ref_params()),
    }


def dump(name: str, obj) -> None:
    (RESULTS / f"{name}.json").write_text(
        json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"[written] results/{name}.json", flush=True)


def main() -> None:
    with Pool(4) as pool:
        t4 = task4(pool)
        dump("task4", t4)
        t5 = task5(pool)
        dump("task5", t5)
        t6 = task6(pool)
        dump("task6", t6)
    t7 = task7()
    dump("task7", t7)
    t8 = task8()
    dump("task8", t8)

    with (RESULTS / "task4_tolerance.csv").open(
            "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(t4["axes"][0].keys()))
        w.writeheader()
        w.writerows(t4["axes"])

    print("=== TASK 4 tolerance limits (reference design) ===")
    for r in t4["axes"]:
        lo = "-" if r["lowest_passing"] is None else f"{r['lowest_passing']:.4g}"
        hi = "-" if r["highest_passing"] is None else f"{r['highest_passing']:.4g}"
        print(f"  {r['axis']:<34} pass in [{lo:>10} , {hi:>10}]")
    if t4["enzymology"]:
        e = t4["enzymology"]
        print(f"  k_TPD0_crit = {e['k_TPD0_crit_per_min']:.5f} /min "
              f"({e['fold_of_baseline']:.3f}x baseline) -> "
              f"k_cat >= {e['k_cat_crit_per_h']:.1f} /h, "
              f"ClpXP availability >= {e['ClpXP_available_fraction_crit']:.3f}, "
              f"K_D <= {e['K_D_crit_uM']:.4f} uM")

    print("\n=== TASK 5 IBC band ===")
    for r in t5["rows"]:
        print(f"  psi_max_in={r['psi_max_in']:.2f}  Cbar={r['C_bar']:8.2f}  "
              f"dlog={r['delta_log']:5.2f}  pass={r['passes']}  "
              f"Cbar/K_in={r['C_bar_over_K_in']:.3f}  "
              f"K_in_max={r['K_in_max_passing']:.1f}  "
              f"residual={r['residual_CFU_bladder_72h']:.3g} CFU")

    print("\n=== TASK 6 burden dependence ===")
    for r in t6["rows"]:
        print(f"  C0={r['C0_ex']:.1e} CFU/mL  MOI={r['moi_optimal']:.4g}  "
              f"D/dose={r.get('D_per_dose', float('nan')):.3e}")
    print("  monotone increasing with burden:",
          t6["monotone_increasing_with_burden"])

    print("\n=== TASK 7 Poisson zero class ===")
    print(f"  zero class at t0 = {t7['zero_class_at_t0']:.3%}, "
          f"max over course = {t7['zero_class_max']:.3%}, "
          f"fraction of course above 1% = "
          f"{t7['fraction_of_course_above_1pct']:.1%}")
    print(f"  MOI mean-field {t7['moi_meanfield']:.4g} -> with zero class "
          f"{t7['moi_with_zero_class']:.4g} "
          f"(x{t7['moi_inflation']:.3f})")
    for s in t7["samples"]:
        print(f"    MOI={s['MOI']:6.3f} m0={s['m0']:7.2f} "
              f"e^-phi*m={s['zero_class_at_t0']:.3%} "
              f"mf_pass={s['meanfield_passes']} z_pass={s['zeroclass_passes']}")

    print("\n=== TASK 8 cross-checks ===")
    e = t8["eta_reconciliation"]
    print(f"  eta: quoted form {e['eta_at_psi_zero_quoted_form']:.4f}, "
          f"implemented {e['eta_at_psi_zero_implemented_form']:.4f}, "
          f"manuscript 0.7575 -> reconciles={e['reconciles_with_0.7575']}, "
          f"forms agree={e['forms_agree']}")
    inv = t8["immune_replacement_invariance"]
    print(f"  immune replacement unchanged: {inv['unchanged']}")
    cal = t8["calibration_structure"]
    print(f"  parent/worker identical: {cal['identical']}")
    ind = t8["induced_immunity_contribution"]
    for k, v in ind.items():
        print(f"  {k:<22} dlog={v['delta_log']:5.2f} Cbar={v['C_bar']:8.2f} "
              f"pass={v['passes']} N_PMN_final={v['N_PMN_final_per_mL']:.2e}")


if __name__ == "__main__":
    main()
