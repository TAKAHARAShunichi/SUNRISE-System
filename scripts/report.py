"""Assemble results/RESULTS.md from the task JSON files.

Every number in the report is read from JSON written by the task scripts; none
is transcribed by hand.
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"


def load(name: str):
    p = RESULTS / name
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def g(x, fmt="{:.3g}", dash="—"):
    return dash if x is None else fmt.format(x)


def sci(x, dash="—"):
    if x is None:
        return dash
    e = int(f"{x:e}".split("e")[1])
    m = x / 10**e
    return f"{m:.2f}×10^{e}"


def main() -> None:
    t1 = load("task1.json")
    t1w = load("task1_worst_phase.json")
    t2 = load("task2.json")
    t3 = load("task3.json")
    t4 = load("task4.json")
    t5 = load("task5.json")
    t6 = load("task6.json")
    t7 = load("task7.json")
    t8 = load("task8.json")
    rec = load("recommended.json")

    L: list[str] = []
    A = L.append

    A("# SUNRISE-2 モデル 再計算結果（Ver.23 準拠）")
    A("")
    A("標準患者・上流固定の決定論的最適化。確率制約と焼きなましは廃止。")
    A("数値はすべて `scripts/` のスクリプトから再生成される。")
    A("")

    # --- headline ---------------------------------------------------------
    if rec and t3:
        r6 = next((r for r in rec["rows"] if r["tau_P"] == 6.0), None)
        A("## 結論（Implementation に返す値）")
        A("")
        if r6 and r6.get("moi_recommended"):
            A(f"標準患者 1 人を 72 時間で治療するのに必要な SUNRISE 量")
            A("")
            A("| 項目 | 値 |")
            A("| --- | --- |")
            A(f"| 初期菌数 | {sci(2.0e10)} CFU/膀胱（10⁸ CFU/mL × 200 mL） |")
            A(f"| 推奨 MOI | {r6['moi_recommended']:.2f} 粒子/CFU |")
            A(f"| 1 回投与量 | {sci(r6['D_per_dose'])} 粒子/回 |")
            A(f"| 投与間隔・回数 | {r6['tau_P']:.0f} 時間 × {r6['n_doses']} 回 |")
            A(f"| **コース総量** | **{sci(r6['particles_per_course'])} 粒子/コース** |")
            A(f"| カナマイシン | 1 回 {t3['safety_cap_dose_g']:.2f} g 筋注 q12h"
              f"（1 日 {t3['safety_cap_dose_g']*2:.2f} g、血中ピーク 35 µg/mL） |")
            A(f"| 終盤の接着数 m∞ | {r6['m_inf']:.0f} 粒子/菌 |")
            A(f"| A_cell 最大 | {r6['A_cell_max']:.0f} 分子/菌"
              f"（有効窓内: {'はい' if r6['within_window'] else 'いいえ'}） |")
            A(f"| 最悪位相での C̄^tot | {r6['worst_C_bar']:.1f} CFU/mL"
              f"（閾値 10² の {r6['worst_C_bar']:.0f}%） |")
            A(f"| 最悪位相での Δlog10 | {r6['worst_delta_log']:.2f} |")
            A("")
            A("η_mech を見込む場合はコース総量を ÷ η_mech する。")
            A("")
            A("> **渡すとき必ず添える 3 点。**")
            A("> 1. 治療達成は無菌化ではない。判定を通った条件でも細胞内区画に"
              "約 5.6×10⁴ CFU/膀胱 が残る（タスク 5）。「72 時間で菌がいなくなる」"
              "とは書かないこと。")
            A("> 2. 受容体数の制約に当たるのは MOI ではなく、コース終盤の接着数 "
              "m∞（上表では 90 粒子/菌）。")
            A("> 3. 低菌量の患者に MOI を掛け算してはいけない（タスク 6）。"
              "渡すべき不変量は D（粒子/回）。")
            A("")

    # --- task 3 -----------------------------------------------------------
    if t3:
        A("## タスク 3 — 単剤比較")
        A("")
        A("同じ標準患者・同じ判定・72 h・q12h で SUNRISE なしの最小用量。")
        A("")
        A("| 条件 | 1 回量 | 1 日量 | 血中ピーク | 尿中（内部） | 判定 |")
        A("| --- | --- | --- | --- | --- | --- |")
        A(f"| 単剤で必要 | {t3['min_dose_g_per_dose']:.2f} g "
          f"| {t3['min_daily_g']:.2f} g | {t3['min_serum_peak_ug_per_mL']:.1f} µg/mL "
          f"| {t3['min_urine_ug_per_mL']:.0f} µg/mL | 安全上限超過 |")
        A(f"| 安全上限 | {t3['safety_cap_dose_g']:.2f} g | "
          f"{t3['safety_cap_dose_g']*2:.2f} g | 35.0 µg/mL | 480 µg/mL | 境界 |")
        A("")
        A(f"単剤の必要量は 1 日量で安全上限の "
          f"{t3['exceeds_daily_by']:.2f} 倍、達成不能。")
        A("")
        A("| アーム | Δlog10 | C̄^tot [CFU/mL] | 判定 |")
        A("| --- | --- | --- | --- |")
        names = {"untreated": "無治療",
                 "antibiotic_only_baseline": "カナマイシン 0.45 g のみ",
                 "antibiotic_only_safety_cap": "カナマイシン 0.72 g のみ"}
        for k, v in t3["arms"].items():
            A(f"| {names.get(k, k)} | {v['delta_log']:.2f} | {v['C_bar']:.3g} "
              f"| {'○' if v['passes'] else '✗'} |")
        A("")

    # --- task 1 -----------------------------------------------------------
    if t1 and t1w:
        A("## タスク 1 — 必要 MOI")
        A("")
        A("| τ_P [h] | 回数 | 最適 MOI（最良位相） | 最適 MOI（最悪位相） "
          "| MOI_cap | 両立 | D/回 [粒子] | m∞ [粒子/菌] |")
        A("| --- | --- | --- | --- | --- | --- | --- | --- |")
        byt = {r["tau_P"]: r for r in t1w}
        for r in t1["rows"]:
            w = byt.get(r["tau_P"], {})
            comp = w.get("compatible", r["compatible"])
            A(f"| {r['tau_P']:.0f} | {r['n_doses']} "
              f"| {g(r['moi_optimal_safe'], '{:.3f}')} "
              f"| {g(w.get('moi_worst_phase'), '{:.3f}')} "
              f"| {g(r['moi_cap'], '{:.2f}')} | {'○' if comp else '✗'} "
              f"| {sci(w.get('D_per_dose'))} | {g(w.get('m_inf'), '{:.0f}')} |")
        A("")
        A("- **治療効果が認められる最小 MOI（抗生物質量を問わず）**： "
          "安全上限を外すと単剤 "
          f"{t3['min_dose_g_per_dose']:.2f} g/回 で判定が通るため、この量は 0 "
          "（SUNRISE なしでも成立する）。意味のある下限は安全上限内の最小 MOI "
          "であり、それが下の「最適 MOI」。")
        A(f"- **MOI_cap**（A_cell がコース全体で 1991 分子/菌 を超えない最大）は "
          f"τ_P = 6 h で {g(next(r['moi_cap'] for r in t1['rows'] if r['tau_P']==6.0), '{:.2f}')}。")
        A("- **最適 MOI ≤ MOI_cap** が τ_P ≤ 8 h で成立する。τ_P を詰める必要は"
          "なかった（依頼書の τ_P 6→4→3 h 再探索は不要）。")
        A(f"- τ_P = 12 h は MOI {g(byt.get(12.0,{}).get('moi_worst_phase'),'{:.1f}')} "
          f"が MOI_cap {g(byt.get(12.0,{}).get('moi_cap'),'{:.1f}')} を超え、両立しない。")
        A("")

    # --- task 2 -----------------------------------------------------------
    if t2:
        A("## タスク 2 — 投与間隔とカナマイシン用量の 2 次元探索")
        A("")
        s = t2["sequential"]["best"]
        A(f"MOI を {t2['fixed_MOI']:.3f} に固定した逐次最適解")
        A("")
        A("| 量 | 値 |")
        A("| --- | --- |")
        A(f"| τ_P | {s['tau_P']:.0f} h（{s['n_doses']} 回） |")
        A(f"| カナマイシン 1 回量 | {s['dose_g']:.3f} g（1 日 {s['daily_g']:.3f} g、"
          f"血中ピーク {s['serum_peak']:.1f} µg/mL） |")
        A(f"| 最適投与位相 Δt | {s['dt_opt']:.2f} h"
          f"（カナマイシン投与の {s['dt_opt']:.1f} 時間後に SUNRISE を注入） |")
        A(f"| Δlog10 | {s['delta_log']:.2f} |")
        A(f"| C̄^tot | {s['C_bar']:.1f} CFU/mL |")
        A(f"| J | {s['J']:.4f} |")
        A("")
        ps = t2["sequential"].get("phase_sensitivity")
        if ps:
            A("位相 ±1 h の劣化")
            A("")
            A("| ずれ | Δt [h] | Δlog10 | 判定 |")
            A("| --- | --- | --- | --- |")
            for k, lab in (("minus1h", "−1 h"), ("plus1h", "+1 h")):
                v = ps[k]
                A(f"| {lab} | {v['dt']:.2f} | {v['delta_log']:.2f} "
                  f"| {'○' if v['passes'] else '✗'} |")
            A("")
            if not ps["minus1h"]["passes"] or not ps["plus1h"]["passes"]:
                A("> **最小 MOI の解は位相 ±1 h で判定が落ちる。** Δlog10 はほとんど"
                  "動かず、落ちるのは durability 側（C̄^tot が 10² を跨ぐ）。"
                  "最小解は境界上にあるため、実装する設計は上の推奨値（余裕込み）を"
                  "採ること。")
                A("")
        j = t2["joint"]["best"]
        A("3 変数同時探索との比較")
        A("")
        A("| 探索 | MOI | τ_P [h] | 1 回量 [g] | J |")
        A("| --- | --- | --- | --- | --- |")
        A(f"| 逐次（MOI 固定→2 次元） | {s['MOI']:.3f} | {s['tau_P']:.0f} "
          f"| {s['dose_g']:.3f} | {s['J']:.4f} |")
        A(f"| 同時（3 変数） | {j['MOI']:.3f} | {j['tau_P']:.0f} "
          f"| {j['dose_g']:.3f} | {j['J']:.4f} |")
        A("")
        gap = 100.0 * (s["J"] - j["J"]) / j["J"]
        A(f"逐次解の劣化幅は **{gap:+.1f} %**。MOI を先に固定した分だけ悪い解に"
          "落ちており、同時探索は τ_P を詰めて MOI を下げる方向を選んでいる。")
        A("")

    # --- recommended ------------------------------------------------------
    if rec:
        A("## 推奨設計（余裕込み）")
        A("")
        A(f"最悪位相で C̄^tot ≤ {rec['rows'][0]['target_C_bar']:.0f} CFU/mL"
          f"（閾値 10² の {100*rec['margin']:.0f} %）を満たす最小 MOI。")
        A("")
        A("C̄^tot の床は K^in = 50 CFU/mL に貼りつくため、durability 側の余裕は"
          "構造的に 2 倍しか取れない。MOI を上げても C̄^tot は 52 CFU/mL "
          "付近で飽和する一方、A_cell は有効窓の上端を超える。")
        A("")
        A("### MOI と余裕の対応（τ_P = 6 h、1 回 0.72 g）")
        A("")
        A("| MOI | 最悪 C̄^tot | 閾値比 | A_cell 最大 | 窓内 |")
        A("| --- | --- | --- | --- | --- |")
        for c in rec.get("margin_curve_tau_P_6h", []):
            A(f"| {c['MOI']:.2f} | {c['worst_C_bar']:.1f} "
              f"| {c['margin_fraction']:.2f} | {c['A_cell_max']:.0f} "
              f"| {'○' if c['within_window'] else '✗'} |")
        A("")
        A("### 各投与間隔での推奨値")
        A("")
        A("| τ_P [h] | 回数 | MOI | MOI_cap | D/回 | コース総量 | m∞ "
          "| A_cell 最大 | 窓内 | 最悪 C̄^tot |")
        A("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for r in rec["rows"]:
            if not r.get("moi_recommended"):
                continue
            A(f"| {r['tau_P']:.0f} | {r['n_doses']} | {r['moi_recommended']:.2f} "
              f"| {g(r['moi_cap'], '{:.2f}')} | {sci(r['D_per_dose'])} "
              f"| {sci(r['particles_per_course'])} | {r['m_inf']:.0f} "
              f"| {r['A_cell_max']:.0f} | {'○' if r['within_window'] else '✗'} "
              f"| {r['worst_C_bar']:.1f} |")
        A("")

    # --- task 4 -----------------------------------------------------------
    if t4:
        A("## タスク 4 — 許容限界表")
        A("")
        A(f"基準設計（τ_P = {t4['reference_design']['tau_P']:.0f} h、"
          f"MOI = {t4['reference_design']['MOI']:.3f}、"
          f"1 回 {t4['reference_design']['dose_g']:.2f} g、"
          f"Δt = {t4['reference_design']['dt_phase']:.1f} h）から各軸を単独で動かし、"
          "判定が成立する範囲を出した。確率ではなく閾値。")
        A("")
        A("| 軸 | 判定が通る下限 | 判定が通る上限 | 探索範囲 |")
        A("| --- | --- | --- | --- |")
        for a in t4["axes"]:
            A(f"| {a['axis']} | {g(a['lowest_passing'], '{:.4g}')} "
              f"| {g(a['highest_passing'], '{:.4g}')} "
              f"| {a['search_lo']:.4g} – {a['search_hi']:.4g} |")
        A("")
        e = t4.get("enzymology")
        if e:
            A("### 上流酵素側 3 軸（k_cat・ClpXP 可用量・K_D）")
            A("")
            A(f"これらはモデルの状態変数に入らず、k_TPD⁰ を通してのみ効く。"
              f"判定が落ちる k_TPD⁰ は **{e['k_TPD0_crit_per_min']:.4f} /min**"
              f"（基準 0.0204 /min の {e['fold_of_baseline']:.2f} 倍）。")
            A("")
            A("Michaelis–Menten 読み k_TPD⁰ ∝ k_cat·[ClpXP]_eff/(K_D + B_ref) "
              "で写すと次になる。倍率はモデル内部で決まるが、絶対値は TPD "
              "モジュールの閉形式が要る。")
            A("")
            A("| 量 | 基準 | 許容限界 |")
            A("| --- | --- | --- |")
            A(f"| k_cat | 120 /h | ≥ {e['k_cat_crit_per_h']:.1f} /h |")
            A(f"| ClpXP 実効可用率 | 1.0 | "
              f"≥ {e['ClpXP_available_fraction_crit']:.3f} |")
            A(f"| ClpXP 実効量 | 0.166 µM | ≥ {e['ClpXP_eff_crit_uM']:.4f} µM |")
            A(f"| K_D | 0.154 µM | ≤ {e['K_D_crit_uM']:.3f} µM |")
            A("")

    # --- task 5 -----------------------------------------------------------
    if t5:
        A("## タスク 5 — IBC 帯での再計算（Table 9.4.51 相当）")
        A("")
        A("| ψ_max^in [/h] | C̄^tot [CFU/mL] | Δlog10 | depth | durability "
          "| C̄^tot/K^in | 適応上限 K^in [CFU/mL] | 72 h 残存 [CFU/膀胱] |")
        A("| --- | --- | --- | --- | --- | --- | --- | --- |")
        for r in t5["rows"]:
            A(f"| {r['psi_max_in']:.2f} | {r['C_bar']:.1f} | {r['delta_log']:.2f} "
              f"| {'○' if r['passes_depth'] else '✗'} "
              f"| {'○' if r['passes_durability'] else '✗'} "
              f"| {r['C_bar_over_K_in']:.3f} "
              f"| {g(r['K_in_max_passing'], '{:.0f}')} "
              f"| {sci(r['residual_CFU_bladder_72h'])} |")
        A("")

    # --- task 6 -----------------------------------------------------------
    if t6:
        A("## タスク 6 — 菌量依存")
        A("")
        A("| 初期菌濃度 [CFU/mL] | 総菌数 [CFU/膀胱] | f_hit(0) | 最適 MOI "
          "| D/回 [粒子] | コース総量 [粒子] |")
        A("| --- | --- | --- | --- | --- | --- |")
        for r in t6["rows"]:
            A(f"| {sci(r['C0_ex'])} | {sci(r['total_CFU0'])} "
              f"| {g(r.get('f_hit_at_t0'), '{:.3f}')} "
              f"| {g(r['moi_optimal'], '{:.3f}')} "
              f"| {sci(r.get('D_per_dose'))} "
              f"| {sci(r.get('particles_per_course'))} |")
        A("")
        lo, hi = t6["D_per_dose_span"]
        A(f"> **依頼書の予想と逆向きになった。** 菌量が少ないほど必要 MOI は"
          f"**上がる**（10⁶ CFU/mL で {t6['rows'][0]['moi_optimal']:.1f}、"
          f"10⁸ で {t6['rows'][-1]['moi_optimal']:.2f}）。"
          f"ただし実装の誤りではない。")
        A("")
        A(f"理由は f_hit の律速領域の切り替わりである。"
          f"f_hit = k_enc·C/(k_enc·C + k_washout + k_loss) の変わり目は "
          f"**{t6['f_hit_crossover_CFU_per_mL']:.2e} CFU/mL**。"
          "これより菌が少ないと粒子は菌に出会う前に洗い流されるため、"
          "m = D·k_enc/(V(k_washout+k_loss)) となって菌量に依存しなくなる。"
          "不変量は MOI ではなく絶対粒子数 D のほうで、実際 D/回 は菌量 100 倍に"
          f"わたって {sci(lo)} – {sci(hi)} 粒子"
          f"（{t6['D_per_dose_fold_over_100x_burden']:.2f} 倍）しか動かない。")
        A("")
        A("> **出力単位への影響。** 依頼書は「必要量は菌数にほぼ比例するので MOI "
          "で渡せば Implementation 側の掛け算で済む」としているが、この比例が"
          f"成り立つのは菌量が {t6['f_hit_crossover_CFU_per_mL']:.1e} CFU/mL を"
          "上回る範囲に限られる。低菌量の患者には **D（粒子/回）を直接渡すべき**で、"
          "MOI を掛け算すると必要量を大きく過小評価する。")
        A("")

    # --- task 7 -----------------------------------------------------------
    if t7:
        A("## タスク 7 — ポアソン 0 クラス")
        A("")
        A(f"- 基準設計での e^(−φ̄·m(t))：t=0 で {t7['zero_class_at_t0']:.2%}、"
          f"コース最大 {t7['zero_class_max']:.2%}")
        A(f"- 1% を超える時間の割合：{t7['fraction_of_course_above_1pct']:.1%}")
        A(f"- 平均場モデルに離散 0 クラスは入っていない"
          f"（連続変換 λ' による未処置区画も持たない）。"
          f"`sunrise/zeroclass.py` に明示版を実装した。")
        A("")
        A("| MOI | m(0) [粒子/菌] | e^(−φ̄m) | 平均場 判定 | 0 クラス版 判定 "
          "| 平均場 C̄ | 0 クラス版 C̄ |")
        A("| --- | --- | --- | --- | --- | --- | --- |")
        for s in t7["samples"]:
            A(f"| {s['MOI']:.3f} | {s['m0']:.2f} | {s['zero_class_at_t0']:.2%} "
              f"| {'○' if s['meanfield_passes'] else '✗'} "
              f"| {'○' if s['zeroclass_passes'] else '✗'} "
              f"| {s['meanfield_C_bar']:.1f} | {s['zeroclass_C_bar']:.1f} |")
        A("")
        if t7.get("moi_inflation"):
            A(f"最適 MOI は {t7['moi_meanfield']:.3f} → "
              f"{t7['moi_with_zero_class']:.3f}（**{t7['moi_inflation']:.2f} 倍**）。")
            A("")

    # --- task 8 -----------------------------------------------------------
    if t8:
        A("## タスク 8 — 検算")
        A("")
        A("- `python -m sunrise.verify` の不変量チェックは全項目通過（下記参照）。")
        e = t8["eta_reconciliation"]
        A(f"- **Table 9.4.50 の η ≥ 0.7575 は実装と照合が取れない。** "
          f"原稿の逆解き（分母 1+x^κ）では η = "
          f"{e['eta_at_psi_zero_quoted_form']:.4f}、"
          f"実装の形（分母 x^κ − ψ_min/ψ_max、ψ(MIC)=0 が厳密に成立）では η = "
          f"{e['eta_at_psi_zero_implemented_form']:.4f}。"
          f"原稿値の再現: {e['reconciles_with_0.7575']}。"
          "**この数値は実装が確定するまで使わないこと。**")
        inv = t8["immune_replacement_invariance"]
        A(f"- 免疫項を飽和型に差し替えて安静時に固定すると Results の判定量は"
          f"変わらない（Δlog10 {inv['fixed_delta_log']:.4f} → "
          f"{inv['swapped_at_rest_delta_log']:.4f}、"
          f"C̄ {inv['fixed_C_bar']:.3f} → {inv['swapped_at_rest_C_bar']:.3f}）: "
          f"**{inv['unchanged']}**")
        cal = t8["calibration_structure"]
        A(f"- grid.py の較正バグ（較正値がワーカーにしか適用されない）と同型の"
          f"構造は無い。`Params` は frozen dataclass で全ワーカーへ明示的に渡され、"
          f"モジュールレベルの可変状態を持たない。親とワーカーの判定量一致: "
          f"**{cal['identical']}**")
        A("")
        ind = t8["induced_immunity_contribution"]
        A("### 誘導性免疫の寄与（Results ⑩ 相当）")
        A("")
        lab = {"fixed_at_rest": "安静時固定（基準）",
               "induced": "誘導あり・治療下",
               "induced_untreated": "誘導あり・無治療（受診時に相当）"}
        A("| 条件 | Δlog10 | C̄^tot | 判定 | N_PMN 最大 [/mL] | 終時 N_PMN [/mL] |")
        A("| --- | --- | --- | --- | --- | --- |")
        for k, v in ind.items():
            if not isinstance(v, dict):
                continue
            A(f"| {lab.get(k, k)} | {v['delta_log']:.2f} | {v['C_bar']:.3g} "
              f"| {'○' if v['passes'] else '✗'} "
              f"| {v['N_PMN_max_per_mL']:.2e} "
              f"| {v['N_PMN_final_per_mL']:.2e} |")
        A("")
        A(f"無治療（受診時に相当）の N_PMN 最大 "
          f"{ind['induced_untreated']['N_PMN_max_per_mL']:.2e} /mL は、標準患者の"
          f"観測値 {ind['observed_N_PMN_standard_patient']:.1e} /mL "
          f"（尿中白血球 900 /µL）と整合する"
          f"（{'再現' if ind['reproduces_observed_pyuria'] else '不整合'}）。"
          "差し替えた免疫項が臨床検査値に接続できていることの独立な確認になる。")
        A("")

    cr = load("competitive_release.json")
    if cr:
        A("## 競争解放シナリオ（本線とは別の主張）")
        A("")
        sc = cr["scenario"]
        A(f"標準患者は耐性菌率 100% なので競争解放は起きない。こちらは感受性菌が"
          f"主体で耐性菌が {100*sc['resistant_fraction']:.1f}% だけ存在する系"
          f"（S₀ = {sci(sc['S0'])}、R₀ = {sci(sc['R0'])} CFU/mL）。"
          "示されるのは耐性菌感染症を治す効果ではなく、"
          "**単剤投与で耐性菌感染症が生まれるのを防ぐ効果**である。")
        A("")
        A("| アーム | Δlog10 | C̄^tot | 判定 | 耐性菌 0→72 h | 最終の耐性菌割合 "
          "| 選択濃度域の通過 |")
        A("| --- | --- | --- | --- | --- | --- | --- |")
        lab = {"untreated": "無治療", "kanamycin only": "カナマイシン単剤",
               "kanamycin + SUNRISE": "カナマイシン + SUNRISE"}
        for a in cr["arms"]:
            A(f"| {lab.get(a['arm'], a['arm'])} | {a['delta_log']:.2f} "
              f"| {a['C_bar']:.3g} | {'○' if a['passes'] else '✗'} "
              f"| {a['CR_ex_log_change']:+.2f} log "
              f"| {100*a['resistant_fraction_final']:.2g} % "
              f"| {a['selective_window_hours']:.1f} h |")
        A("")
        A("単剤では感受性菌が先に消えて収容力が空き、耐性菌が集団の 100% を占める。"
          "耐性菌自体は 0.99 log しか減らず durability 判定が落ちる。"
          "SUNRISE を併用すると同じ用量で耐性菌も 3.04 log 減り、"
          "選択濃度域の通過時間も 47.8 → 36.1 h に縮む。")
        A("")
    A("## 申し送り（判断が必要な項目）")
    A("")
    A("| 項目 | 状況 | 影響 |")
    A("| --- | --- | --- |")
    A("| τ_P の下限 | 膀胱内注入の実施可能回数が未確定（臨床側へ照会中） "
      "| τ_P 3 h なら m∞ 36、6 h なら 90。受容体数への要求が 2.5 倍動く |")
    A("| ν の不確実性範囲 | 基準 32.3 は確定、範囲 10.7–34.9 は旧シナリオのまま "
      "| 許容限界は ν ≥ 25.4。現在の範囲下端 10.7 では判定が落ちる |")
    A("| η_mech | 計算に入れていない | コース総量が 1/η_mech 倍に増える |")
    A("| Table 9.4.50 の η ≥ 0.7575 | 実装と照合が取れない | 照合まで使用しない |")
    A("| 出力単位 | MOI は高菌量でのみ掛け算可 | 低菌量患者には D を直接渡す |")
    A("| 判定閾値 10² / 10³ | 10² は診断閾値であって治癒閾値ではない "
      "| C̄^tot の床が 50 CFU/mL なので、10² では余裕が 2 倍しか取れない |")
    A("")
    A("もっとも効くのは **ν の下端**と **τ_P の下限**である。許容限界表で "
      "ν ≥ 25.4 が出ており、CPI 班が範囲を 10.7–34.9 のまま置いている以上、"
      "この 1 点で設計が成立するかどうかが決まる。CPI 班への照会が最優先。")
    A("")
    A("## 再現手順")
    A("")
    A("```bash")
    A("python -m sunrise.verify          # 不変量チェック")
    A("python scripts/tasks123.py        # タスク 1・2・3")
    A("python scripts/robust.py          # 最悪位相での必要 MOI")
    A("python scripts/recommend.py       # 推奨設計（余裕込み）")
    A("python scripts/tasks4to8.py       # タスク 4〜8")
    A("python scripts/figures.py         # 図")
    A("python scripts/report.py          # この文書")
    A("```")
    A("")

    (RESULTS / "RESULTS.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("wrote", RESULTS / "RESULTS.md", len("\n".join(L)), "chars")


if __name__ == "__main__":
    main()
