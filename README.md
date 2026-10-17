# SUNRISE-System

2026ウィキを最後に作るときのリポジトリ。

## SUNRISE-2 モデル（局所投与時の感染消退）

`SUNRISE-2_Ver23.md` の数式仕様から再実装したモデルと、依頼書
（[モデル 申し送り 2026-10-17](https://claude.ai/code/artifact/98de8643-6f7e-4280-9700-30e0703df57d)）
のタスク 1〜8 の計算結果。

**結果は [`results/RESULTS.md`](results/RESULTS.md) にまとまっている。**

### 構成

| パス | 内容 |
| --- | --- |
| `sunrise/params.py` | 全パラメータ（Ver.23 Table 9.4.5 / 7 / 10 / 14 / 22）と用量換算 |
| `sunrise/model.py` | 9 変数 ODE（CPI 体内動態・薬力学・自然免疫・集団動態） |
| `sunrise/metrics.py` | 判定量 Δlog10 と C̄^tot、有効窓占有率 |
| `sunrise/analytic.py` | 積分なしで決まる閉形式（Results ① 相当） |
| `sunrise/design.py` | 決定論的設計探索（二分法・目的関数 J） |
| `sunrise/zeroclass.py` | ポアソン 0 クラスを明示した変種（タスク 7） |
| `sunrise/verify.py` | 不変量チェック 25 項目 |
| `scripts/` | タスク実行スクリプトと図・報告書の生成 |
| `results/` | 数値（JSON / CSV）と `RESULTS.md` |
| `figures/` | 図 8 枚 |

### 実行

```bash
pip install numpy scipy matplotlib

python -m sunrise.verify          # 不変量チェック（25/25）
python scripts/tasks123.py        # タスク 1・2・3        （約 4 分）
python scripts/robust.py          # 最悪位相での必要 MOI  （約 1 分）
python scripts/recommend.py       # 推奨設計（余裕込み）  （約 3 分）
python scripts/tasks4to8.py       # タスク 4〜8           （約 6 分）
python scripts/figures.py         # 図
python scripts/report.py          # results/RESULTS.md
```

数値はすべてスクリプトから再生成される。`RESULTS.md` の表は
`results/*.json` から組み立てており、手写しは無い。

### Ver.22 からの差し替え

依頼書どおり、期待送達量を ν = 32.3 分子/接着粒子（貫通確率 φ̄ = 0.7697 は
λ' と 0 クラスにのみ使用）、ψ_max^in = 1.3 /h（早期 IBC）、初期条件を耐性菌率
100%（S₀ = 0、R₀ = 10⁸ CFU/mL）、抗生物質の制約を用量側（1 回 0.72 g、
1 日 1.43 g、血中ピーク 35 µg/mL）に置き換えた。免疫項は
k_imm = ε_urine·k_PMN·N_PMN/(1 + C^ex/K_cap) に差し替え、安静時
N_PMN = 10⁵ /mL で従来の 0.10 /h を再現することを不変量に入れている。

確率制約（成立確率 90%）と焼きなましは廃止した。
