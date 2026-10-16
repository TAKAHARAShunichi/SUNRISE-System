# Parameters: where every number came from

Generated from `sunrise/parameters.py`; that file is the single source of
truth and this document explains the entries that need more than one line.

Provenance tags are used strictly:

| Tag | Meaning |
|---|---|
| `measured` | measured by this project's wet lab |
| `literature` | a published measurement, cited below |
| `derived` | computed from other entries, with the derivation given |
| `assumed` | no measurement exists; treated as a sweep axis |

To print the full table:

```
python -m sunrise.parameters
```

## The parameters that decide the result

### `psi_in` — intracellular growth rate

Tagged `assumed` at 0.20 /h, and this is the single most consequential choice
in the model. The two things the intracellular compartment might represent
have growth rates three orders of magnitude apart.

Early intracellular bacterial communities are measured at a doubling time of
30–35 minutes over the first eight hours, by septum counting in time-lapse
imaging of explanted bladders — 1.19–1.39 /h
([PNAS 2004](https://doi.org/10.1073/pnas.0308125100)). An independent route
agrees: micropipette-purified communities at six hours post-infection hold a
median of about 10³ bacteria, which from a single founder implies 36 minutes
([Infection and Immunity 2018](https://doi.org/10.1128/IAI.00740-17)). That is
*faster* than growth in urine, so the intracellular niche is not growth
limiting in the acute phase.

Quiescent reservoirs are a different object: 4–12 bacteria, stable for months
([Microbiology Spectrum 2019](https://doi.org/10.1128/microbiolspec.BAI-0014-2019)),
which bounds sustained net growth below 0.002 /h.

Middle and later stages are bounded only — "more than an hour", so below
0.693 /h — and no number is published for the late or filamentation stages.

The current 0.20 /h corresponds to a 208-minute doubling time. It is not a
measured value for either compartment; it is a midpoint. The sign of `d`, and
therefore whether the model says eradication is possible, flips between the
two interpretations.

### `k_in` — intracellular carrying capacity

The identifiable quantity is the total intracellular burden per bladder, not
the product of a per-reservoir count and a reservoir count. Those two have
never been measured in the same animals, and multiplying the published
figures for each gives 4 × 10⁵, which exceeds the whole-bladder count measured
in the same study.

Direct measurements of the burden itself give 10⁴–10⁵ CFU per bladder at the
acute stage and about 10³ in the chronic phase with a sterile lumen
([Infection and Immunity 2011](https://doi.org/10.1128/IAI.05339-11)). The
dominant uncertainty is the intracellular fraction at six hours, reported
anywhere from 3% to 80% — a thirtyfold spread that propagates straight into
this parameter.

In a 200 mL bladder those become 50 CFU/mL (acute central), 15 (acute low) and
5 (chronic). The default is 50, swept 15–500.

A first-principles geometric check does not constrain this. The ceiling from
cell volume, available cytoplasm fraction, packing fraction and the measured
volume of an intracellular coccoid bacterium is about 1.5 × 10⁴ bacteria per
umbrella cell, and a human bladder has 1.4–3 × 10⁶ umbrella cells, so the
product ceiling is around 2 × 10¹⁰ — five orders of magnitude above the value
in use. The binding constraints are biological, not physical: 58% of
community-bearing cells exfoliate within one cycle
([eLife 2021](https://doi.org/10.7554/eLife.66481)), and the cycle itself runs
16–24 hours to filamentous release.

That same check does adjudicate one literature disagreement. The widely quoted
"up to 10⁶ bacteria per mature community" would need a sphere about 124 µm
across, 25–50 times the volume of its host cell, and the study that reports
that figure also measures roughly 160 CFU per community directly.

### `k_invade` and `k_emerge`

Neither has ever been measured as a rate.

Invasion is reported as a fraction at a timepoint, always normalised to
adherent or cell-associated bacteria rather than to the luminal pool. The
conversions: 4.4% of adherent bacteria invaded by one hour in mouse gives
0.045 /h; 0.5–10% of cell-associated over two hours in human 5637 cells gives
0.0025–0.053 /h
([EMBO Journal 2000](https://doi.org/10.1093/emboj/19.12.2803)); a
block-face electron microscopy census of mouse organoids at six hours gives
0.0065 /h for solitary cytoplasmic bacteria. All three denominators make these
upper bounds. The default is 0.01 /h.

Emergence is inferred from residence time. Mouse: invasion at about one hour,
flux at 16 hours, so 0.065 /h. Human bladder chip: first release at 15.6
hours, slow cases at 28 hours, giving 0.074 and 0.038 /h
([eLife 2021](https://doi.org/10.7554/eLife.66481)). The default is 0.065 /h.

Two cautions. Dispersal is a burst per community, not a continuous flux, and
it competes with exfoliation — of 100 communities tracked, 28 released
bacteria only, 27 exfoliated without releasing, 31 did both and 14 did
neither. A model that routes all intracellular bacteria back to the lumen
overstates re-seeding. And a constant rate collapses a sharply peaked
residence-time distribution into an exponential, which lets some bacteria
leave within an hour; if residence time matters to your conclusion, a two- to
four-stage Erlang chain is the right fix.

Equality between the two rates has no precedent in the published models and
is not used here.

### `k_exf` — intracellular clearance by exfoliation and autophagy

Zero, and no study reports this as a rate. The only derivable anchor is 0.043
/h, from 58 of 100 community-bearing cells exfoliating over a 20-hour imaging
window on a human bladder chip — which is a rate for community-bearing cells,
not for the umbrella layer, over a window the study does not define as an
exposure time.

Zero is both the conservative choice and arguably the best point estimate for
wild-type bacteria, because three independent findings converge. Autophagy
captures efficiently — over 80% in LAMP1-positive compartments at 24 hours —
but the compartments do not acidify, because the haemolysin reduces V-ATPase
recruitment. TRPML3 senses that failure and triggers lysosomal exocytosis, so
bacteria leave alive and membrane-encased
([Cell 2015](https://doi.org/10.1016/j.cell.2015.05.009)). And Atg16L1
deficiency is *protective* in vivo
([PNAS 2012](https://doi.org/10.1073/pnas.1203952109)).

Two structural consequences. Expulsion moves bacteria to the lumen rather than
killing them, so in a two-compartment model it is a transfer term and not a
death term; implementing it as death overstates total clearance. And
exfoliation is a discrete, inflammasome-driven burst around 24–48 hours that
is self-limiting once the umbrella layer is gone, so a constant rate is a
scenario axis rather than an uncertainty interval.

### `eta` — intracellular drug access

0.05, swept 0.005–0.1. The band is defensible; the functional form is not.

Accumulation pushes it up but slowly: whole-cell concentration ratios of
2.4–2.6 for kanamycin in rat embryo fibroblasts, reached over 2–4 days
([Antimicrob Agents Chemother 2001](https://doi.org/10.1128/AAC.45.11.2977-2986.2001)).
On a 24-hour urinary timescale that contribution is about 1, not 2.6.
Compartmentalisation pushes it down hard: the drug is largely confined to
lysosomes, which are 3–5% of cell volume, bounding the cytosolic ratio at
0.1–0.5 — and that is an assay sensitivity limit, not a measurement. Nobody
has measured drug concentration at the bacterium, in any system.

The reason the form is wrong: intracellularly the EC50 is roughly preserved
while the maximum effect collapses, and the Hill slope flattens about tenfold
([J Antimicrob Chemother 2007](https://doi.org/10.1093/jac/dkl489);
[Antimicrob Agents Chemother 2021](https://doi.org/10.1128/aac.01468-21)).
Low membrane potential does not block uptake or ribosome dissociation — it
turns killing bacteriostatic, and killing resumes on repolarisation
([eLife 2020](https://doi.org/10.7554/eLife.58706)). A concentration
multiplier can be overcome by raising the dose; a switch to bacteriostasis
cannot.

A higher-fidelity replacement splits the intracellular compartment into a
killable fraction and a refractory one, and applies the extracellular
pharmacodynamics to the killable fraction only.

Two unresolved items are worth stating. The mechanism literature conflicts:
the same tool (fluorescent gentamicin with a protonophore) has supported both
"depolarisation suppresses uptake" and "depolarisation does not block entry".
And the compartment itself is contested — one study describes the bacteria as
sequestered in late endosomes, where aminoglycosides concentrate 100–250×,
which would predict `eta > 1` rather than `eta ≪ 1`.

## Where this model and the standalone modules disagree

Running `python scripts/audit.py` reproduces this comparison.

| Quantity | Here | Standalone modules | Status |
|---|---|---|---|
| Maximum kill rate | 1.700 /h | 1.702 /h | agree |
| Resistance factor baseline | 0.3233 µM | 0.3233 µM | integrated |
| Degron degradation rate | 9.8185 × 10⁻⁴ /min | same | integrated |
| Extracellular capacity | 1 × 10⁹ CFU/mL | 6.29 × 10⁸ measured | the measured value is lower |
| Payload per particle | 32 molecules | 0.20 µM ≈ 120 molecules | two independent estimates support 30–35 |
| Injection success | 0.98 | 0.98 | integrated |
| Hill exponent | 2.95 | 1.722 after recalibration | see below |
| Resistant MIC | 125 µg/mL | 163 implied | see below |

The Hill exponent and the MIC are the substantive gap. The two formulations
differ in structure: this model shifts the MIC with the resistance factor
level and uses a four-parameter form in A/MIC, while the standalone model
routes the drug through a detox capacity and a Michaelis constant to an
intracellular active concentration. The exponents are therefore exponents on
different variables and cannot be substituted for one another.

The numerical coincidence is informative, though: 2.95 is exactly the
pre-recalibration exponent, so this model inherited the older fit. Refitting
this model's functional form to the recalibrated curve gives MIC 1.60 with
exponent 1.73 for the susceptible strain, and MIC 191 with exponent 3.05 for
the resistant strain. **The two strains need different exponents**, which a
single parameter cannot express — the detox mechanism changes the shape of the
curve, not just its position.

In the middle of the concentration range the two formulations disagree by up
to 2.2× in kill rate (worst around 100 µg/mL), converging above 300 µg/mL.
Since trough concentrations fall in exactly that middle range, this model kills
somewhat too fast at the trough.

Propagating the discrepancy does not change the conclusion, and sharpens it:
with the MIC derived from the recalibrated parameters, the minimum peak
concentration for the antibiotic alone rises from 516 to 674–795 µg/mL — further
above the clinical ceiling — while the combination moves from 92 to 118–138
µg/mL, still under a third of the ceiling.

## Present in the standalone modules, absent here

- The protease competition window after a kanamycin dose.
- The structurally non-injectable subpopulation, which the standalone
  integration identifies as its largest unresolved risk. A rough estimate
  suggests a 72-hour course demands a fraction below about 3 × 10⁻¹⁰, three
  orders of magnitude tighter than the 18-hour in-vitro figure.
- The detox capacity and Michaelis constant route to intracellular active
  concentration, replaced here by the MIC shift.
- The closed form from degrader turnover to effective degradation rate,
  replaced here by a linear map with a calibrated constant.

## Endpoints

The 10² CFU/mL threshold is a **diagnostic** criterion for symptomatic
cystitis with pyuria
([N Engl J Med 2013](https://doi.org/10.1056/NEJMoa1302186)), not a cure
criterion; no guideline or regulatory document retrieved uses it as one. The
regulatory analogue is below 10³ CFU/mL of the baseline uropathogen at a
post-therapy test-of-cure visit, paired with symptom resolution. Both
thresholds are available in `parameters.py`; running the model at each and
reporting the sensitivity is the defensible approach.

Judging at the end of the course, with drug still present, measures
suppression rather than cure, and is the opposite of the regulatory design.
Averaging over a window has no precedent at all; it is defensible as a way to
remove the dependence on dosing phase, not as a clinical mimic.

The log-reduction endpoint has better support: stasis, 1-log, 2-log and 3-log
kill is the standard pharmacodynamic ladder, and the argument that host
clearance finishes the job below some density is itself expressed in
log-reduction terms. Note that the model has no mechanism for that handoff —
its immune term is linear and does not strengthen at low density — so using
the argument means importing a claim the model does not make. Adding a
saturating immune term does not fix this: that term acts on the extracellular
compartment, which does not enter `d`.
