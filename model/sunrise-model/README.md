# SUNRISE within-host model

A within-host model of a urinary tract infection treated with kanamycin
together with a delivered protein degrader that destroys the resistance factor
inside the bacterium.

Four bacterial compartments — susceptible or resistant, extracellular or
intracellular — plus the drug, the delivery particles, and the degrader load
inside each bacterium, over a treatment course in the bladder.

```python
from sunrise import metrics

metrics.evaluate()                 # the baseline course
metrics.evaluate(dose=0.0)         # antibiotic alone
metrics.compare_regimens()         # all four arms side by side
```

## Why you might want this

Most of this repository is not the ODE. It is the bookkeeping around it:
where every parameter came from, which ones are guesses, and what the model
structurally cannot tell you regardless of how you tune it. Three pieces are
worth borrowing even if your project has nothing to do with urinary tract
infections.

**Parameters carry their provenance.** Every entry in `parameters.py` is
tagged `measured`, `literature`, `derived` or `assumed`, with the measurement
conditions written out and a plausible range attached. `parameters.table()`
prints the whole thing as markdown for a wiki page. A test asserts that no
parameter can be added without a source and a unit.

**The structural analysis runs before the dose sweep.** `analysis.py` answers
"can this treatment work at all" from the parameters alone, in closed form,
without integrating anything. For us the answer was no — and it would have
stayed hidden under a dose sweep that looked like it was converging. See
*What the analysis found*, below.

**Two endpoints, because one hides the failure.** A course can drive the
population a long way down and then lose it. Reporting only the final density,
or only the log reduction, makes those two look the same. `metrics.py` reports
both and a verdict on each.

## Install and check

```
pip install -r requirements.txt
python scripts/validate.py       # reproduces the published figures
python tests/test_model.py       # 19 structural and consistency checks
```

`validate.py` fixes the one quantity with no measurement behind it — the
degradation rate — against a single published result, then checks ten other
published numbers that were not used in that fit. All ten agree to within 5%.

## What the analysis found

In the low-density limit the resistant population is a linear two-compartment
system, so whether the infection can be eradicated is an eigenvalue question
on two numbers:

```
a_bar = psi_max − kill(A) − k_imm − k_wash − k_invade      (extracellular)
d     = psi_in  − kill(eta·A) − k_exf − k_emerge           (intracellular)
```

Both eigenvalues are negative only when `a_bar + d < 0` **and**
`a_bar·d > k_emerge·k_invade`. Under any useful regimen `a_bar < 0`, so the
second condition forces `d < 0`. **If the intracellular compartment grows on
its own, no amount of extracellular killing can drive the system to
extinction** — the determinant is negative and one eigenvalue stays positive
whatever the dose.

That has three consequences we would not have found by sweeping doses.

Strengthening the immune term does nothing. We replaced the linear immune
clearance with a saturating form that rises as density falls, calibrated so
the untreated dynamics were unchanged. Raising the low-density clearance rate
to 5,000 /h left the dominant eigenvalue unmoved in the fourth decimal, because
the term does not appear in `d`.

Shrinking the reservoir does not help either. The intracellular carrying
capacity is a capacity, not a rate: the growth exponent is identical across
three orders of magnitude in `k_in`. What changes is where the population
settles, not whether it comes back.

And that settling point is what a durability endpoint measures. When `d > 0`
the system does not grow without bound — it fills the intracellular
compartment and stops. The final-cycle total is then 0.95 `k_in` across three
orders of magnitude, so judging the course on total density is, arithmetically,
judging `k_in`. That is why the endpoint moved to a log reduction at the nadir,
with the density kept as a second, separately reported criterion.

```python
from sunrise import analysis
analysis.decay_condition()   # a_bar, d, the eigenvalue, and the verdict
analysis.floor_level()       # where it settles if it cannot be eradicated
analysis.levers()            # what value each lever would need to flip d
```

## Adapting it to your project

Rewrite `parameters.py` and nothing else. `model.py` reads every number from
there; there are no literals in the equations. If your system has an
intracellular or otherwise protected compartment — a biofilm, a phagosome, a
persister state — the decay condition applies to it unchanged, with `d` built
from whatever terms govern that compartment.

Things worth checking in your own model, which cost us several revisions:

- Is your endpoint sensitive to the dose, or is it pinned by a capacity? Sweep
  the capacity over two orders of magnitude and see whether the verdict tracks
  it.
- Is your growth exponent really the asymptotic one? Seed the simulation at a
  density small relative to **every** capacity, not just the obvious one. We
  twice measured an exponent that turned out to be an artefact of the seed.
- Does your calibration point move when you change an unrelated parameter?
  Ours did, until the calibration was pinned to the conditions its target came
  from.

## Layout

```
sunrise/
  parameters.py    every parameter, with provenance, units and ranges
  model.py         the ODE system, drug profile, and closed forms
  metrics.py       endpoints and verdicts
  analysis.py      eigenvalue condition, floor level, levers
scripts/
  validate.py      reproduce the published figures
  figures.py       regenerate the figures
  audit.py         compare against the standalone module parameters
docs/
  parameters.md    the provenance table and the open questions
tests/
  test_model.py    structural and consistency checks
```

## Caveats we would want a reader to know

The degradation rate is calibrated, not measured, and the calibration is a
single point. The closed form from the standalone degradation module gives
0.041–0.091 /min over the relevant load range against a placeholder of 0.14,
so the placeholder is optimistic by roughly 1.5×.

The intracellular access fraction `eta` is a lumped approximation that the
literature does not support as a functional form. Intracellularly the EC50 is
roughly preserved while the maximum effect collapses and the Hill slope
flattens about tenfold, so a single multiplier on concentration will
mis-extrapolate across a wide sweep — which is exactly how it is used here.

The intracellular compartment is a single pool, and the two things it might
represent behave oppositely. Acute intracellular communities grow at a measured
1.19–1.39 /h; quiescent reservoirs are below 0.002 /h. The first gives `d > 0`
and no eradication; the second gives `d < 0`. A single pool must
mis-specify one of them.

The clinical concentration ceiling is derived from systemic dosing. For
intravesical instillation it needs re-deriving from the fraction absorbed, and
the published instillation literature suggests a substantially higher ceiling.

## Licence

Code is MIT. Documentation and figures are CC BY 4.0. Parameter values taken
from the literature belong to their original authors and are cited in
`docs/parameters.md`.
