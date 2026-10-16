# Reusing this for your own system

Written for another iGEM team arriving with a different organism, drug or
mechanism. The equations here are specific to our project; the three habits
below are not, and they are what we would actually recommend copying.

## Start by asking what your model cannot do

Before any dose sweep, work out whether the thing you are modelling can reach
the outcome you want **at all**. For a two-compartment system this is an
eigenvalue question and takes a page of algebra. Ours:

```
a_bar = (growth) − (kill) − (clearance) − (entry into the protected pool)
d     = (growth inside) − (kill inside) − (death inside) − (exit)
```

with eradication possible only when `a_bar + d < 0` and
`a_bar·d > (exit)·(entry)`. Since any useful regimen makes `a_bar` negative,
the second condition forces `d < 0`, and no amount of work on the exposed
compartment can substitute.

We found this after weeks of dose sweeps that looked like they were
converging on a feasible region. The sweeps were real; they just answered a
question we had not checked was the right one. `analysis.decay_condition()`
is thirty lines and would have told us on day one.

The generalisation: if your system has a compartment your intervention cannot
reach — a biofilm, a phagosome, a persister state, a tissue reservoir — write
down its net rate as a single scalar and check the sign before anything else.

## Make your endpoint prove it is measuring the regimen

Our first endpoint was the bacterial density at the end of the course. It
looked sensible and it was almost meaningless: whenever the protected
compartment grows, the system settles at a level set by that compartment's
carrying capacity, and the endpoint becomes `0.95 × capacity` no matter what
dose you give. Sweeping the dose a thousandfold moved it by 8%.

A cheap test, worth running on any endpoint:

```python
for k in (0.01, 0.1, 1.0, 10.0):          # scale a capacity, not a rate
    print(k, my_endpoint(capacity=k * baseline_capacity))
```

If your endpoint tracks the capacity proportionally, it is reporting that
parameter rather than your intervention.

What we changed to: report two criteria, because a regimen can pass one and
fail the other, and the gap between them is the interesting part.

- **Depth** — log reduction at the nadir. This is the form the
  pharmacodynamic literature uses, and it does not depend on where you
  truncate the course.
- **Durability** — density over the final dosing cycle. This is closer to
  what regulators assess, though they do it after therapy ends.

Our antibiotic-alone arm passes depth (3.3 log) and fails durability
(rebounds to 7 × 10⁷). Reporting either number alone hides that completely.

## Tag every parameter with its provenance, and let the code enforce it

`parameters.py` gives every entry a `source` of `measured`, `literature`,
`derived` or `assumed`, the measurement conditions, and a plausible range. A
test fails if anything lacks a source or a unit, and `parameters.table()`
prints the whole thing as markdown for a wiki page.

This is dull and it repeatedly earned its keep. Two examples from our own
repository. A rate we had cited to a published model turned out not to exist
in that paper — it has no compartment of that kind at all, and the value had
propagated through revisions as a citation nobody re-opened. And a parameter
we had treated as a mild assumption turned out to have measurements spanning
two orders of magnitude, with our value sitting at the extreme optimistic end.

Writing the conditions out is what catches these. "0.03 /h, from
[a paper]" survives review; "0.03 /h, extracellular to intracellular, cited to
[a paper] — but that model has no intracellular compartment" does not.

## Smaller things that cost us revisions

**Seed the low-density limit relative to every capacity.** We twice measured
a growth exponent that turned out to be an artefact: the initial density was
small compared to the obvious capacity but not the other one. Scaling the
seed to the capacity under test made the exponent constant across three
orders of magnitude, as it should have been.

**Pin your calibration to the conditions its target came from.** We calibrate
one free constant against one published result. When we later changed an
unrelated parameter, the calibration silently moved with it and every
downstream number shifted. The fix is three lines: save, override, restore.

**Check that a parallel worker and the parent agree.** A figure annotation of
ours was computed in the parent process while the grid behind it was computed
in workers that had received a calibration the parent never got. The figure
was internally inconsistent for several revisions. If you use
`multiprocessing`, pass state through `initargs` rather than relying on
`fork` inheritance.

**Write a test for each claim you intend to publish.** Ours include "immune
clearance cannot make eradication possible", "the exponent does not depend on
the capacity" and "depth and durability can disagree". One of them failed
when written, which is how we learnt that the antibiotic concentration is not
a purely extracellular lever — it reaches the protected compartment through
the access fraction, and does enter `d`, just at an impractical dose. That
correction came from the test, not from reading the equations again.

## If you want to run ours

```
pip install -r requirements.txt
python scripts/validate.py       # reproduces our published figures
python tests/test_model.py
python -m sunrise.parameters     # the provenance table
```

Then rewrite `parameters.py`. `model.py` reads every number from there, so
there are no literals in the equations to hunt down. If your mechanism
differs structurally — a different number of compartments, a different
resistance mechanism — `model.rhs` is about sixty lines and is meant to be
edited rather than configured.
