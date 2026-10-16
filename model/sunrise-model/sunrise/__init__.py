"""SUNRISE: a within-host model of targeted degradation of a resistance factor.

Four bacterial compartments (susceptible / resistant, extracellular /
intracellular) under combined kanamycin and a delivered protein degrader,
over a treatment course in the bladder.

Quick start::

    from sunrise import parameters as prm, metrics

    metrics.evaluate()                        # the baseline course
    metrics.evaluate(dose=0.0)                # antibiotic alone
    metrics.compare_regimens()                # all four arms

Before reading a dose sweep, check what the model can do at all::

    from sunrise import analysis
    analysis.decay_condition()                # can it be eradicated?
    analysis.floor_level()                    # if not, where does it settle?

Adapting this to another organism, drug or resistance factor should only
require rewriting ``parameters.py``.
"""

from . import analysis, metrics, model, parameters

__all__ = ["analysis", "metrics", "model", "parameters"]
__version__ = "1.0.0"
