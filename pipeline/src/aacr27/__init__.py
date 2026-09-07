"""aacr27 — pan-TCGA audit of H&E-inferred tumour immune-microenvironment signatures.

Research question
-----------------
Do pathology models that appear to predict tumour immune-microenvironment
signatures from H&E actually carry immune-specific information, or are they
reading tumour purity, tissue source site, and cancer type?

Two controls carry the paper, and neither exists in the current image-to-
expression literature:

1. The **Venet null** — size-matched random gene sets pushed through the
   identical scoring and prediction pipeline. Venet et al. (PLoS Comput Biol
   2011) showed >90% of random gene sets predict breast cancer outcome; nobody
   has asked the equivalent question of image-predicted signatures.

2. The **purity decomposition** — incremental R^2 of the image prediction over
   tumour purity, cancer type, tissue source site and stage.

Design constraints this package is built around
-----------------------------------------------
- Fully open, ungated, permissively-licensed data only (no dbGaP, no EGA,
  no non-commercial model weights).
- CPU-only for the core analysis; slide-level embeddings, ridge heads.
- Patient-level everything: splits, bootstraps, and unit of analysis.

Entry point: `aacr27.experiment.run_audit`.
"""

from . import (  # noqa: F401
    barcodes,
    data,
    decomposition,
    experiment,
    globalaxis,
    models,
    outcome,
    signatures,
    splits,
    stats,
)
from .experiment import AuditConfig, AuditResult, run_audit  # noqa: F401

__version__ = "0.1.0"

__all__ = [
    "barcodes",
    "data",
    "decomposition",
    "experiment",
    "globalaxis",
    "models",
    "outcome",
    "signatures",
    "splits",
    "stats",
    "AuditConfig",
    "AuditResult",
    "run_audit",
]
