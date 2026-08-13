"""dNK chemokine profiling — dataset-agnostic pipeline.

Rules enforced in code (see docs/PREREGISTRATION.md):
  R1  every p value consumes one value per donor
  R2  min_achievable_p reported alongside every test
  R3  detection rate is the primary readout
  R10 no scanpy
  R11 panel match rate below threshold is a hard failure
"""
__version__ = "1.0.0"
PANEL_VERSION = "chemokine_panel_v1"
