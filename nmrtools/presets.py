"""Fixed settings for specific experiments, shared by the scripts and the interactive shell."""

from __future__ import annotations

# Molar masses in g/mol.
EXPT4_MW = {"linalool": 154.25, "ketone": 126.20, "side product": 128.21, "EtOAc": 88.11}

# One well-resolved, non-overlapping resonance per species (handout, section 1).
# (label, window low, window high, number of H)
EXPT4_QUANT_SIGNALS = {
    "linalool": ("CH=CH2 internal H, δ 5.91", 5.80, 6.00, 1),
    "side product": ("CH(OH), δ 3.81", 3.70, 3.92, 1),
    "EtOAc": ("OCH2, δ 4.12", 4.05, 4.20, 2),
}
