"""Tools for reading processed Bruker 1H NMR spectra and for the routine
analysis of the CM3191 Expt 4 (linalool) spectra."""

from .analysis import Component, coupling_constants, integrate, mass_correction, pick_lines
from .bruker import Spectrum, load_spectrum

__all__ = [
    "Component",
    "Spectrum",
    "coupling_constants",
    "integrate",
    "load_spectrum",
    "mass_correction",
    "pick_lines",
]
