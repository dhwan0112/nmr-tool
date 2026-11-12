"""Tools for reading and analysing processed Bruker 1H NMR spectra."""

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
