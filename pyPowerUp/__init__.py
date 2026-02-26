"""Performing power, sample size and minimum detectable effect size calculations of Multilevel Randomized Experiments.

pyPowerUp is a Python implementation of PowerUpR for calculating statistical power,
minimum detectable effect size (MDES), and minimum required sample size for various
multilevel randomized experiments (MRE) with continuous outcomes.

Modules
-------
power : Power analysis functions
    Calculate statistical power for various experimental designs
mde : Minimum detectable effect size functions
    Calculate minimum detectable effect sizes for various experimental designs
sample_size : Sample size calculation functions
    Calculate minimum required sample sizes for various experimental designs
utils : Utility functions
    Internal utility functions for calculations

Example
-------
>>> from pyPowerUp.power import power_bcra3f2
>>> power = power_bcra3f2(effect_size=.145, rho2=.10, n=20, J=44, K=5, alpha=0.05)
>>> round(power, 3)
0.803
"""

__version__ = "1.0.0"
__author__ = "Conor McNamara"
__email__ = "conor.s.mcnamara@gmail.com"

from pyPowerUp import mde, power, sample_size, utils

__all__ = [
    # Modules
    "mde",
    "power",
    "sample_size",
    "utils",
    # Version info
    "__version__",
    "__author__",
    "__email__",
]
