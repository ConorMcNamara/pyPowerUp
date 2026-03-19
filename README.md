# pyPowerUp

[![CI](https://github.com/ConorMcNamara/pyPowerUp/actions/workflows/python-package.yml/badge.svg)](https://github.com/ConorMcNamara/pyPowerUp/actions/workflows/python-package.yml)
[![Python Version](https://img.shields.io/badge/python-3.13%20%7C%203.14-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A modern Python implementation of [PowerUpR](https://cran.r-project.org/web/packages/PowerUpR/PowerUpR.pdf) for calculating statistical power, sample size, and minimum detectable effect size of multilevel randomized experiments.

## Features

Includes tools to calculate statistical power, minimum detectable effect size (MDES), MDES difference (MDESD), and minimum required sample size for various multilevel randomized experiments (MRE) with continuous outcomes:

- **14 types** of MRE designs to detect main treatment effect
- **7 types** of MRE designs to detect moderated treatment effect (2-1-1, 2-1-2, 2-2-1, 2-2-2, 3-3-1, 3-3-2, and 3-3-3 designs)
- **5 types** of MRE designs to detect mediated treatment effects (2-1-1, 2-2-1, 3-1-1, 3-2-1, and 3-3-1 designs)
- **4 types** of partially nested (PN) designs to detect main treatment effect
- **3 types** of PN designs to detect mediated treatment effects

See the ['PowerUp!' Excel series](https://www.causalevaluation.org/) for more information.

## Installation

### From PyPI (when published)

```bash
pip install pypowerup
```

### From source

```bash
git clone https://github.com/ConorMcNamara/pyPowerUp.git
cd pyPowerUp
pip install -e .
```

### For development

```bash
git clone https://github.com/ConorMcNamara/pyPowerUp.git
cd pyPowerUp
pip install -e ".[dev]"
pre-commit install
```

### Generating requirements.txt (optional)

If you need a `requirements.txt` file for legacy systems or Docker:

```bash
pip install pip-tools
pip-compile pyproject.toml
```

## Requirements

- Python 3.13 or 3.14
- NumPy >= 1.26.0, < 3.0.0
- SciPy >= 1.15.2, < 2.0.0

**Note**: All dependencies are managed in `pyproject.toml` following modern Python packaging standards (PEP 621).

## Quick Example

```python
from pyPowerUp.power import power_bcra3f2

# Calculate power for a three-level blocked cluster randomized trial
power = power_bcra3f2(
    effect_size=0.145,
    rho2=0.10,
    n=20,
    J=44,
    K=5,
    alpha=0.05
)
print(round(power, 3))
# Output: 0.803
```

## Available Functions

### Power Analysis
Calculate statistical power for various experimental designs:
- `power_bcra3f2`, `power_bcra3r2`, `power_bcra4f3`, `power_bcra4r2`, `power_bcra4r3`
- `power_bira2c1`, `power_bira2f1`, `power_bira2r1`, `power_bira3r1`, `power_bira4r1`
- `power_cra2r2`, `power_cra3r3`, `power_cra4r4`
- `power_ira1r1`

### Minimum Detectable Effect Size
Calculate minimum detectable effect sizes:
- `mde_bcra3f2`, `mde_bcra3r2`, `mde_bcra4f3`, `mde_bcra4r2`, `mde_bcra4r3`
- `mde_bira2c1`, `mde_bira2f1`, `mde_bira2r1`, `mde_bira3r1`, `mde_bira4r1`
- `mde_cra2r2`, `mde_cra3r3`, `mde_cra4r4`
- `mde_ira1r1`

### Sample Size Calculation
Calculate minimum required sample sizes:
- `sample_size_bcra3f2`, `sample_size_bcra3r2`, `sample_size_bcra4f3`, `sample_size_bcra4r2`, `sample_size_bcra4r3`
- `sample_size_bira2c1`, `sample_size_bira2f1`, `sample_size_bira2r1`, `sample_size_bira3r1`, `sample_size_bira4r1`
- `sample_size_cra2r2`, `sample_size_cra3r3`, `sample_size_cra4r4`
- `sample_size_ira1r1`

## Design Philosophy

This implementation aims to maintain compatibility with the original R package (PowerUpR) in terms of naming conventions and functionality, while adhering to Python's PEP 8 style guidelines and modern best practices. The package includes full type hints for better IDE support and type checking.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate and ensure all tests pass before submitting.

## Development

```bash
# Run tests
pytest

# Run tests with coverage
pytest --cov=pyPowerUp --cov-report=html

# Format code
ruff format .

# Lint code
ruff check .

# Type check
zuban check pyPowerUp
```

## References

Dong, N. & Maynard, R. A. (2013). PowerUp!: A tool for calculating minimum detectable effect sizes and minimum required sample sizes for experimental and quasi-experimental design studies. *Journal of Research on Educational Effectiveness*, 6(1), 24-67. [doi:10.1080/19345747.2012.673143](https://doi.org/10.1080/19345747.2012.673143)

Bulus, M., Dong, N., Kelcey, B., & Spybrook, J. (2019). PowerUpR: Power Analysis Tools for Multilevel Randomized Experiments. R package version 1.1.0. [CRAN](https://CRAN.R-project.org/package=PowerUpR)

## Related Projects

- [PowerUpR (R implementation)](https://CRAN.R-project.org/package=PowerUpR)
- [pypowerup (alternative Python implementation)](https://github.com/sophiamyang/pypowerup)

## License

This project is licensed under the MIT License - see the LICENSE file for details.
