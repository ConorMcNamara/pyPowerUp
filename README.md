# pyPowerUp
A Python implementation of [PowerUpR](https://cran.r-project.org/web/packages/PowerUpR/PowerUpR.pdf); a library for calculating the power, sample size and minimum detectable effect size of Multilevel Randomized Experiments.

## Quick Example
``` 
from pyPowerUp.power import power_bcra3f2
power = power_bcra3f2(effect_size=.145, rho2=.10, n=20, J=44, K=5, alpha=0.05)
round(power, 3)
0.803
```

## Notes
Whenever possible, I tried to follow the R naming and code-style to ensure as much 1-1 comparison as possible; however, some liberties were taken to ensure the code follows PEP-8 guidelines. 
