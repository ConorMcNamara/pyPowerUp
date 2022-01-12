# pyPowerUp
A Python implementation of [PowerUpR](https://cran.r-project.org/web/packages/PowerUpR/PowerUpR.pdf); a library for calculating the power, sample size and minimum detectable effect size of Multilevel Randomized Experiments.

To quote the documentation
> Includes tools to calculate statistical power, minimum detectable effect size (MDES), MDES difference (MDESD), and minimum required sample size for various multilevel randomized experiments (MRE) with continuous outcomes. Accomodates 14 types of MRE designs to detect main treatment effect, seven types of MRE designs to detect moderated treatment effect (2-1-1, 2-1-2, 2-2-1, 2-2-2, 3-3-1, 3-3-2, and 3-3-3 designs; <total.lev> - <trt.lev> - <mod.lev>), five types of MRE designs to detect mediated treatment effects (2-1-1, 2-2-1, 3-1-1, 3-2-1, and 3-3-1 designs; <trt.lev> - <med.lev> - <out.lev>), four types of partially nested (PN) design to detect main treatment effect, and three types of PN designs to detect mediated treatment effects (2/1, 3/1, 3/2; <trt.arm.lev> / <ctrl.arm.lev>). See 'PowerUp!' Excel series at <https://www.causalevaluation.org/>.

## Quick Example
``` 
from pyPowerUp.power import power_bcra3f2
power = power_bcra3f2(effect_size=.145, rho2=.10, n=20, J=44, K=5, alpha=0.05)
round(power, 3)
0.803
```

## Notes
Whenever possible, I tried to follow the R naming and code-style to ensure as much 1-1 comparison as possible; however, some liberties were taken to ensure the code follows PEP-8 guidelines. 
