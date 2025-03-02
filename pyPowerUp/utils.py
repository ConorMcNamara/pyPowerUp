from typing import Dict
from scipy.stats import t as t_dist, nct
from math import sqrt


def _mde(power: float, alpha: float, sse: float, df: int, two_tailed: bool) -> Dict:
    """Calculates the mde of the test

    Parameters
    ----------
    power: float
        The power of the test
    alpha: float
        The significance level of the test
    sse: float
        The sum of squared errors of the test
    df: int
        The degrees of freedom of the test
    two_tailed: bool
        Whether the test is one-tailed or two-tailed

    Returns
    -------
    A dictionary containing the Minimum Detectable Effect and the Confidence Intervals around said effect
    """
    if sse < 0:
        raise ValueError("Sum of Squared Error cannot be less than 0")
    if df < 1:
        raise ValueError("degrees of freedom must be at least 1")
    t1 = abs(t_dist.ppf(alpha / 2, df)) if two_tailed else abs(t_dist.ppf(alpha, df))
    t2 = abs(t_dist.ppf(power, df))
    m = t1 + t2 if power >= 0.5 else t1 - t2
    mde = m * sse
    lower_bound = mde * (1 - t1 / m)
    upper_bound = mde * (1 + t1 / m)
    return {'minimum_detectable_effect': mde,
            f'{int((1 - round(alpha, 2)) * 100)}% Confidence Interval': [lower_bound, upper_bound]}


def _power(effect_size: float, alpha: float, sse: float, df: float, two_tailed: bool) -> float:
    """Calculates the power of the test

    Parameters
    ----------
    effect_size: float
        The effect size of the test
    alpha: float
        The significance level of the test
    sse: float
        The sum of squared errors of the test
    df: int
        The degrees of freedom of the test
    two_tailed: bool
        Whether the test is one-tailed or two-tailed

    Returns
    -------
    The power of the test
    """
    if sse < 0:
        raise ValueError("Sum of Squared Error cannot be less than 0")
    if df < 1:
        raise ValueError("degrees of freedom must be at least 1")
    lamda = effect_size / sse
    if two_tailed:
        power = 1 - nct.cdf(t_dist.isf(alpha / 2, df), df, lamda) + nct.cdf(-t_dist.isf(alpha / 2, df), df, lamda)
    else:
        power = 1 - nct.cdf(t_dist.isf(alpha, df), df, lamda)
    return power


def _sse_a221(esa: float, r2m2: float, p: float, J: float) -> float:
    var_a221 = (1 - (r2m2 + p * (1 - p) * pow(esa, 2))) / (p * (1 - p) * J)
    if var_a221 < 0:
        raise ValueError("Variance cannot be less than 0")
    return sqrt(var_a221)


def _se_b221(esa: float, esb: float, escp: float, rho2: float, r22: float, r21: float, r2m2: float, p: float, n: float,
             J: float) -> float:
    var_b221 = (rho2 * (1 - (r22 + p * (1 - p) * pow(esa * esb + escp, 2) / rho2 + (pow(esb, 2) / rho2) * (
            1 - r2m2 - p * (1 - p) * pow(esa, 2)))) +
                (1 - rho2) * (1 - r21) / n) / (J * (1 - (r2m2 + p * (1 - p) * pow(esa, 2))))
    if var_b221 < 0:
        raise ValueError("Variance cannot be less than 0")
    return sqrt(var_b221)


def _se_a211(esa: float, rhom2: float, r2m1: float, r2m2: float, n: int, J: float, p: float) -> float:
    t2mbar = rhom2 * (1 - r2m2 - (p * (1 - p) * pow(esa, 2)) / rhom2)
    sig2mbar = (1 - rhom2) * (1 - r2m1)
    var_a211 = (t2mbar + sig2mbar / n) / (J * p * (1 - p))
    if var_a211 < 0:
        raise ValueError("Variance cannot be less than 0")
    return sqrt(var_a211)


def _se_b1211(esb1: float, rho2: float, rhom2: float, r21: float, r2m1: float, n: int, J: float) -> float:
    sig2mbar = (1 - rhom2) * (1 - r2m1)
    sig2ybar = (1 - rho2) * (1 - r21 - (((1 - rhom2) / (1 - rho2)) * pow(esb1, 2) * (1 - r2m1)))
    var_b1211 = sig2ybar / ((J * n - J) * sig2mbar)
    if var_b1211 < 0:
        raise ValueError("Variance cannot be less than 0")
    return sqrt(var_b1211)


def _se_b211(esa: float, esB: float, esb1: float, escp: float, rho2: float, rhom2: float, r22: float, r21: float,
             r2m2: float, r2m1: float, n: float, J: float, p: float) -> float:
    t2mbar = rhom2 * (1 - r2m2 - (p * (1 - p) * pow(esa, 2)) / rhom2)
    sig2mbar = (1 - rhom2) * (1 - r2m1)
    t2ybar = rho2 * (1 - r22) - p * (1 - p) * pow(esa * esB + escp, 2) - \
             ((1 / (p * (1 - p))) * pow(esB, 2) * rhom2 * (1 - r2m2) +
              (1 / (p * (1 - p))) * pow(esB, 2) * (1 - rhom2) * (1 - r2m1) / n - pow(esa, 2) * pow(esB, 2)) / (
                     1 / (p * (1 - p)))
    sig2ybar = (1 - rho2) * (1 - r21 - (((1 - rhom2) / (1 - rho2)) * pow(esb1, 2) * (1 - r2m1)))
    var_b211 = (t2ybar + sig2ybar / n) / (J * (t2mbar + sig2mbar / n))
    if var_b211 < 0:
        raise ValueError("Variance cannot be less than 0")
    return sqrt(var_b211)


def _se_a321(rhom3: float, r2m2: float, r2m3: float, p: float, J: float, K: int):
    var_a321 = (rhom3 * (1 - r2m3) + (1 - rhom3) * (1 - r2m2) / J) / \
               (p * (1 - p) * (K - 5))
    if var_a321 < 0:
        raise ValueError("Variance cannot be less than 0")
    return sqrt(var_a321)


def _se_b321(rhom3: float, rho2: float, rho3: float, r2m2: float, r2m3: float, r21: float, r22: float, r23: float,
             p: float, n: int, J: float, K: int) -> float:
    var_b321 = (rho3 * (1 - r23) + rho2 * (1 - r22) / J + (1 - rho3 - rho2) * (1 - r21) / (n * J)) / \
               ((K - 6) * (rhom3 * (1 - r2m3) + (1 - rhom3) * (1 - r2m2) / J))
    if var_b321 < 0:
        raise ValueError("Variance cannot be less than 0")
    return sqrt(var_b321)


def _se_sobel(x: float, y: float, se_x: float, se_y: float) -> float:
    var_sobel = pow(x, 2) + pow(y, 2) + pow(se_x, 2) + pow(se_y, 2)
    if var_sobel < 0:
        raise ValueError("Variance cannot be less than 0")
    return sqrt(var_sobel)


def _power_sobel(x: float, y: float, se_x: float, se_y: float, alpha: float, two_tailed: bool, df: int = 1e08) -> float:
    se_sobel = _se_sobel(x, y, se_x, se_y)
    power = _power(x * y, alpha, se_sobel, df, two_tailed)
    return power


def _power_jt(x: float, y: float, se_x: float, se_y: float, alpha: float, two_tailed: bool, df_x: int,
              df_y: int) -> float:
    power_x = _power(x, alpha, se_x, df_x, two_tailed)
    power_y = _power(y, alpha, se_y, df_y, two_tailed)
    return power_x * power_y
