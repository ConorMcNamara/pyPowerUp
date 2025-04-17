from math import sqrt, ceil
from typing import Dict

from pyPowerUp.utils import _mde


def mde_bcra3f2(
    rho2: float,
    n: float,
    J: float,
    K: int,
    power: float = 0.80,
    alpha: float = 0.10,
    two_tailed: bool = True,
    p: float = 0.50,
    g2: int = 0,
    r21: int = 0,
    r22: int = 0,
    print_pretty: bool = True,
) -> Dict:
    """Calculates the Minimum Detectable Effect of a Three-Level Blocked (Fixed) Cluster-level Random Assignment Design,
    Treatment at Level 2

    Parameters
    ----------
    rho2 : float
        Proportion of variance in the outcome between level 2 units (unconditional ICC2)
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    K : int
        Number of level 3 units
    power : float, default=0.8
        Statistical power of the test.
    alpha : float, default=0.10
        Probability of Type I error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.5
        Average proportion of level 2 units randomly assigned to treatment within level 3 units
    g2 : int, default=0
        Number of covariates at level 2
    r21 : float, default=0
        Proportion of level 1 variance in the outcome explained by level 1 covariates
    r22 : float, default=0
        Proportion of level 2 variance in the outcome explained by level 2 covariates
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    A dictionary containing the minimum_detectable effect as well as confidence intervals for said effect
    """
    df = ceil(K * (J - 2) - g2)
    sse = sqrt(rho2 * (1 - r22) / (p * (1 - p) * J * K) + (1 - rho2) * (1 - r21) / (p * (1 - p) * J * K * n))
    mde = _mde(power, alpha, sse, df, two_tailed)
    if print_pretty:
        confidence_intervals = [round(i, 3) for i in mde[f"{int((1 - round(alpha, 2)) * 100)}% Confidence Interval"]]
        str_print = (
            "Minimum Detectable Effect Size"
            + "\n"
            + "-" * 39
            + "\n"
            + f" {round(mde['minimum_detectable_effect'], 3)} {int((1 - round(alpha, 2)) * 100)}% CI {confidence_intervals}"
            + "\n"
            + "-" * 39
            + "\n"
            + f"Degrees of Freedom: {df}"
            + "\n"
            + f"Standardized Standard Error: {round(sse, 3)}"
            + "\n"
            + f"Type I Error Rate: {round(alpha, 2)}"
            + "\n"
            + f"Type II Error Rate: {round(1 - power, 2)}"
            + "\n"
            + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return mde


def mde_bcra3r2(
    rho2: float,
    rho3: float,
    omega3: float,
    n: float,
    J: float,
    K: int,
    power: float = 0.80,
    alpha: float = 0.10,
    two_tailed: bool = True,
    p: float = 0.50,
    g3: int = 0,
    r21: int = 0,
    r22: int = 0,
    r2t3: float = 0,
    print_pretty: bool = True,
) -> Dict:
    """Calculates the Minimum Detectable Effect of a Three-Level Blocked Cluster-level Random Assignment Design,
    Treatment at Level 2

    Parameters
    ----------
    rho2 : float
        Proportion of variance in the outcome between level 2 units (unconditional ICC2)
    rho3 : float
        Proportion of variance in the outcome between level 3 units (unconditional ICC3)
    omega3 : float
        Treatment effect heterogeneity as ratio of treatment effect variance among level 3 units to the residual
        variance at level 3
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    K : int
        Number of level 3 units
    power : float, default=0.8
        Statistical power
    alpha : float, default=0.1
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.5
        Average proportion of level 2 units randomly assigned to the treatment within level 3 units
    g3 : int, default=0
        Number of covariates at level 3
    r21 : float, default=0
        Proportion of level 1 variance in the outcome explained by level 1 covariates
    r22 : float, default=0
        Proportion of level 2 variance in the outcome explained by level 2 covariates
    r2t3 : float, default=0
        Proportion of treatment effect variance among level 3 units explained by level 3 covariates
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    A dictionary containing the minimum_detectable effect as well as confidence intervals for said effect
    """
    df = K - g3 - 1
    sse = sqrt(
        rho3 * omega3 * (1 - r2t3) / K
        + rho2 * (1 - r22) / (p * (1 - p) * J * K)
        + (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * n)
    )
    mde = _mde(power, alpha, sse, df, two_tailed)
    if print_pretty:
        confidence_intervals = [round(i, 3) for i in mde[f"{int((1 - round(alpha, 2)) * 100)}% Confidence Interval"]]
        str_print = (
            "Minimum Detectable Effect Size"
            + "\n"
            + "-" * 39
            + "\n"
            + f" {round(mde['minimum_detectable_effect'], 3)} {int((1 - round(alpha, 2)) * 100)}% CI {confidence_intervals}"
            + "\n"
            + "-" * 39
            + "\n"
            + f"Degrees of Freedom: {df}"
            + "\n"
            + f"Standardized Standard Error: {round(sse, 3)}"
            + "\n"
            + f"Type I Error Rate: {round(alpha, 2)}"
            + "\n"
            + f"Type II Error Rate: {round(1 - power, 2)}"
            + "\n"
            + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return mde


def mde_bcra4f3(
    rho2: float,
    rho3: float,
    n: int,
    J: int,
    K: int,
    L: int,
    power: float = 0.80,
    alpha: float = 0.10,
    two_tailed: bool = True,
    p: float = 0.50,
    r21: int = 0,
    r22: int = 0,
    r23: int = 0,
    g3: int = 0,
    print_pretty: bool = True,
) -> Dict:
    """Calculates the Minimum Detectable Effect of a Four-Level Blocked (Fixed) Cluster-level Random Assignment Design,
    Treatment at Level 3

    Parameters
    ----------
    rho2 : float
        Proportion of variance in the outcome between level 2 units (unconditional ICC2)
    rho3 : float
        Proportion of variance in the outcome between level 3 units (unconditional ICC3)
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    K : float
        Harmonic mean of level 3 units across level 4 units (or simple average)
    L : int
        Number of level 4 units
    power : float, default=0.8
        Statistical power
    alpha : float, default=0.1
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.5
        Average proportion of level 3 units randomly assigned to the treatment within level 4 units
    r21 : float, default=0
        Proportion of level 1 variance in the outcome explained by level 1 covariates
    r22 : float, default=0
        Proportion of level 2 variance in the outcome explained by level 2 covariates
    r23 : float, default=0
        Proportion of level 3 variance in the outcome explained by level 3 covariates
    g3 : int, default=0
        Number of covariates at level 3
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    A dictionary containing the minimum_detectable effect as well as confidence intervals for said effect
    """
    df = L * (K - 2) - g3
    sse = sqrt(
        rho3 * (1 - r23) / (p * (1 - p) * K * L)
        + rho2 * (1 - r22) / (p * (1 - p) * J * K * L)
        + (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * L * n)
    )
    mde = _mde(power, alpha, sse, df, two_tailed)
    if print_pretty:
        confidence_intervals = [round(i, 3) for i in mde[f"{int((1 - round(alpha, 2)) * 100)}% Confidence Interval"]]
        str_print = (
            "Minimum Detectable Effect Size"
            + "\n"
            + "-" * 39
            + "\n"
            + f" {round(mde['minimum_detectable_effect'], 3)} {int((1 - round(alpha, 2)) * 100)}% CI {confidence_intervals}"
            + "\n"
            + "-" * 39
            + "\n"
            + f"Degrees of Freedom: {df}"
            + "\n"
            + f"Standardized Standard Error: {round(sse, 3)}"
            + "\n"
            + f"Type I Error Rate: {round(alpha, 2)}"
            + "\n"
            + f"Type II Error Rate: {round(1 - power, 2)}"
            + "\n"
            + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return mde


def mde_bcra4r2(
    rho2: float,
    rho3: float,
    rho4: float,
    omega3: float,
    omega4: float,
    n: int,
    J: int,
    K: int,
    L: int,
    power: float = 0.80,
    alpha: float = 0.10,
    two_tailed: bool = True,
    p: float = 0.50,
    r21: int = 0,
    r22: int = 0,
    r2t3: int = 0,
    r2t4: int = 0,
    g4: int = 0,
    print_pretty: bool = True,
) -> Dict:
    """Calculates the Minimum Detectable Effect of a Four-Level Blocked Cluster-level Random Assignment Design, Treatment
    at Level 2

    Parameters
    ----------
    rho2 : float
        Proportion of variance in the outcome between level 2 units (unconditional ICC2)
    rho3 : float
        Proportion of variance in the outcome between level 3 units (unconditional ICC3)
    rho4 : float
        Proportion of variance in the outcome between level 4 untis (unconditional ICC4)
    omega3 : float
        Treatment effect heterogeneity as ratio of treatment effect variance among level 3 units to the residual
        variance at level 3
    omega4 : float
        Treatment effect heterogeneity as ratio of treatment effect variance among level 4 units to the residual
        variance at level 4
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    K : float
        Harmonic mean of level 3 units across level 4 units (or simple average)
    L : int
        Number of level 4 units
    power : float, default=0.8
        Statistical power
    alpha : float, default=0.1
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.5
        Average proportion of level 2 units randomly assignewd to treatment within level 3 units
    r21 : float, default=0
        Proportion of level 1 variance in the outcome explained by level 1 covariates
    r22 : float, default=0
        Proportion of level 2 variance in the outcome explained by level 2 covariates
    r2t3 : float, default=0
        Proportion of treatment effect variance among level 3 units explained by level 3 covariates
    r2t4 : float, default=0
        Proportion of treatment effect variance among level 4 units explained by level 4 covariates
    g4 : int
        Number of covariates at level 4
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    A dictionary containing the minimum_detectable effect as well as confidence intervals for said effect
    """
    df = L - g4 - 1
    sse = sqrt(
        rho4 * omega4 * (1 - r2t4) / L
        + rho3 * omega3 * (1 - r2t3) / (K * L)
        + rho2 * (1 - r22) / (p * (1 - p) * J * K * L)
        + (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * L * n)
    )
    mde = _mde(power, alpha, sse, df, two_tailed)
    if print_pretty:
        confidence_intervals = [round(i, 3) for i in mde[f"{int((1 - round(alpha, 2)) * 100)}% Confidence Interval"]]
        str_print = (
            "Minimum Detectable Effect Size"
            + "\n"
            + "-" * 39
            + "\n"
            + f" {round(mde['minimum_detectable_effect'], 3)} {int((1 - round(alpha, 2)) * 100)}% CI {confidence_intervals}"
            + "\n"
            + "-" * 39
            + "\n"
            + f"Degrees of Freedom: {df}"
            + "\n"
            + f"Standardized Standard Error: {round(sse, 3)}"
            + "\n"
            + f"Type I Error Rate: {round(alpha, 2)}"
            + "\n"
            + f"Type II Error Rate: {round(1 - power, 2)}"
            + "\n"
            + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return mde


def mde_bcra4r3(
    rho2: float,
    rho3: float,
    rho4: float,
    omega4: float,
    n: float,
    J: float,
    K: int,
    L: int,
    power: float = 0.8,
    alpha: float = 0.10,
    two_tailed: bool = True,
    p: float = 0.50,
    r21: int = 0,
    r22: int = 0,
    r23: int = 0,
    r2t4: int = 0,
    g4: int = 0,
    print_pretty: bool = True,
) -> Dict:
    """Calculates the Minimum Detectable Effect of a Four-Level Blocked Cluster-level Random Assignment Design, Treatment
     at Level 3

    Parameters
    ----------
    rho2 : float
        Proportion of variance in the outcome between level 2 units (unconditional ICC2)
    rho3 : float
        Proportion of variance in the outcome between level 3 units (unconditional ICC3)
    rho4 : float
        Proportion of variance in the outcome between level 4 units (unconditional ICC4)
    omega4 : float
        Treatment effect heterogeneity as ratio of treatment effect variance among level 4 units to the residual
        variance at level 4
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    K : float
        Harmonic mean of level 3 units across level 4 units (or simple average)
    L : int
        Number of level 4 units
    power : float, default=0.8
        Statistical power
    alpha : float, default=0.1
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0
        Average proportion of level 3 units randomly assigned to treatment within level 4 units
    r21 : float, default=0
        Proportion of level 1 variance in the outcome explained by level 1 covariates
    r22 : float, default=0
        Proportion of level 2 variance in the outcome explained by level 2 covariates
    r23 : float, default=0
        Proportion of level 3 variance in the outcome explained by level 3 covariates
    r2t4 : float, default=0
        Proportion of treatment effect variance among level 4 units explained by level 4 covariates
    g4 : int, default=0
        Number of covariates at level 4
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    A dictionary containing the minimum_detectable effect as well as confidence intervals for said effect
    """
    df = L - g4 - 1
    sse = sqrt(
        rho4 * omega4 * (1 - r2t4) / L
        + rho3 * (1 - r23) / (p * (1 - p) * K * L)
        + rho2 * (1 - r22) / (p * (1 - p) * J * K * L)
        + (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * L * n)
    )
    mde = _mde(power, alpha, sse, df, two_tailed)
    if print_pretty:
        confidence_intervals = [round(i, 3) for i in mde[f"{int((1 - round(alpha, 2)) * 100)}% Confidence Interval"]]
        str_print = (
            "Minimum Detectable Effect Size"
            + "\n"
            + "-" * 39
            + "\n"
            + f" {round(mde['minimum_detectable_effect'], 3)} {int((1 - round(alpha, 2)) * 100)}% CI {confidence_intervals}"
            + "\n"
            + "-" * 39
            + "\n"
            + f"Degrees of Freedom: {df}"
            + "\n"
            + f"Standardized Standard Error: {round(sse, 3)}"
            + "\n"
            + f"Type I Error Rate: {round(alpha, 2)}"
            + "\n"
            + f"Type II Error Rate: {round(1 - power, 2)}"
            + "\n"
            + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return mde


def mde_bira2c1(
    n: float,
    J: float,
    power: float = 0.80,
    alpha: float = 0.10,
    two_tailed: bool = True,
    p: float = 0.5,
    g1: int = 0,
    r21: int = 0,
    print_pretty: bool = True,
) -> Dict:
    """Calculates the Minimum Detectable Effect of a Two-Level Blocked (Constant Treatment Effect) Individual-level
    Random Assignment Design, Treatment at Level 1

    Parameters
    ----------
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    power : float, default=0.8
        Statistical power
    alpha : float, default=0.1
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.5
        Average proportion of level 1 units randomly assigned to treatment within level 2 units
    g1 : int, default=0
        Number of covariates at level 1
    r21 : float, default=0
        Proportion of level 1 variance in the outcome explained by level 1 covariates
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    A dictionary containing the minimum_detectable effect as well as confidence intervals for said effect
    """
    df = ceil(J * (n - 1) - g1 - 1)
    sse = sqrt((1 - r21) / (p * (1 - p) * J * n))
    mde = _mde(power, alpha, sse, df, two_tailed)
    if print_pretty:
        confidence_intervals = [round(i, 3) for i in mde[f"{int((1 - round(alpha, 2)) * 100)}% Confidence Interval"]]
        str_print = (
            "Minimum Detectable Effect Size"
            + "\n"
            + "-" * 39
            + "\n"
            + f" {round(mde['minimum_detectable_effect'], 3)} {int((1 - round(alpha, 2)) * 100)}% CI {confidence_intervals}"
            + "\n"
            + "-" * 39
            + "\n"
            + f"Degrees of Freedom: {df}"
            + "\n"
            + f"Standardized Standard Error: {round(sse, 3)}"
            + "\n"
            + f"Type I Error Rate: {round(alpha, 2)}"
            + "\n"
            + f"Type II Error Rate: {round(1 - power, 2)}"
            + "\n"
            + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return mde


def mde_bira2f1(
    n: int,
    J: int,
    power: float = 0.80,
    alpha: float = 0.10,
    two_tailed: bool = True,
    p: float = 0.50,
    g1: int = 0,
    r21: int = 0,
    print_pretty: bool = True,
) -> Dict:
    """Calculates the Minimum Detectable Effect of a Two-Level Blocked (Fixed) Individual-level Random Assignment Design,
    Treatment at Level 1

    Parameters
    ----------
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    power : float, default=0.8
        Statistical power
    alpha : float, default=0.1
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.5
        Average proportion of level 1 units randomly assigned to treatment within level 2 units
    g1 : int, default=0
        Number of covariates at level 1
    r21 : float, default=0
        Proportion of level 1 variance in the outcome explained by level 1 covariates
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    A dictionary containing the minimum_detectable effect as well as confidence intervals for said effect
    """
    df = J * (n - 2) - g1
    sse = sqrt((1 - r21) / (p * (1 - p) * J * n))
    mde = _mde(power, alpha, sse, df, two_tailed)
    if print_pretty:
        confidence_intervals = [round(i, 3) for i in mde[f"{int((1 - round(alpha, 2)) * 100)}% Confidence Interval"]]
        str_print = (
            "Minimum Detectable Effect Size"
            + "\n"
            + "-" * 39
            + "\n"
            + f" {round(mde['minimum_detectable_effect'], 3)} {int((1 - round(alpha, 2)) * 100)}% CI {confidence_intervals}"
            + "\n"
            + "-" * 39
            + "\n"
            + f"Degrees of Freedom: {df}"
            + "\n"
            + f"Standardized Standard Error: {round(sse, 3)}"
            + "\n"
            + f"Type I Error Rate: {round(alpha, 2)}"
            + "\n"
            + f"Type II Error Rate: {round(1 - power, 2)}"
            + "\n"
            + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return mde


def mde_bira2r1(
    rho2: float,
    omega2: float,
    n: float,
    J: float,
    power: float = 0.80,
    alpha: float = 0.10,
    two_tailed: bool = True,
    p: float = 0.50,
    g2: int = 0,
    r21: int = 0,
    r2t2: int = 0,
    print_pretty: bool = True,
) -> Dict:
    """Calculates the Minimum Detectable Effect of a Two-Level Blocked Individual-level Random Assignment Design,
    Treatment at Level 1

    Parameters
    ----------
    rho2 : float
        Proportion of variance in the outcome between level 2 units (unconditional ICC2)
    omega2 : float
        Treatment effect heterogeneity as ratio of treatment effect variance among level 2 units to the residual
        variance at level 2
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    power : float, default=0.8
        Statistical power
    alpha : float, default=0.1
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.5
        Average proportion of level 1 units randomly assigned to treatment within level 2 units
    g2 : int, default=0
        Number of covariates at level 2
    r21 : float, default=0
        Proportion of level 1 variance in the outcome explained by level 1 covariates
    r2t2 : float, default=0
        Proportion of treatment effect variance among level 2 units explained by level 2 covariates
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    A dictionary containing the minimum_detectable effect as well as confidence intervals for said effect
    """
    df = ceil(J - g2 - 1)
    sse = sqrt(rho2 * omega2 * (1 - r2t2) / J + (1 - rho2) * (1 - r21) / (p * (1 - p) * J * n))
    mde = _mde(power, alpha, sse, df, two_tailed)
    if print_pretty:
        confidence_intervals = [round(i, 3) for i in mde[f"{int((1 - round(alpha, 2)) * 100)}% Confidence Interval"]]
        str_print = (
            "Minimum Detectable Effect Size"
            + "\n"
            + "-" * 39
            + "\n"
            + f" {round(mde['minimum_detectable_effect'], 3)} {int((1 - round(alpha, 2)) * 100)}% CI {confidence_intervals}"
            + "\n"
            + "-" * 39
            + "\n"
            + f"Degrees of Freedom: {df}"
            + "\n"
            + f"Standardized Standard Error: {round(sse, 3)}"
            + "\n"
            + f"Type I Error Rate: {round(alpha, 2)}"
            + "\n"
            + f"Type II Error Rate: {round(1 - power, 2)}"
            + "\n"
            + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return mde


def mde_bira3r1(
    rho2: float,
    rho3: float,
    omega2: float,
    omega3: float,
    n: float,
    J: float,
    K: int,
    power: float = 0.8,
    alpha: float = 0.10,
    two_tailed: bool = True,
    p: float = 0.50,
    r21: int = 0,
    r2t2: int = 0,
    r2t3: int = 0,
    g3: int = 0,
    print_pretty: bool = True,
) -> Dict:
    """Calculates the Minimum Detectable Effect of a

    Parameters
    ----------
    rho2 : float
        Proportion of variance in the outcome between level 2 units (unconditional ICC2)
    rho3 : float
        Proportion of variance in the outcome between level 3 units (unconditional ICC3)
    omega2 : float
        Treatment effect heterogeneity as ratio of treatment effect variance among level 2 units to the residual
        variance at level 2
    omega3 : float
        Treatment effect heterogeneity as ratio of treatment effect variance among level 3 units to the residual
        variance at level 3
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    K : int
        Number of level 3 units
    power : float, default=0.8
        Statistical power
    alpha : float, default=0.1
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0
        Average proportion of level 1 units randomly assigned to treatment within level 2 units
    r21 : float, default=0
        Proportion of level 1 variance in the outcome explained by level 1 covariates
    r2t2 : float, default=0
        Proportion of treatment effect variance among level 2 units explained by level 2 covariates
    r2t3 : float, default=0
        Proportion of treatment effect variance among level 3 units explaiend by level 3 covariates
    g3 : int, default=0
        Number of covariates at level 3
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    A dictionary containing the minimum_detectable effect as well as confidence intervals for said effect
    """
    df = K - g3 - 1
    sse = sqrt(
        rho3 * omega3 * (1 - r2t3) / K
        + rho2 * omega2 * (1 - r2t2) / (J * K)
        + (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * n)
    )
    mde = _mde(power, alpha, sse, df, two_tailed)
    if print_pretty:
        confidence_intervals = [round(i, 3) for i in mde[f"{int((1 - round(alpha, 2)) * 100)}% Confidence Interval"]]
        str_print = (
            "Minimum Detectable Effect Size"
            + "\n"
            + "-" * 39
            + "\n"
            + f" {round(mde['minimum_detectable_effect'], 3)} {int((1 - round(alpha, 2)) * 100)}% CI {confidence_intervals}"
            + "\n"
            + "-" * 39
            + "\n"
            + f"Degrees of Freedom: {df}"
            + "\n"
            + f"Standardized Standard Error: {round(sse, 3)}"
            + "\n"
            + f"Type I Error Rate: {round(alpha, 2)}"
            + "\n"
            + f"Type II Error Rate: {round(1 - power, 2)}"
            + "\n"
            + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return mde


def mde_bira4r1(
    rho2: float,
    rho3: float,
    rho4: float,
    omega2: float,
    omega3: float,
    omega4: float,
    n: float,
    J: float,
    K: int,
    L: int,
    power: float = 0.80,
    alpha: float = 0.10,
    two_tailed: bool = True,
    p: float = 0.50,
    r21: int = 0,
    r2t2: int = 0,
    r2t3: int = 0,
    r2t4: int = 0,
    g4: int = 0,
    print_pretty: bool = True,
) -> Dict:
    """Calculates the Minimum Detectable Effect of a Four-Level Blocked Individual-level Random Assignment Design,
    Treatment at Level 1

    Parameters
    ----------
    rho2 : float
        Proportion of variance in the outcome between level 2 units (unconditional ICC2)
    rho3 : float
        Proportion of variance in the outcome between level 3 units (unconditional ICC3)
    rho4 : float
        Proportion of variance in the outcome between level 4 units (unconditional ICC4)
    omega2 : float
        Treatment effect heterogeneity as ratio of treatment effect variance among level 2 units to the residual
        variance at level 2
    omega3 : float
        Treatment effect heterogeneity as ratio of treatment effect variance among level 3 units to the residual
        variance at level 3
    omega4 : float
        Treatment effect heterogeneity as ratio of treatment effect variance among level 4 units to the residual
        variance at level 4
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    K : float
        Harmonic mean of level 3 units across level 4 units (or simple average)
    L : int
        Number of level 4 units
    power : float, default=0.8
        Statistical power
    alpha : float, default=0.1
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.5
        Average proportion of level 1 units randomly assigned to treatment within level 2 units
    r21 : float, default=0
        Proportion of level 1 variance in the outcome explained by level 1 covariates
    r2t2 : float, default=0
        Proportion of treatment effect variance among level 2 units explained by level 2 covariates
    r2t3 : float, default=0
        Proportion of treatment effect variance among level 3 units explained by level 3 covariates
    r2t4 : float, default=0
        Proportion of treatment effect variance among level 4 units explained by level 4 covariates
    g4 : int, default=0
        Number of covariates at level 4
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    A dictionary containing the minimum_detectable effect as well as confidence intervals for said effect
    """
    df = L - g4 - 1
    sse = sqrt(
        rho4 * omega4 * (1 - r2t4) / L
        + rho3 * omega3 * (1 - r2t3) / (K * L)
        + rho2 * omega2 * (1 - r2t2) / (J * K * L)
        + (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * L * n)
    )
    mde = _mde(power, alpha, sse, df, two_tailed)
    if print_pretty:
        confidence_intervals = [round(i, 3) for i in mde[f"{int((1 - round(alpha, 2)) * 100)}% Confidence Interval"]]
        str_print = (
            "Minimum Detectable Effect Size"
            + "\n"
            + "-" * 39
            + "\n"
            + f" {round(mde['minimum_detectable_effect'], 3)} {int((1 - round(alpha, 2)) * 100)}% CI {confidence_intervals}"
            + "\n"
            + "-" * 39
            + "\n"
            + f"Degrees of Freedom: {df}"
            + "\n"
            + f"Standardized Standard Error: {round(sse, 3)}"
            + "\n"
            + f"Type I Error Rate: {round(alpha, 2)}"
            + "\n"
            + f"Type II Error Rate: {round(1 - power, 2)}"
            + "\n"
            + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return mde


def mde_cra2r2(
    rho2: float,
    n: float,
    J: float,
    power: float = 0.80,
    alpha: float = 0.05,
    two_tailed: bool = True,
    p: float = 0.50,
    g2: int = 0,
    r21: int = 0,
    r22: int = 0,
    print_pretty: bool = True,
) -> Dict:
    """Calculates the Minimum Detectable Effect of a Two-level Cluster-randomized Trials to Detect Main, Moderation and
    Mediation Effects

    Parameters
    ----------
    rho2 : float
        Proportion of variance in the outcome between level 2 units (unconditional ICC2)
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    power : float, default=0.8
        Statistical power
    alpha : float, default=0.1
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.5
        Proportion of level 2 units randomly assigned to treatment
    g2 : int, default=0
        Number of covariates at level 2
    r21 : float, default=0
        Proportion of level 1 variance in the outcome explained by level 1 covariates
    r22 : float, default=0
        Proportion of level 2 variance in the outcome explained by level 2 covariates
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    A dictionary containing the minimum_detectable effect as well as confidence intervals for said effect
    """
    df = ceil(J - g2 - 2)
    sse = sqrt(rho2 * (1 - r22) / (p * (1 - p) * J) + (1 - rho2) * (1 - r21) / (p * (1 - p) * J * n))
    mde = _mde(power, alpha, sse, df, two_tailed)
    if print_pretty:
        confidence_intervals = [round(i, 3) for i in mde[f"{int((1 - round(alpha, 2)) * 100)}% Confidence Interval"]]
        str_print = (
            "Minimum Detectable Effect Size"
            + "\n"
            + "-" * 39
            + "\n"
            + f" {round(mde['minimum_detectable_effect'], 3)} {int((1 - round(alpha, 2)) * 100)}% CI {confidence_intervals}"
            + "\n"
            + "-" * 39
            + "\n"
            + f"Degrees of Freedom: {df}"
            + "\n"
            + f"Standardized Standard Error: {round(sse, 3)}"
            + "\n"
            + f"Type I Error Rate: {round(alpha, 2)}"
            + "\n"
            + f"Type II Error Rate: {round(1 - power, 2)}"
            + "\n"
            + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return mde


def mde_cra3r3(
    rho2: float,
    rho3: float,
    n: float,
    J: float,
    K: int,
    power: float = 0.80,
    alpha: float = 0.10,
    two_tailed: bool = True,
    p: float = 0.50,
    g3: int = 0,
    r21: int = 0,
    r22: int = 0,
    r23: int = 0,
    print_pretty: bool = True,
) -> Dict:
    """Calculates the Minimum Detectable Effect of a Three-level Cluster-randomized Trials to Detect Main, Moderation,
    and Mediation Effects

    Parameters
    ----------
    rho2 : float
        Proportion of variance in the outcome between level 2 units (unconditional ICC2)
    rho3 : float
        Proportion of variance in the outcome between level 3 units (unconditional ICC3)
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    K : int
        Level 3 sample size
    power : float, default=0.8
        Statistical power
    alpha : float, default=0.1
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0
        Proportion of level 3 units randomly assigned to treatment
    g3 : int, default=0
        Number of covariates at level 3
    r21 : float, default=0
        Proportion of level 1 variance in the outcome explained by level 1 covariates
    r22 : float, default=0
        Proportion of level 2 variance in the outcome explained by level 2 covariates
    r23 : float, default=0
        Proportion of level 3 variance in the outcome explained by level 3 covariates
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    A dictionary containing the minimum_detectable effect as well as confidence intervals for said effect
    """
    df = K - g3 - 2
    sse = sqrt(
        rho3 * (1 - r23) / (p * (1 - p) * K)
        + rho2 * (1 - r22) / (p * (1 - p) * J * K)
        + (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * n)
    )
    mde = _mde(power, alpha, sse, df, two_tailed)
    if print_pretty:
        confidence_intervals = [round(i, 3) for i in mde[f"{int((1 - round(alpha, 2)) * 100)}% Confidence Interval"]]
        str_print = (
            "Minimum Detectable Effect Size"
            + "\n"
            + "-" * 39
            + "\n"
            + f" {round(mde['minimum_detectable_effect'], 3)} {int((1 - round(alpha, 2)) * 100)}% CI {confidence_intervals}"
            + "\n"
            + "-" * 39
            + "\n"
            + f"Degrees of Freedom: {df}"
            + "\n"
            + f"Standardized Standard Error: {round(sse, 3)}"
            + "\n"
            + f"Type I Error Rate: {round(alpha, 2)}"
            + "\n"
            + f"Type II Error Rate: {round(1 - power, 2)}"
            + "\n"
            + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return mde


def mde_cra4r4(
    rho2: float,
    rho3: float,
    rho4: float,
    n: float,
    J: float,
    K: float,
    L: int,
    power: float = 0.80,
    alpha: float = 0.10,
    two_tailed: bool = True,
    p: float = 0.50,
    r21: float = 0,
    r22: float = 0,
    r23: float = 0,
    r24: float = 0,
    g4: int = 0,
    print_pretty: bool = True,
) -> Dict:
    """Calculates the Minimum Detectable Effect of a Four-Level Cluster-randomized Trial

    Parameters
    ----------
    rho2 : float
        Proportion of variance in the outcome between level 2 units (unconditional ICC2)
    rho3 : float
        Proportion of variance in the outcome between level 3 units (unconditional ICC3)
    rho4 : float
        Proportion of variance in the outcome between level 4 units (unconditional ICC4)
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    K : float
        Harmonic mean of level 3 units across level 4 units (or simple average)
    L : int
        Number of level 4 units
    power : float, default=0.8
        Statistical power
    alpha : float, default=0.1
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.5
        Proportion of level 4 units randomly assigned to treatment
    r21 : float, default=0
        Proportion of level 1 variance in the outcome explained by level 1 covariates
    r22 : float, default=0
        Proportion of level 2 variance in the outcome explained by level 2 covariates
    r23 : float, default=0
        Proportion of level 3 variance in the outcome explained by level 3 covariates
    r24 : float, default=0
        Proportion of level 4 variance in the outcome explained by level 4 covariates
    g4 : int, default=0
        Number of covariates at level 4
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    A dictionary containing the minimum_detectable effect as well as confidence intervals for said effect
    """
    df = L - g4 - 2
    sse = sqrt(
        rho4 * (1 - r24) / (p * (1 - p) * L)
        + rho3 * (1 - r23) / (p * (1 - p) * K * L)
        + rho2 * (1 - r22) / (p * (1 - p) * J * K * L)
        + (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * L * n)
    )
    mde = _mde(power, alpha, sse, df, two_tailed)
    if print_pretty:
        confidence_intervals = [round(i, 3) for i in mde[f"{int((1 - round(alpha, 2)) * 100)}% Confidence Interval"]]
        str_print = (
            "Minimum Detectable Effect Size"
            + "\n"
            + "-" * 39
            + "\n"
            + f" {round(mde['minimum_detectable_effect'], 3)} {int((1 - round(alpha, 2)) * 100)}% CI {confidence_intervals}"
            + "\n"
            + "-" * 39
            + "\n"
            + f"Degrees of Freedom: {df}"
            + "\n"
            + f"Standardized Standard Error: {round(sse, 3)}"
            + "\n"
            + f"Type I Error Rate: {round(alpha, 2)}"
            + "\n"
            + f"Type II Error Rate: {round(1 - power, 2)}"
            + "\n"
            + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return mde


def mde_ira1r1(
    n: int,
    power: float = 0.80,
    alpha: float = 0.10,
    two_tailed: bool = True,
    p: float = 0.50,
    g1: int = 0,
    r21: float = 0,
    print_pretty: bool = True,
) -> Dict:
    """Calculates the Minimum Detectable Effect of a Individual-level Random Assignment Design

    Parameters
    ----------
    n : int
        Sample size
    power : float, default=0.8
        Statistical power
    alpha : float, default=0.1
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0
        Proportion of units randomly assigned to treatment
    g1 : int, default=0
        Number of covariates
    r21 : float, default=0
        Proportion of variance in the outcome explained by covariates
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    A dictionary containing the minimum_detectable effect as well as confidence intervals for said effect
    """
    df = n - g1 - 2
    sse = sqrt((1 - r21) / (p * (1 - p) * n))
    mde = _mde(power, alpha, sse, df, two_tailed)
    if print_pretty:
        confidence_intervals = [round(i, 3) for i in mde[f"{int((1 - round(alpha, 2)) * 100)}% Confidence Interval"]]
        str_print = (
            "Minimum Detectable Effect Size"
            + "\n"
            + "-" * 39
            + "\n"
            + f" {round(mde['minimum_detectable_effect'], 3)} {int((1 - round(alpha, 2)) * 100)}% CI {confidence_intervals}"
            + "\n"
            + "-" * 39
            + "\n"
            + f"Degrees of Freedom: {df}"
            + "\n"
            + f"Standardized Standard Error: {round(sse, 3)}"
            + "\n"
            + f"Type I Error Rate: {round(alpha, 2)}"
            + "\n"
            + f"Type II Error Rate: {round(1 - power, 2)}"
            + "\n"
            + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return mde
