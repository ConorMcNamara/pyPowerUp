from math import ceil
from typing import Optional

from scipy.stats import t as t_dist

import numpy as np


def sample_size_bcra3f2(
        rho2: float,
        n: float,
        J: float,
        effect_size: float = 0.25,
        power: float = 0.80,
        alpha: float = 0.10,
        two_tailed: bool = True,
        K0: int = 10,
        tol: float = 0.10,
        p: float = 0.50,
        g2: int = 0,
        r21: float = 0,
        r22: float = 0,
        print_pretty: bool = True,
) -> Optional[int]:
    """Calculates the minimum sample size of a Three-Level Blocked (Fixed) Cluster-level Random Assignment Design,
    Treatment at Level 2

    Parameters
    ----------
    rho2 : float
        Proportion of variance in the outcome between level 2 units (unconditional ICC2)
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    effect_size : float, default=0.25
        The effecet size
    power : float, default=0.80
        The power of our test
    alpha : float, default=0.10
        Probability of Type I error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    K0 : int, default=10
        Our initial guess for the minimum sample size
    tol : float, default=0.10
        Tolerance to end iterative process for finding the minimum sample size
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
    The minimum sample size of our test
    """
    df = 0
    for i in range(100):
        df = K0 * (J - 2) - g2
        if df < 0 or np.isinf(df):
            break
        T1 = (
            abs(t_dist.ppf(alpha / 2, df)) if two_tailed else abs(t_dist.ppf(alpha, df))
        )
        T2 = abs(t_dist.ppf(power, df))
        M = T1 + T2 if power >= 0.5 else T1 - T2
        K1 = pow(M / effect_size, 2) * (
                rho2 * (1 - r22) / (p * (1 - p) * J)
                + (1 - rho2) * (1 - r21) / (p * (1 - p) * J * n)
        )
        if abs(K1 - K0) < tol:
            break
        K0 = (K1 + K0) / 2
    sample_size = ceil(K0) if df > 0 else None
    if print_pretty:
        print(f"K = {sample_size}")
    return sample_size


def sample_size_bcra3r2(
        rho2: float,
        rho3: float,
        omega3: float,
        n: float,
        J: float,
        effect_size: float = 0.25,
        power: float = 0.80,
        alpha: float = 0.10,
        two_tailed: bool = True,
        K0: int = 10,
        tol: float = 0.10,
        p: float = 0.50,
        g3: int = 0,
        r21: float = 0,
        r22: float = 0,
        r2t3: float = 0,
        print_pretty: bool = True,
) -> Optional[int]:
    """Calculates the power of a Three-Level Blocked Cluster-level Random Assignment Design,
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
    effect_size : float, default=0.25
        Effect size
    power : float, default=0.80
        The power of our test
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    K0 : int, default=10
        Our initial guess for the minimum sample size
    tol : float, default=0.10
        Tolerance to end iterative process for finding the minimum sample size
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
    The minimum sample size of our test
    """
    df = 0
    for i in range(100):
        df = K0 - g3 - 1
        if df < 0 or np.isinf(df):
            break
        T1 = (
            abs(t_dist.ppf(alpha / 2, df)) if two_tailed else abs(t_dist.ppf(alpha, df))
        )
        T2 = abs(t_dist.ppf(power, df))
        M = T1 + T2 if power >= 0.5 else T1 - T2
        K1 = pow(M / effect_size, 2) * (
                rho3 * omega3 * (1 - r2t3)
                + rho2 * (1 - r22) / (p * (1 - p) * J)
                + (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * n)
        )
        if abs(K1 - K0) < tol:
            break
        K0 = (K1 + K0) / 2
    sample_size = ceil(K0) if df > 0 else None
    if print_pretty:
        print(f"K = {sample_size}")
    return sample_size


def sample_size_bcra4f3(
        rho2: float,
        rho3: float,
        n: float,
        J: float,
        K: float,
        effect_size: float = 0.25,
        power: float = 0.80,
        alpha: float = 0.05,
        two_tailed: bool = True,
        L0: int = 10,
        tol: float = 0.10,
        p: float = 0.50,
        g3: int = 0,
        r21: float = 0,
        r22: float = 0,
        r23: float = 0,
        print_pretty: bool = True,
) -> Optional[int]:
    """Calculates the minimum sample size of a Four-Level Blocked (Fixed) Cluster-level Random Assignment Design,
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
    effect_size : float, default=0.25
        Effect size
    power : float, default=0.80
        The power of our test
    alpha : float, default=0.1
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    L0 : int, default=10
        Our initial guess of our minimum sample size
    tol : float, default=0.10
        Tolerance to end iterative process for finding the minimum sample size
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
    The minimum sample size of our test
    """
    df = 0
    for i in range(100):
        df = L0 * (K - 2) - g3
        if df < 0 or np.isinf(df):
            break
        T1 = (
            abs(t_dist.ppf(alpha / 2, df)) if two_tailed else abs(t_dist.ppf(alpha, df))
        )
        T2 = abs(t_dist.ppf(power, df))
        M = T1 + T2 if power >= 0.5 else T1 - T2
        L1 = pow(M / effect_size, 2) * (
                rho3 * (1 - r23) / (p * (1 - p) * K)
                + rho2 * (1 - r22) / (p * (1 - p) * K * J)
                + (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * n)
        )
        if abs(L1 - L0) < tol:
            break
        L0 = (L1 + L0) / 2
    sample_size = ceil(L0) if df > 0 else None
    if print_pretty:
        print(f"L = {sample_size}")
    return sample_size


def sample_size_bcra4r2(
        rho2: float,
        rho3: float,
        rho4: float,
        omega3: float,
        omega4: float,
        n: float,
        J: float,
        K: float,
        effect_size: float = 0.25,
        power: float = 0.80,
        alpha: float = 0.10,
        two_tailed: bool = True,
        L0: int = 10,
        tol: float = 0.10,
        p: float = 0.50,
        r21: float = 0,
        r22: float = 0,
        r2t3: float = 0,
        r2t4: float = 0,
        g4: int = 0,
        print_pretty: bool = True,
) -> Optional[int]:
    """Calculates the minimum sample size of a Four-Level Blocked Cluster-level Random Assignment Design, Treatment
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
    effect_size : float, default=0.25
        Effect size
    power : float, default=0.80
        The power of our test
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    L0 : int, default=10
        Our initial guess of our minimum sample size
    tol : float, default=0.10
        Tolerance to end iterative process for finding the minimum sample size
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
    g4 : int, default=0
        Number of covariates at level 4
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    The minimum sample size of our test
    """
    df = 0
    for i in range(100):
        df = L0 - g4 - 1
        if df < 0 or np.isinf(df):
            break
        T1 = (
            abs(t_dist.ppf(alpha / 2, df)) if two_tailed else abs(t_dist.ppf(alpha, df))
        )
        T2 = abs(t_dist.ppf(power, df))
        M = T1 + T2 if power >= 0.5 else T1 - T2
        L1 = pow(M / effect_size, 2) * (
                rho4 * omega4 * (1 - r2t4)
                + rho3 * omega3 * (1 - r2t3) / K
                + rho2 * (1 - r22) / (p * (1 - p) * K * J)
                + (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * n)
        )
        if abs(L1 - L0) < tol:
            break
        L0 = (L1 + L0) / 2
    sample_size = ceil(L0) if df > 0 else None
    if print_pretty:
        print(f"L = {sample_size}")
    return sample_size


def sample_size_bcra4r3(
        rho2: float,
        rho3: float,
        rho4: float,
        omega4: float,
        n: float,
        J: float,
        K: float,
        effect_size: float = 0.25,
        power: float = 0.80,
        alpha: float = 0.10,
        two_tailed: bool = True,
        L0: int = 10,
        tol: float = 0.10,
        p: float = 0.50,
        r21: float = 0,
        r22: float = 0,
        r23: float = 0,
        r2t4: float = 0,
        g4: int = 0,
        print_pretty: bool = True,
) -> Optional[int]:
    """Calculates the minimum sample size of a Four-Level Blocked Cluster-level Random Assignment Design, Treatment
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
    effect_size : float, default=0.25
        The effect size
    power : float, default=0.8
        The power of our test
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    L0 : int, default=10
        Our initial guess of our minimum sample size
    tol : float, default=0.10
        Tolerance to end iterative process for finding the minimum sample size
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
    The minimum sample size of our test
    """
    df = 0
    for i in range(100):
        df = L0 - g4 - 1
        if df < 0 or np.isinf(df):
            break
        T1 = (
            abs(t_dist.ppf(alpha / 2, df)) if two_tailed else abs(t_dist.ppf(alpha, df))
        )
        T2 = abs(t_dist.ppf(power, df))
        M = T1 + T2 if power >= 0.5 else T1 - T2
        L1 = pow(M / effect_size, 2) * (
                rho4 * omega4 * (1 - r2t4)
                + rho3 * (1 - r23) / (p * (1 - p) * K)
                + rho2 * (1 - r22) / (p * (1 - p) * K * J)
                + (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * K * J * n)
        )
        if abs(L1 - L0) < tol:
            break
        L0 = (L1 + L0) / 2
    sample_size = ceil(L0) if df > 0 else None
    return sample_size


def sample_size_bira2c1(
        n: float,
        effect_size: float = 0.25,
        power: float = 0.80,
        alpha: float = 0.05,
        two_tailed: bool = True,
        J0: int = 10,
        tol: float = 0.10,
        p: float = 0.50,
        g1: int = 0,
        r21: float = 0,
        print_pretty: bool = True,
) -> Optional[int]:
    """Calculates the minimum sample size of a Two-Level Blocked Individual-level Random Assignment Design,
    Treatment at Level 1

    Parameters
    ----------
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    effect_size : float, default=0.25
        The effect size
    power : float, default=0.80
        The power of our test
    alpha : float, default=0.1
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    J0 : int, default=10
        Our initial guess of the minimum sample size
    tol : float, default=0.10
        Tolerance to end iterative process for finding the minimum sample size
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
    The minimum sample size of our test
    """
    df = 0
    for i in range(100):
        df = J0 * (n - 1) - g1 - 1
        if df < 0 or np.isinf(df):
            break
        T1 = (
            abs(t_dist.ppf(alpha / 2, df)) if two_tailed else abs(t_dist.ppf(alpha, df))
        )
        T2 = abs(t_dist.ppf(power, df))
        M = T1 + T2 if power >= 0.5 else T1 - T2
        J1 = pow(M / effect_size, 2) * ((1 - r21) / (p * (1 - p) * n))
        if abs(J1 - J0) < tol:
            break
        J0 = (J1 + J0) / 2
    sample_size = ceil(J0) if df > 0 else None
    if print_pretty:
        print(f"J = {sample_size}")
    return sample_size


def sample_size_bira2f1(
        n: float,
        effect_size: float = 0.25,
        power: float = 0.80,
        alpha: float = 0.05,
        two_tailed: bool = True,
        J0: int = 10,
        tol: bool = 0.10,
        p: float = 0.50,
        g1: int = 0,
        r21: float = 0,
        print_pretty: bool = True,
) -> Optional[int]:
    """Calculates the minimum sample size of a Two-Level Blocked (Fixed) Individual-level Random Assignment Design,
    Treatment at Level 1

    Parameters
    ----------
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    effect_size : float, default=0.25
        The effect size
    power : float, default=0.80
        The power of our test
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    J0 : int, default=10
        Our initial guess of the minimum sample size
    tol : float, default=0.10
        Tolerance to end iterative process for finding the minimum sample size
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
    The minimum sample size of our test
    """
    df = 0
    for i in range(100):
        df = J0 * (n - 2) - g1
        if df < 0 or np.isinf(df):
            break
        T1 = (
            abs(t_dist.ppf(alpha / 2, df)) if two_tailed else abs(t_dist.ppf(alpha, df))
        )
        T2 = abs(t_dist.ppf(power, df))
        M = T1 + T2 if power >= 0.5 else T1 - T2
        J1 = pow(M / effect_size, 2) * ((1 - r21) / (p * (1 - p) * n))
        if abs(J1 - J0) < tol:
            break
        J0 = (J1 + J0) / 2
    sample_size = ceil(J0) if df > 0 else None
    print(f"J = {sample_size}")
    return sample_size


def sample_size_bira2r1(
        rho2: float,
        omega2: float,
        n: float,
        effect_size: float = 0.25,
        power: float = 0.80,
        alpha: float = 0.05,
        two_tailed: bool = True,
        J0: int = 10,
        tol: bool = 0.10,
        p: float = 0.50,
        g2: int = 0,
        r21: float = 0,
        r2t2: float = 0,
        print_pretty: bool = True,
) -> Optional[int]:
    """Calculates the minimum sample size of a Two-Level Individual-level Random Assignment Design,
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
    effect_size : float, default=0.25
        The effect size
    power : float, default=0.80
        The power of our test
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    J0 : int, default=10
        Our initial guess of the minimum sample size
    tol : float, default=0.10
        Tolerance to end iterative process for finding the minimum sample size
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
    The minimum sample size of our test
    """
    df = 0
    for i in range(100):
        df = J0 - g2 - 1
        if df < 0 or np.isinf(df):
            break
        T1 = (
            abs(t_dist.ppf(alpha / 2, df)) if two_tailed else abs(t_dist.ppf(alpha, df))
        )
        T2 = abs(t_dist.ppf(power, df))
        M = T1 + T2 if power >= 0.5 else T1 - T2
        J1 = pow(M / effect_size, 2) * (
                rho2 * omega2 * (1 - r2t2) + (1 - rho2) * (1 - r21) / (p * (1 - p) * n)
        )
        if abs(J1 - J0) < tol:
            break
        J0 = (J1 + J0) / 2
    sample_size = ceil(J0) if df > 0 else None
    if print_pretty:
        print(f"J = {sample_size}")
    return sample_size


def sample_size_bira3r1(
        rho2: float,
        rho3: float,
        omega2: float,
        omega3: float,
        n: float,
        J: float,
        effect_size: float = 0.25,
        power: float = 0.80,
        alpha: float = 0.10,
        two_tailed: bool = True,
        K0: int = 10,
        tol: float = 0.10,
        p: float = 0.50,
        r21: float = 0,
        r2t2: float = 0,
        r2t3: float = 0,
        g3: int = 0,
        print_pretty: bool = True,
) -> Optional[int]:
    """Calculates the minimum sample size of a Three-Level Blocked Individual-level Random Assignment Design,
    Treatment at Level 1

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
    effect_size : float, default=0.25
        The effect size
    power : float, default=0.8
        The power of our test
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    K0 : int, default=10
        Our initial guess to the minimum sample size
    tol : float, default=0.10
        Tolerance to end iterative process for finding the minimum sample size
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
    The minimum sample size of our test
    """
    df = 0
    for i in range(100):
        df = K0 - g3 - 1
        if df < 0 or np.isinf(df):
            break
        T1 = (
            abs(t_dist.ppf(alpha / 2, df)) if two_tailed else abs(t_dist.ppf(alpha, df))
        )
        T2 = abs(t_dist.ppf(power, df))
        M = T1 + T2 if power >= 0.5 else T1 - T2
        K1 = pow(M / effect_size, 2) * (
                rho3 * omega3 * (1 - r2t3)
                + rho2 * omega2 * (1 - r2t2) / J
                + (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * n)
        )
        if abs(K1 - K0) < tol:
            break
        K0 = (K1 + K0) / 2
    sample_size = ceil(K0) if df > 0 else None
    if print_pretty:
        print(f"K = {sample_size}")
    return sample_size


def sample_size_bira4r1(
        rho2: float,
        rho3: float,
        rho4: float,
        omega2: float,
        omega3: float,
        omega4: float,
        n: float,
        J: float,
        K: float,
        effect_size: float = 0.25,
        power: float = 0.80,
        alpha: float = 0.10,
        two_tailed: bool = True,
        L0: int = 10,
        tol: float = 0.10,
        p: float = 0.50,
        r21: float = 0,
        r2t2: float = 0,
        r2t3: float = 0,
        r2t4: float = 0,
        g4: int = 0,
        print_pretty: bool = True,
) -> Optional[float]:
    """Calculates the minimum sample size of a Four-Level Blocked Individual-level Random Assignment Design,
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
    effect_size : float, default=0.25
        Effect size
    power : float, default=0.80
        The power of our test
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    L0 : int, default=10
        Our initial guess to the minimum sample size
    tol : float, default=0.10
        Tolerance to end iterative process for finding the minimum sample size
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
    The minimum sample size of our test
    """
    df = 0
    for i in range(100):
        df = L0 - g4 - 1
        if df < 0 or np.isinf(df):
            break
        T1 = (
            abs(t_dist.ppf(alpha / 2, df)) if two_tailed else abs(t_dist.ppf(alpha, df))
        )
        T2 = abs(t_dist.ppf(power, df))
        M = T1 + T2 if power >= 0.5 else T1 - T2
        L1 = pow(M / effect_size, 2) * (
                rho4 * omega4 * (1 - r2t4)
                + rho3 * omega3 * (1 - r2t3) / K
                + rho2 * omega2 * (1 - r2t2) / (K * J)
                + (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * K * J * n)
        )
        if abs(L1 - L0) < tol:
            break
        L0 = (L1 + L0) / 2
    sample_size = ceil(L0) if df > 0 else None
    if print_pretty:
        print(f"L = {sample_size}")
    return sample_size


def sample_size_cra2r2(
        rho2: float,
        n: float,
        effect_size: float = 0.25,
        power: float = 0.80,
        alpha: float = 0.10,
        two_tailed: bool = True,
        J0: int = 10,
        tol: float = 0.10,
        p: float = 0.50,
        g2: int = 0,
        r21: float = 0,
        r22: float = 0,
        print_pretty: bool = True,
) -> Optional[float]:
    """Calculates the minimum sample size of a Two-level Cluster-randomized Trials to Detect Main, Moderation and
    Mediation Effects

    Parameters
    ----------
    rho2 : float
        Proportion of variance in the outcome between level 2 units (unconditional ICC2)
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    effect_size : float, default=0.25
        Effect size
    power : float, default=0.80
        The power of our test
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    J0 : int, default=10
        Our initial guess to the minimum sample size
    tol : float, default=0.10
        Tolerance to end iterative process for finding the minimum sample size
    p : float, default=0.50
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
    The minimum sample size of our test
    """
    df = 0
    for i in range(100):
        df = J0 - g2 - 2
        if df < 0 or np.isinf(df):
            break
        T1 = (
            abs(t_dist.ppf(alpha / 2, df)) if two_tailed else abs(t_dist.ppf(alpha, df))
        )
        T2 = abs(t_dist.ppf(power, df))
        M = T1 + T2 if power >= 0.5 else T1 - T2
        J1 = pow(M / effect_size, 2) * (
                rho2 * (1 - r22) / (p * (1 - p))
                + (1 - rho2) * (1 - r21) / (p * (1 - p) * n)
        )
        if abs(J1 - J0) < tol:
            break
        J0 = (J1 + J0) / 2
    sample_size = ceil(J0) if df > 0 else None
    if print_pretty:
        print(f"J = {sample_size}")
    return sample_size


def sample_size_cra3r3(
        rho2: float,
        rho3: float,
        n: float,
        J: float,
        effect_size: float = 0.25,
        power: float = 0.80,
        alpha: float = 0.10,
        two_tailed: bool = True,
        K0: int = 10,
        tol: float = 0.10,
        p: float = 0.50,
        g3: int = 0,
        r21: float = 0,
        r22: float = 0,
        r23: float = 0,
        print_pretty: bool = True,
) -> Optional[float]:
    """Calculates the minimum sample size of a Three-level Cluster-randomized Trials to Detect Main, Moderation,
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
    effect_size : float, default=0.25
        Effect size
    power : float, default=0.80
        The power of our test
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    K0 : int, default=10
        Our initial guess to the minimum sample size
    tol : float, default=0.10
        Tolerance to end iterative process for finding the minimum sample size
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
    The minimum sample size of our test
    """
    df = 0
    for i in range(100):
        df = K0 - g3 - 2
        if df < 0 or np.isinf(df):
            break
        T1 = (
            abs(t_dist.ppf(alpha / 2, df)) if two_tailed else abs(t_dist.ppf(alpha, df))
        )
        T2 = abs(t_dist.ppf(power, df))
        M = T1 + T2 if power >= 0.5 else T1 - T2
        K1 = pow(M / effect_size, 2) * (
                rho3 * (1 - r23) / (p * (1 - p))
                + rho2 * (1 - r22) / (p * (1 - p) * J)
                + (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * n)
        )
        if abs(K1 - K0) < tol:
            break
        K0 = (K1 + K0) / 2
    sample_size = ceil(K0) if df > 0 else None
    if print_pretty:
        print(f"K = {sample_size}")
    return sample_size


def sample_size_cra4r4(
        rho2: float,
        rho3: float,
        rho4: float,
        n: float,
        J: float,
        K: float,
        effect_size: float = 0.25,
        power: float = 0.80,
        alpha: float = 0.10,
        two_tailed: bool = True,
        L0: int = 10,
        tol: float = 0.10,
        p: float = 0.50,
        r21: float = 0,
        r22: float = 0,
        r23: float = 0,
        r24: float = 0,
        g4: int = 0,
        print_pretty: bool = True,
) -> Optional[float]:
    """Calculates the minimum sample size of a Four-Level Cluster-randomized Trial

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
    effect_size : float, default=0.25
        Effect size
    power : float, default=0.80
        The power of our test
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    L0 : int, default=10
        Our initial guess to the minimum sample size
    tol : float, default=0.10
        Tolerance to end iterative process for finding the minimum sample size
    p : float, default=0.50
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
    The minimum sample size of our test
    """
    df = 0
    for i in range(100):
        df = L0 - g4 - 2
        if df < 0 or np.isinf(df):
            break
        T1 = (
            abs(t_dist.ppf(alpha / 2, df)) if two_tailed else abs(t_dist.ppf(alpha, df))
        )
        T2 = abs(t_dist.ppf(power, df))
        M = T1 + T2 if power >= 0.5 else T1 - T2
        L1 = pow(M / effect_size, 2) * (
                rho4 * (1 - r24) / (p * (1 - p))
                + rho3 * (1 - r23) / (p * (1 - p) * K)
                + rho2 * (1 - r22) / (p * (1 - p) * K * J)
                + (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * n)
        )
        if abs(L1 - L0) < tol:
            break
        L0 = (L1 + L0) / 2
    sample_size = ceil(L0) if df > 0 else None
    if print_pretty:
        print(f"L = {sample_size}")
    return sample_size


def sample_size_ira1r1(
        effect_size: float = 0.25,
        power: float = 0.80,
        alpha: float = 0.10,
        two_tailed=True,
        n0: int = 10,
        tol: float = 0.10,
        p: float = 0.50,
        g1: int = 0,
        r21: float = 0,
        print_pretty: bool = True,
) -> Optional[float]:
    """Calculates the minimum sample size of a Individual-level Random Assignment Design

    Parameters
    ----------
    effect_size : float, default=0.25
        Effect size
    power : float, default=0.80
        The power of our test
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    n0 : int, default=10
        Our initial guess to the minimum sample size
    tol : float, default=0.10
        Tolerance to end iterative process for finding the minimum sample size
    p : float, default=0.50
        Proportion of units randomly assigned to treatment
    g1 : int, default=0
        Number of covariates
    r21 : float, default=0
        Proportion of variance in the outcome explained by covariates
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    The minimum sample size of our test
    """
    df = 0
    for i in range(100):
        df = n0 - g1 - 1
        if df < 0 or np.isinf(df):
            break
        T1 = (
            abs(t_dist.ppf(alpha / 2, df)) if two_tailed else abs(t_dist.ppf(alpha, df))
        )
        T2 = abs(t_dist.ppf(power, df))
        M = T1 + T2 if power >= 0.5 else T1 - T2
        n1 = pow(M / effect_size, 2) * ((1 - r21) / (p * (1 - p)))
        if abs(n1 - n0) < tol:
            break
        n0 = (n1 + n0) / 2
    sample_size = ceil(n0) if df > 0 else None
    if print_pretty:
        print(f"n = {sample_size}")
    return sample_size
