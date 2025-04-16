from math import sqrt, ceil

from pyPowerUp.utils import _power


def power_bcra3f2(
        rho2: float,
        n: float,
        J: float,
        K: float,
        effect_size: float = 0.25,
        alpha: float = 0.10,
        two_tailed: bool = True,
        p: float = 0.50,
        g2: int = 0,
        r21: float = 0,
        r22: float = 0,
        print_pretty: bool = True,
) -> float:
    """Calculates the power of a Three-Level Blocked (Fixed) Cluster-level Random Assignment Design,
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
    effect_size : float, default=0.25
        Effect size
    alpha : float, default=0.10
        Probability of Type I error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.50
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
    The power of the test
    """

    df = ceil(K * (J - 2) - g2)
    sse = sqrt(
        rho2 * (1 - r22) / (p * (1 - p) * J * K)
        + (1 - rho2) * (1 - r21) / (p * (1 - p) * J * K * n)
    )
    power = _power(effect_size, alpha, sse, df, two_tailed)
    if print_pretty:
        str_print = (
                "Statistical Power"
                + "\n"
                + "-" * 39
                + "\n"
                + f" {round(power, 3)}"
                + "\n"
                + "-" * 39
                + "\n"
                + f"Degrees of Freedom: {df}"
                + "\n"
                + f"Standardized Standard Error: {round(sse, 3)}"
                + "\n"
                + f"Type I Error Rate: {round(alpha, 2)}"
                + "\n"
                + f"Type II Error Rate: {round(1 - power, 3)}"
                + "\n"
                + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return power


def power_bcra3r2(
        rho2: float,
        rho3: float,
        omega3: float,
        n: float,
        J: float,
        K: int,
        effect_size: float = 0.25,
        alpha: float = 0.10,
        two_tailed: bool = True,
        p: float = 0.50,
        g3: int = 0,
        r21: float = 0,
        r22: float = 0,
        r2t3: float = 0,
        print_pretty: bool = True,
) -> float:
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
    K : int
        Number of level 3 units
    effect_size : float, default=0.25
        Effect size
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.50
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
    The power of our test
    """
    df = K - g3 - 1
    sse = sqrt(
        rho3 * omega3 * (1 - r2t3) / K
        + rho2 * (1 - r22) / (p * (1 - p) * J * K)
        + (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * n)
    )
    power = _power(effect_size, alpha, sse, df, two_tailed)
    if print_pretty:
        str_print = (
                "Statistical Power"
                + "\n"
                + "-" * 39
                + "\n"
                + f" {round(power, 3)}"
                + "\n"
                + "-" * 39
                + "\n"
                + f"Degrees of Freedom: {df}"
                + "\n"
                + f"Standardized Standard Error: {round(sse, 3)}"
                + "\n"
                + f"Type I Error Rate: {round(alpha, 2)}"
                + "\n"
                + f"Type II Error Rate: {round(1 - power, 3)}"
                + "\n"
                + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return power


def power_bcra4f3(
        rho2: float,
        rho3: float,
        n: int,
        J: int,
        K: int,
        L: int,
        effect_size: float = 0.25,
        alpha: float = 0.10,
        two_tailed: bool = True,
        p: float = 0.50,
        r21: float = 0,
        r22: float = 0,
        r23: float = 0,
        g3: int = 0,
        print_pretty: bool = True,
) -> float:
    """Calculates the power of a Four-Level Blocked (Fixed) Cluster-level Random Assignment Design,
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
    effect_size : float, default=0.25
        Effect size
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.50
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
    The power of our test
    """
    df = L * (K - 2) - g3
    sse = sqrt(
        rho3 * (1 - r23) / (p * (1 - p) * K * L)
        + rho2 * (1 - r22) / (p * (1 - p) * J * K * L)
        + (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * L * n)
    )
    power = _power(effect_size, alpha, sse, df, two_tailed)
    if print_pretty:
        str_print = (
                "Statistical Power"
                + "\n"
                + "-" * 39
                + "\n"
                + f" {round(power, 3)}"
                + "\n"
                + "-" * 39
                + "\n"
                + f"Degrees of Freedom: {df}"
                + "\n"
                + f"Standardized Standard Error: {round(sse, 3)}"
                + "\n"
                + f"Type I Error Rate: {round(alpha, 2)}"
                + "\n"
                + f"Type II Error Rate: {round(1 - power, 3)}"
                + "\n"
                + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return power


def power_bcra4r2(
        rho2: float,
        rho3: float,
        rho4: float,
        omega3: float,
        omega4: float,
        n: int,
        J: int,
        K: int,
        L: int,
        effect_size: float = 0.25,
        alpha: float = 0.10,
        two_tailed: bool = True,
        p: float = 0.50,
        r21: float = 0,
        r22: float = 0,
        r2t3: float = 0,
        r2t4: float = 0,
        g4: int = 0,
        print_pretty: bool = True,
) -> float:
    """Calculates the power of a Four-Level Blocked Cluster-level Random Assignment Design, Treatment
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
    effect_size : float, default=0.25
        Effect size
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.50
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
    The power of our test
    """
    df = L - g4 - 1
    sse = sqrt(
        rho4 * omega4 * (1 - r2t4) / L
        + rho3 * omega3 * (1 - r2t3) / (K * L)
        + rho2 * (1 - r22) / (p * (1 - p) * J * K * L)
        + (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * L * n)
    )
    power = _power(effect_size, alpha, sse, df, two_tailed)
    if print_pretty:
        str_print = (
                "Statistical Power"
                + "\n"
                + "-" * 39
                + "\n"
                + f" {round(power, 3)}"
                + "\n"
                + "-" * 39
                + "\n"
                + f"Degrees of Freedom: {df}"
                + "\n"
                + f"Standardized Standard Error: {round(sse, 3)}"
                + "\n"
                + f"Type I Error Rate: {round(alpha, 2)}"
                + "\n"
                + f"Type II Error Rate: {round(1 - power, 3)}"
                + "\n"
                + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return power


def power_bcra4r3(
        rho2: float,
        rho3: float,
        rho4: float,
        omega4: float,
        n: float,
        J: float,
        K: int,
        L: int,
        effect_size: float = 0.25,
        alpha: float = 0.10,
        two_tailed: bool = True,
        p: float = 0.50,
        r21: float = 0,
        r22: float = 0,
        r23: float = 0,
        r2t4: float = 0,
        g4: int = 0,
        print_pretty: bool = True,
) -> float:
    """Calculates the power of a Four-Level Blocked Cluster-level Random Assignment Design, Treatment
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
    effect_size : float, default=0.25
        Effect size
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.50
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
    The power of our test
    """
    df = L - g4 - 1
    sse = sqrt(
        rho4 * omega4 * (1 - r2t4) / L
        + rho3 * (1 - r23) / (p * (1 - p) * K * L)
        + rho2 * (1 - r22) / (p * (1 - p) * J * K * L)
        + (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * L * n)
    )
    power = _power(effect_size, alpha, sse, df, two_tailed)
    if print_pretty:
        str_print = (
                "Statistical Power"
                + "\n"
                + "-" * 39
                + "\n"
                + f" {round(power, 3)}"
                + "\n"
                + "-" * 39
                + "\n"
                + f"Degrees of Freedom: {df}"
                + "\n"
                + f"Standardized Standard Error: {round(sse, 3)}"
                + "\n"
                + f"Type I Error Rate: {round(alpha, 2)}"
                + "\n"
                + f"Type II Error Rate: {round(1 - power, 3)}"
                + "\n"
                + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return power


def power_bira2c1(
        n: float,
        J: float,
        effect_size: float = 0.25,
        alpha: float = 0.10,
        two_tailed: bool = True,
        p: float = 0.50,
        g1: int = 0,
        r21: float = 0,
        print_pretty: bool = True,
) -> float:
    """Calculates the power of a Two-Level Blocked (Constant Treatment Effect) Individual-level
    Random Assignment Design, Treatment at Level 1

    Parameters
    ----------
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    effect_size : float, default=0.25
        Statistical power
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.50
        Average proportion of level 1 units randomly assigned to treatment within level 2 units
    g1 : int, default=0
        Number of covariates at level 1
    r21 : float, default=0
        Proportion of level 1 variance in the outcome explained by level 1 covariates
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    The power of our test
    """
    df = ceil(J * (n - 1) - g1 - 1)
    sse = sqrt((1 - r21) / (p * (1 - p) * J * n))
    power = _power(effect_size, alpha, sse, df, two_tailed)
    if print_pretty:
        str_print = (
                "Statistical Power"
                + "\n"
                + "-" * 39
                + "\n"
                + f" {round(power, 3)}"
                + "\n"
                + "-" * 39
                + "\n"
                + f"Degrees of Freedom: {df}"
                + "\n"
                + f"Standardized Standard Error: {round(sse, 3)}"
                + "\n"
                + f"Type I Error Rate: {round(alpha, 2)}"
                + "\n"
                + f"Type II Error Rate: {round(1 - power, 3)}"
                + "\n"
                + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return power


def power_bira2f1(
        n: int,
        J: int,
        effect_size: float = 0.25,
        alpha: float = 0.10,
        two_tailed: bool = True,
        p: float = 0.50,
        g1: int = 0,
        r21: float = 0,
        print_pretty: bool = True,
) -> float:
    """Calculates the power of a Two-Level Blocked (Fixed) Individual-level Random Assignment Design,
    Treatment at Level 1

    Parameters
    ----------
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    effect_size : float, default=0.25
        Effect size
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.50
        Average proportion of level 1 units randomly assigned to treatment within level 2 units
    g1 : int, default=0
        Number of covariates at level 1
    r21 : float, default=0
        Proportion of level 1 variance in the outcome explained by level 1 covariates
    print_pretty : bool, default=True
        Whether we wish to print the results similar to PowerUpR's output

    Returns
    -------
    The power of our test
    """
    df = J * (n - 2) - g1
    sse = sqrt((1 - r21) / (p * (1 - p) * J * n))
    power = _power(effect_size, alpha, sse, df, two_tailed)
    if print_pretty:
        str_print = (
                "Statistical Power"
                + "\n"
                + "-" * 39
                + "\n"
                + f" {round(power, 3)}"
                + "\n"
                + "-" * 39
                + "\n"
                + f"Degrees of Freedom: {df}"
                + "\n"
                + f"Standardized Standard Error: {round(sse, 3)}"
                + "\n"
                + f"Type I Error Rate: {round(alpha, 2)}"
                + "\n"
                + f"Type II Error Rate: {round(1 - power, 3)}"
                + "\n"
                + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return power


def power_bira2r1(
        rho2: float,
        omega2: float,
        n: float,
        J: float,
        effect_size: float = 0.25,
        alpha: float = 0.10,
        two_tailed: bool = True,
        p: float = 0.50,
        g2: int = 0,
        r21: float = 0,
        r2t2: float = 0,
        print_pretty: bool = True,
) -> float:
    """Calculates the power of a Two-Level Blocked Individual-level Random Assignment Design,
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
    effect_size : float, default=0.25
        Effecet size
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.50
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
    The power of our test
    """
    df = ceil(J - g2 - 1)
    sse = sqrt(
        rho2 * omega2 * (1 - r2t2) / J + (1 - rho2) * (1 - r21) / (p * (1 - p) * J * n)
    )
    power = _power(effect_size, alpha, sse, df, two_tailed)
    if print_pretty:
        str_print = (
                "Statistical Power"
                + "\n"
                + "-" * 39
                + "\n"
                + f" {round(power, 3)}"
                + "\n"
                + "-" * 39
                + "\n"
                + f"Degrees of Freedom: {df}"
                + "\n"
                + f"Standardized Standard Error: {round(sse, 3)}"
                + "\n"
                + f"Type I Error Rate: {round(alpha, 2)}"
                + "\n"
                + f"Type II Error Rate: {round(1 - power, 3)}"
                + "\n"
                + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return power


def power_bira3r1(
        rho2: float,
        rho3: float,
        omega2: float,
        omega3: float,
        n: float,
        J: float,
        K: int,
        effect_size: float = 0.25,
        alpha: float = 0.10,
        two_tailed: bool = True,
        p: float = 0.50,
        r21: float = 0,
        r2t2: float = 0,
        r2t3: float = 0,
        g3: int = 0,
        print_pretty: bool = True,
) -> float:
    """Calculates the power of a Three-Level Blocked Individual-level Random Assignment Design, Treatment at Level 1

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
    effect_size : float, default=0.25
        Effect size
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.50
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
    The power of our test
    """
    df = K - g3 - 1
    sse = sqrt(
        rho3 * omega3 * (1 - r2t3) / K
        + rho2 * omega2 * (1 - r2t2) / (J * K)
        + (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * n)
    )
    power = _power(effect_size, alpha, sse, df, two_tailed)
    if print_pretty:
        str_print = (
                "Statistical Power"
                + "\n"
                + "-" * 39
                + "\n"
                + f" {round(power, 3)}"
                + "\n"
                + "-" * 39
                + "\n"
                + f"Degrees of Freedom: {df}"
                + "\n"
                + f"Standardized Standard Error: {round(sse, 3)}"
                + "\n"
                + f"Type I Error Rate: {round(alpha, 2)}"
                + "\n"
                + f"Type II Error Rate: {round(1 - power, 3)}"
                + "\n"
                + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return power


def power_bira4r1(
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
        effect_size: float = 0.25,
        alpha: float = 0.10,
        two_tailed: bool = True,
        p: float = 0.50,
        r21: float = 0,
        r2t2: float = 0,
        r2t3: float = 0,
        r2t4: float = 0,
        g4: int = 0,
        print_pretty: bool = True,
) -> float:
    """Calculates the power of a Four-Level Blocked Individual-level Random Assignment Design,
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
    effect_size : float, default=0.25
        Effect size
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
    p : float, default=0.50
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
    The power of our test
    """
    df = L - g4 - 1
    sse = sqrt(
        rho4 * omega4 * (1 - r2t4) / L
        + rho3 * omega3 * (1 - r2t3) / (K * L)
        + rho2 * omega2 * (1 - r2t2) / (J * K * L)
        + (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * L * n)
    )
    power = _power(effect_size, alpha, sse, df, two_tailed)
    if print_pretty:
        str_print = (
                "Statistical Power"
                + "\n"
                + "-" * 39
                + "\n"
                + f" {round(power, 3)}"
                + "\n"
                + "-" * 39
                + "\n"
                + f"Degrees of Freedom: {df}"
                + "\n"
                + f"Standardized Standard Error: {round(sse, 3)}"
                + "\n"
                + f"Type I Error Rate: {round(alpha, 2)}"
                + "\n"
                + f"Type II Error Rate: {round(1 - power, 3)}"
                + "\n"
                + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return power


def power_cra2r2(
        rho2: float,
        n: float,
        J: float,
        effect_size: float = 0.25,
        alpha: float = 0.10,
        two_tailed: bool = True,
        p: float = 0.50,
        g2: int = 0,
        r21: float = 0,
        r22: float = 0,
        print_pretty: bool = True,
) -> float:
    """Calculates the power of a Two-level Cluster-randomized Trials to Detect Main, Moderation and
    Mediation Effects

    Parameters
    ----------
    rho2 : float
        Proportion of variance in the outcome between level 2 units (unconditional ICC2)
    n : float
        Harmonic mean of level 1 units across level 2 units (or simple average)
    J : float
        Harmonic mean of level 2 units across level 3 units (or simple average)
    effect_size : float, default=0.25
        Effect size
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
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
    The power of our test
    """
    df = ceil(J - g2 - 2)
    sse = sqrt(
        rho2 * (1 - r22) / (p * (1 - p) * J)
        + (1 - rho2) * (1 - r21) / (p * (1 - p) * J * n)
    )
    power = _power(effect_size, alpha, sse, df, two_tailed)
    if print_pretty:
        str_print = (
                "Statistical Power"
                + "\n"
                + "-" * 39
                + "\n"
                + f" {round(power, 3)}"
                + "\n"
                + "-" * 39
                + "\n"
                + f"Degrees of Freedom: {df}"
                + "\n"
                + f"Standardized Standard Error: {round(sse, 3)}"
                + "\n"
                + f"Type I Error Rate: {round(alpha, 2)}"
                + "\n"
                + f"Type II Error Rate: {round(1 - power, 3)}"
                + "\n"
                + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return power


def power_cra3r3(
        rho2: float,
        rho3: float,
        n: float,
        J: float,
        K: int,
        effect_size: float = 0.25,
        alpha: float = 0.10,
        two_tailed: bool = True,
        p: float = 0.50,
        g3: int = 0,
        r21: float = 0,
        r22: float = 0,
        r23: float = 0,
        print_pretty: bool = True,
) -> float:
    """Calculates the power of a Three-level Cluster-randomized Trials to Detect Main, Moderation,
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
    effect_size : float, default=0.25
        Effect size
    alpha : float, default=0.10
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
    The power of our test
    """
    df = K - g3 - 2
    sse = sqrt(
        rho3 * (1 - r23) / (p * (1 - p) * K)
        + rho2 * (1 - r22) / (p * (1 - p) * J * K)
        + (1 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * n)
    )
    power = _power(effect_size, alpha, sse, df, two_tailed)
    if print_pretty:
        str_print = (
                "Statistical Power"
                + "\n"
                + "-" * 39
                + "\n"
                + f" {round(power, 3)}"
                + "\n"
                + "-" * 39
                + "\n"
                + f"Degrees of Freedom: {df}"
                + "\n"
                + f"Standardized Standard Error: {round(sse, 3)}"
                + "\n"
                + f"Type I Error Rate: {round(alpha, 2)}"
                + "\n"
                + f"Type II Error Rate: {round(1 - power, 3)}"
                + "\n"
                + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return power


def power_cra4r4(
        rho2: float,
        rho3: float,
        rho4: float,
        n: float,
        J: float,
        K: float,
        L: int,
        effect_size: float = 0.25,
        alpha: float = 0.10,
        two_tailed: bool = True,
        p: float = 0.50,
        r21: float = 0,
        r22: float = 0,
        r23: float = 0,
        r24: float = 0,
        g4: int = 0,
        print_pretty: bool = True,
) -> float:
    """Calculates the power of a Four-Level Cluster-randomized Trial

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
    effect_size : float, default=0.25
        Effect size
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
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
    The power of our test
    """
    df = L - g4 - 2
    sse = sqrt(
        rho4 * (1 - r24) / (p * (1 - p) * L)
        + rho3 * (1 - r23) / (p * (1 - p) * K * L)
        + rho2 * (1 - r22) / (p * (1 - p) * J * K * L)
        + (1 - rho4 - rho3 - rho2) * (1 - r21) / (p * (1 - p) * J * K * L * n)
    )
    power = _power(effect_size, alpha, sse, df, two_tailed)
    if print_pretty:
        str_print = (
                "Statistical Power"
                + "\n"
                + "-" * 39
                + "\n"
                + f" {round(power, 3)}"
                + "\n"
                + "-" * 39
                + "\n"
                + f"Degrees of Freedom: {df}"
                + "\n"
                + f"Standardized Standard Error: {round(sse, 3)}"
                + "\n"
                + f"Type I Error Rate: {round(alpha, 2)}"
                + "\n"
                + f"Type II Error Rate: {round(1 - power, 3)}"
                + "\n"
                + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return power


def power_ira1r1(
        n: int,
        effect_size: float = 0.25,
        alpha: float = 0.10,
        two_tailed: bool = True,
        p: float = 0.50,
        g1: int = 0,
        r21: float = 0,
        print_pretty: bool = True,
) -> float:
    """Calculates the power of a Individual-level Random Assignment Design

    Parameters
    ----------
    n : int
        Sample size
    effect_size : float, default=0.25
        Effect size
    alpha : float, default=0.10
        Probability of Type 1 error
    two_tailed : bool, default=True
        Whether our hypothesis is one tailed or two tailed
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
    The power of our test
    """
    df = n - g1 - 2
    sse = sqrt((1 - r21) / (p * (1 - p) * n))
    power = _power(effect_size, alpha, sse, df, two_tailed)
    if print_pretty:
        str_print = (
                "Statistical Power"
                + "\n"
                + "-" * 39
                + "\n"
                + f" {round(power, 3)}"
                + "\n"
                + "-" * 39
                + "\n"
                + f"Degrees of Freedom: {df}"
                + "\n"
                + f"Standardized Standard Error: {round(sse, 3)}"
                + "\n"
                + f"Type I Error Rate: {round(alpha, 2)}"
                + "\n"
                + f"Type II Error Rate: {round(1 - power, 3)}"
                + "\n"
                + f"Two-Tailed Test: {two_tailed}"
        )
        print(str_print)
    return power
