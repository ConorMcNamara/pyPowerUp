import pytest

from pyPowerUp import power


def test_bcra3f2() -> None:
    result = power.power_bcra3f2(effect_size=.145, rho2=.10, n=20, J=44, K=5, alpha=0.05)
    # power.bcra3f2(es = .145, rho2=.10, n=20, J=44, K=5)
    #
    # Statistical power:
    # ---------------------------------------
    #  0.803
    # ---------------------------------------
    # Degrees of freedom: 210
    # Standardized standard error: 0.051
    # Type I error rate: 0.05
    # Type II error rate: 0.197
    # Two-tailed test: TRUE
    assert result == pytest.approx(0.803, abs=0.001)


def test_bcra3r2() -> None:
    result = power.power_bcra3r2(effect_size=.246, rho3=.13, rho2=.10, omega3=.4, n=10, J=6, K=24, alpha=0.05)
    # power.bcra3r2(es = .246, rho3=.13, rho2=.10, omega3=.4, n=10, J=6, K=24)
    #
    # Statistical power:
    # ---------------------------------------
    #  0.799
    # ---------------------------------------
    # Degrees of freedom: 23
    # Standardized standard error: 0.084
    # Type I error rate: 0.05
    # Type II error rate: 0.201
    # Two-tailed test: TRUE
    assert result == pytest.approx(0.799, abs=0.001)


def test_bcra4f3() -> None:
    result = power.power_bcra4f3(effect_size=0.339, rho3=.15, rho2=.15, n=10, J=4, K=4, L=15, alpha=0.05)
    # power.bcra4f3(es=0.339, rho3=.15, rho2=.15, n=10, J=4, K=4, L=15)
    #
    # Statistical power:
    # ---------------------------------------
    #  0.801
    # ---------------------------------------
    # Degrees of freedom: 30
    # Standardized standard error: 0.117
    # Type I error rate: 0.05
    # Type II error rate: 0.199
    # Two-tailed test: TRUE
    assert result == pytest.approx(0.801, abs=0.001)


def test_bcra4r2() -> None:
    result = power.power_bcra4r2(effect_size=.206, rho4=.05, rho3=.15, rho2=.15, omega4=.50, omega3=.50, n=10, J=4, K=4,
                                 L=20, alpha=0.05)
    # power.bcra4r2(es = .206, rho4=.05, rho3=.15, rho2=.15, omega4=.50, omega3=.50, n=10, J=4, K=4, L=20)
    #
    # Statistical power:
    # ---------------------------------------
    #  0.799
    # ---------------------------------------
    # Degrees of freedom: 19
    # Standardized standard error: 0.07
    # Type I error rate: 0.05
    # Type II error rate: 0.201
    # Two-tailed test: TRUE
    assert result == pytest.approx(0.799, abs=0.001)


def test_bcra4r3() -> None:
    result = power.power_bcra4r3(effect_size=.316, rho4=.05, rho3=.15, rho2=.15, omega4=.50, n=10, J=4, K=4, L=20,
                                 alpha=0.05)
    # power.bcra4r3(es = .316, rho4=.05, rho3=.15, rho2=.15, omega4=.50, n=10, J=4, K=4, L=20)
    #
    # Statistical power:
    # ---------------------------------------
    #  0.8
    # ---------------------------------------
    # Degrees of freedom: 19
    # Standardized standard error: 0.107
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result == pytest.approx(0.800, abs=0.001)


def test_bira2c1() -> None:
    result = power.power_bira2c1(effect_size=.325, n=15, J=20, alpha=0.05)
    # power.bira2c1(es=.325, n=15, J=20, alpha=0.05)
    #
    # Statistical power:
    # ---------------------------------------
    #  0.801
    # ---------------------------------------
    # Degrees of freedom: 279
    # Standardized standard error: 0.115
    # Type I error rate: 0.05
    # Type II error rate: 0.199
    # Two-tailed test: TRUE
    assert result == pytest.approx(0.801, abs=0.001)


def test_bira2f1() -> None:
    result = power.power_bira2f1(effect_size=.325, n=15, J=20, alpha=0.05)
    # power.bira2f1(es=.325, n=15, J=20)
    #
    # Statistical power:
    # ---------------------------------------
    #  0.801
    # ---------------------------------------
    # Degrees of freedom: 260
    # Standardized standard error: 0.115
    # Type I error rate: 0.05
    # Type II error rate: 0.199
    # Two-tailed test: TRUE
    assert result == pytest.approx(0.801, abs=0.001)


def test_bira2r1() -> None:
    result = power.power_bira2r1(effect_size=.366, rho2=.17, omega2=.50, n=15, J=20, alpha=0.05)
    # power.bira2r1(es=.366, rho2=.17, omega2=.50, n=15, J=20)
    #
    # Statistical power:
    # ---------------------------------------
    #  0.801
    # ---------------------------------------
    # Degrees of freedom: 19
    # Standardized standard error: 0.124
    # Type I error rate: 0.05
    # Type II error rate: 0.199
    # Two-tailed test: TRUE
    assert result == pytest.approx(0.801, abs=0.001)


def test_bira3r1() -> None:
    result = power.power_bira3r1(effect_size=.045, rho3=.20, rho2=.15, omega3=.10, omega2=.10, n=69, J=10, K=100,
                                 alpha=0.05)
    # power.bira3r1(es = .045, rho3=.20, rho2=.15, omega3=.10, omega2=.10, n=69, J=10, K=100)
    #
    # Statistical power:
    # ---------------------------------------
    #  0.8
    # ---------------------------------------
    # Degrees of freedom: 99
    # Standardized standard error: 0.016
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result == pytest.approx(0.800, abs=0.001)


def test_bira4r1() -> None:
    result = power.power_bira4r1(effect_size=0.142, rho4=.05, rho3=.15, rho2=.15, omega4=.50, omega3=.50, omega2=.50,
                                 n=10, J=4, K=4, L=27, alpha=0.05)
    # power.bira4r1(es = 0.142, rho4=.05, rho3=.15, rho2=.15, omega4=.50, omega3=.50, omega2=.50,n=10, J=4, K=4, L=27)
    #
    # Statistical power:
    # ---------------------------------------
    #  0.797
    # ---------------------------------------
    # Degrees of freedom: 26
    # Standardized standard error: 0.049
    # Type I error rate: 0.05
    # Type II error rate: 0.203
    # Two-tailed test: TRUE
    assert result == pytest.approx(0.797, abs=0.001)


def test_cra2r2() -> None:
    result = power.power_cra2r2(effect_size=.629, rho2=.17, n=15, J=20, alpha=0.05)
    # power.cra2r2(es=.629, rho2=.17, n=15, J=20)
    #
    # Statistical power:
    # ---------------------------------------
    #  0.8
    # ---------------------------------------
    # Degrees of freedom: 18
    # Standardized standard error: 0.212
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result == pytest.approx(0.800, abs=0.001)


def test_cra3r3() -> None:
    result = power.power_cra3r3(effect_size=.269, rho3=.06, rho2=.17, n=15, J=3, K=60, alpha=0.05)
    # power.cra3r3(es=.269, rho3=.06, rho2=.17, n=15, J=3, K=60)
    #
    # Statistical power:
    # ---------------------------------------
    #  0.8
    # ---------------------------------------
    # Degrees of freedom: 58
    # Standardized standard error: 0.094
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result == pytest.approx(0.800, abs=0.001)


def test_ira1r1() -> None:
    result = power.power_ira1r1(effect_size=0.356, n=250, alpha=0.05)
    # power.ira1r1(es=.356, n=250)
    #
    # Statistical power:
    # ---------------------------------------
    #  0.801
    # ---------------------------------------
    # Degrees of freedom: 248
    # Standardized standard error: 0.126
    # Type I error rate: 0.05
    # Type II error rate: 0.199
    # Two-tailed test: TRUE
    assert result == pytest.approx(0.801, abs=0.001)
