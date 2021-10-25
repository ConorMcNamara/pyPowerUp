import pytest

from pyPowerUp import mde


def test_bcra3f2() -> None:
    result = mde.mde_bcra3f2(rho2=.10, n=20, J=44, K=5, alpha=0.05)
    # mdes.bcra3f2(rho2=.10, n=20, J=44, K=5)
    #
    # Minimum detectable effect size:
    # ---------------------------------------
    #  0.145 95% CI [0.043,0.246]
    # ---------------------------------------
    # Degrees of freedom: 210
    # Standardized standard error: 0.051
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result['minimum_detectable_effect'] == pytest.approx(0.145, abs=0.001)
    assert result['95% Confidence Interval'][0] == pytest.approx(0.043, abs=0.001)
    assert result['95% Confidence Interval'][1] == pytest.approx(0.246, abs=0.001)


def test_bcra3r2() -> None:
    result = mde.mde_bcra3r2(rho3=.13, rho2=.10, omega3=.4, n=10, J=6, K=24, alpha=0.05)
    # mdes.bcra3r2(rho3=.13, rho2=.10, omega3=.4, n=10, J=6, K=24)
    #
    # Minimum detectable effect size:
    # ---------------------------------------
    #  0.246 95% CI [0.072,0.42]
    # ---------------------------------------
    # Degrees of freedom: 23
    # Standardized standard error: 0.084
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result['minimum_detectable_effect'] == pytest.approx(0.246, abs=0.001)
    assert result['95% Confidence Interval'][0] == pytest.approx(0.072, abs=0.001)
    assert result['95% Confidence Interval'][1] == pytest.approx(0.420, abs=0.001)


def test_bcra4f3() -> None:
    result = mde.mde_bcra4f3(rho3=.15, rho2=.15, n=10, J=4, K=4, L=15, alpha=0.05)
    # mdes.bcra4f3(rho3=.15, rho2=.15, n=10, J=4, K=4, L=15)
    #
    # Minimum detectable effect size:
    # ---------------------------------------
    #  0.339 95% CI [0.1,0.577]
    # ---------------------------------------
    # Degrees of freedom: 30
    # Standardized standard error: 0.117
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result['minimum_detectable_effect'] == pytest.approx(0.339, abs=0.001)
    assert result['95% Confidence Interval'][0] == pytest.approx(0.100, abs=0.001)
    assert result['95% Confidence Interval'][1] == pytest.approx(0.577, abs=0.001)


def test_bcra4r2() -> None:
    result = mde.mde_bcra4r2(rho4=.05, rho3=.15, rho2=.15, omega4=.50, omega3=.50, n=10, J=4, K=4, L=20, alpha=0.05)
    # mdes.bcra4r2(rho4=.05, rho3=.15, rho2=.15, omega4=.50, omega3=.50, n=10, J=4, K=4, L=20)
    #
    # Minimum detectable effect size:
    # ---------------------------------------
    #  0.206 95% CI [0.06,0.352]
    # ---------------------------------------
    # Degrees of freedom: 19
    # Standardized standard error: 0.07
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result['minimum_detectable_effect'] == pytest.approx(0.206, abs=0.001)
    assert result['95% Confidence Interval'][0] == pytest.approx(0.060, abs=0.001)
    assert result['95% Confidence Interval'][1] == pytest.approx(0.352, abs=0.001)


def test_bcra4r3() -> None:
    result = mde.mde_bcra4r3(rho4=.05, rho3=.15, rho2=.15, omega4=.50, n=10, J=4, K=4, L=20, alpha=0.05)
    # mdes.bcra4r3(rho4=.05, rho3=.15, rho2=.15, omega4=.50, n=10, J=4, K=4, L=20)
    #
    # Minimum detectable effect size:
    # ---------------------------------------
    #  0.316 95% CI [0.092,0.54]
    # ---------------------------------------
    # Degrees of freedom: 19
    # Standardized standard error: 0.107
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result['minimum_detectable_effect'] == pytest.approx(0.316, abs=0.001)
    assert result['95% Confidence Interval'][0] == pytest.approx(0.092, abs=0.001)
    assert result['95% Confidence Interval'][1] == pytest.approx(0.540, abs=0.001)


def test_bira2c1() -> None:
    result = mde.mde_bira2c1(n=15, J=20, alpha=0.05)
    # mdes.bira2c1(n=15, J=20)
    #
    # Minimum detectable effect size:
    # ---------------------------------------
    #  0.325 95% CI [0.097,0.552]
    # ---------------------------------------
    # Degrees of freedom: 279
    # Standardized standard error: 0.115
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result['minimum_detectable_effect'] == pytest.approx(0.325, abs=0.001)
    assert result['95% Confidence Interval'][0] == pytest.approx(0.097, abs=0.001)
    assert result['95% Confidence Interval'][1] == pytest.approx(0.552, abs=0.001)


def test_bira2f1() -> None:
    result = mde.mde_bira2c1(n=15, J=20, alpha=0.05)
    # mdes.bira2f1(n=15, J=20)
    #
    # Minimum detectable effect size:
    # ---------------------------------------
    #  0.325 95% CI [0.097,0.552]
    # ---------------------------------------
    # Degrees of freedom: 260
    # Standardized standard error: 0.115
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result['minimum_detectable_effect'] == pytest.approx(0.325, abs=0.001)
    assert result['95% Confidence Interval'][0] == pytest.approx(0.097, abs=0.001)
    assert result['95% Confidence Interval'][1] == pytest.approx(0.552, abs=0.001)


def test_bira2r1() -> None:
    result = mde.mde_bira2r1(rho2=.17, omega2=.50, n=15, J=20, alpha=0.05)
    # mdes.bira2r1(rho2=.17, omega2=.50, n=15, J=20)
    #
    # Minimum detectable effect size:
    # ---------------------------------------
    #  0.366 95% CI [0.107,0.625]
    # ---------------------------------------
    # Degrees of freedom: 19
    # Standardized standard error: 0.124
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result['minimum_detectable_effect'] == pytest.approx(0.366, abs=0.001)
    assert result['95% Confidence Interval'][0] == pytest.approx(0.107, abs=0.001)
    assert result['95% Confidence Interval'][1] == pytest.approx(0.625, abs=0.001)


def test_bira3r1() -> None:
    result = mde.mde_bira3r1(rho3=.20, rho2=.15, omega3=.10, omega2=.10, n=69, J=10, K=100, alpha=0.05)
    # mdes.bira3r1(rho3=.20, rho2=.15, omega3=.10, omega2=.10, n=69, J=10, K=100)
    #
    # Minimum detectable effect size:
    # ---------------------------------------
    #  0.045 95% CI [0.013,0.077]
    # ---------------------------------------
    # Degrees of freedom: 99
    # Standardized standard error: 0.016
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result['minimum_detectable_effect'] == pytest.approx(0.045, abs=0.001)
    assert result['95% Confidence Interval'][0] == pytest.approx(0.013, abs=0.001)
    assert result['95% Confidence Interval'][1] == pytest.approx(0.077, abs=0.001)


def test_bira4r1() -> None:
    result = mde.mde_bira4r1(rho4=.05, rho3=.15, rho2=.15, omega4=.50, omega3=.50, omega2=.50, n=10, J=4, K=4, L=27, alpha=0.05)
    # mdes.bira4r1(rho4=.05, rho3=.15, rho2=.15, omega4=.50, omega3=.50, omega2=.50, n=10, J=4, K=4, L=27)
    #
    # Minimum detectable effect size:
    # ---------------------------------------
    #  0.142 95% CI [0.042,0.243]
    # ---------------------------------------
    # Degrees of freedom: 26
    # Standardized standard error: 0.049
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result['minimum_detectable_effect'] == pytest.approx(0.142, abs=0.001)
    assert result['95% Confidence Interval'][0] == pytest.approx(0.042, abs=0.001)
    assert result['95% Confidence Interval'][1] == pytest.approx(0.243, abs=0.001)


def test_cra2r2() -> None:
    result = mde.mde_cra2r2(rho2=.17, n=15, J=20, alpha=0.05)
    # mdes.cra2r2(rho2=.17, n=15, J=20)
    #
    # Minimum detectable effect size:
    # ---------------------------------------
    #  0.629 95% CI [0.183,1.075]
    # ---------------------------------------
    # Degrees of freedom: 18
    # Standardized standard error: 0.212
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result['minimum_detectable_effect'] == pytest.approx(0.629, abs=0.001)
    assert result['95% Confidence Interval'][0] == pytest.approx(0.183, abs=0.001)
    assert result['95% Confidence Interval'][1] == pytest.approx(1.075, abs=0.001)


def test_cra3r3() -> None:
    result = mde.mde_cra3r3(rho3=.06, rho2=.17, n=15, J=3, K=60, alpha=0.05)
    # mdes.cra3r3(rho3=.06, rho2=.17, n=15, J=3, K=60)
    #
    # Minimum detectable effect size:
    # ---------------------------------------
    #  0.269 95% CI [0.08,0.458]
    # ---------------------------------------
    # Degrees of freedom: 58
    # Standardized standard error: 0.094
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result['minimum_detectable_effect'] == pytest.approx(0.269, abs=0.001)
    assert result['95% Confidence Interval'][0] == pytest.approx(0.080, abs=0.001)
    assert result['95% Confidence Interval'][1] == pytest.approx(0.458, abs=0.001)


def test_ira1r1() -> None:
    result = mde.mde_ira1r1(n=250, alpha=0.05)
    # mdes.ira1r1(n=250)
    #
    # Minimum detectable effect size:
    # ---------------------------------------
    #  0.356 95% CI [0.107,0.605]
    # ---------------------------------------
    # Degrees of freedom: 248
    # Standardized standard error: 0.126
    # Type I error rate: 0.05
    # Type II error rate: 0.2
    # Two-tailed test: TRUE
    assert result['minimum_detectable_effect'] == pytest.approx(0.356, abs=0.001)
    assert result['95% Confidence Interval'][0] == pytest.approx(0.107, abs=0.001)
    assert result['95% Confidence Interval'][1] == pytest.approx(0.605, abs=0.001)