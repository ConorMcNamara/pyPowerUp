import pytest

from pyPowerUp import sample_size


def test_bcra3f2() -> None:
    # One discrepancy: the R code calls round() whereas I call ceil(). What this means is that I will always round up
    # to the nearest whole number whereas R will round down if the decimal is < 0.5.
    result = sample_size.sample_size_bcra3f2(effect_size=.145, rho2=.10, n=20, J=44, alpha=0.05)
    # mrss.bcra3f2(es = .145, rho2=.10, n=20, J=44)
    #
    # K = 5
    assert result == pytest.approx(5, abs=1)


def test_bcra3r2() -> None:
    result = sample_size.sample_size_bcra3r2(effect_size=.246, rho3=.13, rho2=.10, omega3=.4, n=10, J=6, alpha=0.05)
    # mrss.bcra3r2(es = .246, rho3=.13, rho2=.10, omega3=.4, n=10, J=6, alpha = 0.05)
    #
    # K = 24
    assert result == pytest.approx(24, abs=1)


def test_bcra4f3() -> None:
    result = sample_size.sample_size_bcra4f3(effect_size=0.339, rho3=.15, rho2=.15, n=10, J=4, K=4, alpha=0.05)
    # mrss.bcra4f3(es=0.339, rho3=.15, rho2=.15, n=10, J=4, K=4)
    #
    # L = 15
    assert result == pytest.approx(15, abs=1)


def test_bcra4r2() -> None:
    result = sample_size.sample_size_bcra4r2(effect_size=.206, rho4=.05, rho3=.15, rho2=.15, omega4=.50, omega3=.50,
                                             n=10, J=4, K=4, alpha=0.05)
    # mrss.bcra4r2(es = .206, rho4=.05, rho3=.15, rho2=.15, omega4=.50, omega3=.50, n=10, J=4, K=4)
    #
    # L = 20
    assert result == pytest.approx(20, abs=1)


def test_bcra4r3() -> None:
    result = sample_size.sample_size_bcra4r3(effect_size=0.316, rho4=.05, rho3=.15, rho2=.15, omega4=.50, n=10, J=4,
                                             K=4,
                                             alpha=0.05)
    # mrss.bcra4r3(es = .316, rho4=.05, rho3=.15, rho2=.15, omega4=.50, n=10, J=4, K=4)
    #
    # L = 20
    assert result == pytest.approx(20, abs=1)


def test_bira2c1() -> None:
    result = sample_size.sample_size_bira2c1(effect_size=0.325, n=15, alpha=0.05)
    # mrss.bira2c1(es=.325, n=15)
    #
    # J = 20
    assert result == pytest.approx(20, abs=1)


def test_bira2f1() -> None:
    result = sample_size.sample_size_bira2f1(effect_size=0.325, n=15, alpha=0.05)
    # mrss.bira2f1(es=.325, n=15)
    #
    # J = 20
    assert result == pytest.approx(20, abs=1)


def test_bira2r1() -> None:
    result = sample_size.sample_size_bira2r1(effect_size=.366, rho2=.17, omega2=.50, n=15, alpha=0.05)
    # mrss.bira2r1(es=.366, rho2=.17, omega2=.50, n=15)
    #
    # J = 20
    assert result == pytest.approx(20, abs=1)


def test_bira3r1() -> None:
    result = sample_size.sample_size_bira3r1(effect_size=.045, rho3=.20, rho2=.15, omega3=.10, omega2=.10, n=69, J=10,
                                             alpha=0.05)
    # mrss.bira3r1(es = .045, rho3=.20, rho2=.15, omega3=.10, omega2=.10, n=69, J=10)
    #
    # K = 100
    assert result == pytest.approx(100, abs=1)


def test_bira4r1() -> None:
    result = sample_size.sample_size_bira4r1(effect_size=0.142, rho4=.05, rho3=.15, rho2=.15, omega4=.50, omega3=.50,
                                             omega2=.50, n=10, J=4, K=4, alpha=0.05)
    # mrss.bira4r1(es = 0.142, rho4=.05, rho3=.15, rho2=.15, omega4=.50, omega3=.50, omega2=.50, n=10, J=4, K=4)
    #
    # L = 27
    assert result == pytest.approx(27, abs=1)


def test_cra2r2() -> None:
    result = sample_size.sample_size_cra2r2(effect_size=.629, rho2=.17, n=15, alpha=0.05)
    # mrss.cra2r2(es=.629, rho2=.17, n=15)
    #
    # J = 20
    assert result == pytest.approx(20, abs=1)


def test_cra3r3() -> None:
    result = sample_size.sample_size_cra3r3(effect_size=.269, rho3=.06, rho2=.17, n=15, J=3, alpha=0.05)
    # mrss.cra3r3(es=.269, rho3=.06, rho2=.17, n=15, J=3)
    #
    # K = 60
    assert result == pytest.approx(60, abs=1)


def test_cra4r4() -> None:
    result = sample_size.sample_size_cra4r4(effect_size=.412, rho4=.05, rho3=.05, rho2=.10, n=10, J=2, K=3, alpha=0.05)
    # mrss.cra4r4(es = .412, rho4=.05, rho3=.05, rho2=.10, n=10, J=2, K=3)
    #
    # L = 20
    assert result == pytest.approx(20, abs=1)


def test_ira1r1() -> None:
    result = sample_size.sample_size_ira1r1(effect_size=0.356, alpha=0.05)
    # mrss.ira1r1(es=.356)
    #
    # n = 250
    assert result == pytest.approx(250, abs=1)


if __name__ == "__main__":
    pytest.main()