from dataclasses import dataclass
import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2
import pandas as pd

@dataclass(frozen=True)
class FitResult:
    p_th: float
    p_th_err: float
    nu: float
    a: float
    b: float
    c: float
    # D: float
    # mu: float
    chi2_red: float
    p_value: float
    x_window: float
    distances: list[float]

def fss(X, p_th, nu, a, b, c):
    p, d = X
    x = (p - p_th) * d ** (1.0 / nu)
    return a + (b * x) + (c * x * x)

def fss_with_correction(X, p_th, nu, a, b, c, D, mu):
    p, d = X
    x = (p - p_th) * d ** (1.0 / nu)
    return a + (b * x) + (c * x * x) + D * d ** (-1.0 / mu)

def estimate_threshold(points: pd.DataFrame, seed: int = 0, n_boot: int = 5000) -> FitResult:
    ps = points["p"].to_numpy()
    ds = points["d"].to_numpy()
    pls = points["pl"].to_numpy()
    sigs = points["sigma"].to_numpy()
    shots = points["shots"].to_numpy()

    distances = [11, 13, 15]
    keep = np.isin(ds, distances)
    ps, ds, pls, sigs, shots = ps[keep], ds[keep], pls[keep], sigs[keep], shots[keep]

    # p_hi = float(np.max(ps))
    # lower = [0.0,    0.5, -np.inf, -np.inf, -np.inf, -np.inf, 0.3]
    # upper = [p_hi,   5.0,  np.inf,  np.inf,  np.inf,  np.inf, 10.0]

    p0 = [float(np.median(ps)), 1.0, float(np.median(pls)), 0.0, 0.0]
    # p0_with_correction = [float(np.median(ps)), 1.0, float(np.median(pls)), 0.0, 0.0, 0.0, 1.5]

    popt, pcov = curve_fit(
        fss, (ps, ds), pls, sigma=sigs, p0=p0, 
        maxfev=50000, absolute_sigma=True
    )

    x_window = 0.005
    x_abs = np.abs((ps - popt[0]) * ds ** (1.0 / popt[1]))
    in_win = x_abs <= x_window
    ps, ds, pls, sigs, shots = ps[in_win], ds[in_win], pls[in_win], sigs[in_win], shots[in_win]

    popt, pcov = curve_fit(
        fss, (ps, ds), pls, sigma=sigs, p0=p0, 
        maxfev=50000, absolute_sigma=True
    )

    # Reduced chi-square
    resid = (pls - fss((ps, ds), *popt)) / sigs
    dof = max(len(pls) - len(popt), 1)
    chi2_red = float(np.sum(resid ** 2) / dof)
    p_value = chi2.sf(chi2_red * dof, dof)

    # 1-sigma error on p_th from the fit covariance
    p_th_err_cov = float(np.sqrt(pcov[0][0]))

    # Parametric bootstrap
    rng = np.random.default_rng(seed)
    p_clip = np.clip(pls, 0.0, 1.0)
    boot_p_th: list[float] = []
    for _ in range(n_boot):
        boot_pl = rng.binomial(shots, p_clip) / shots
        bopt, _ = curve_fit(
            fss, (ps, ds), boot_pl, sigma=sigs, p0=popt,
            maxfev=50000, absolute_sigma=True,
        )
        boot_p_th.append(float(bopt[0]))
    p_th_err_boot = float(np.std(boot_p_th))
    p_th_err = max(p_th_err_cov, p_th_err_cov * np.sqrt(chi2_red), p_th_err_boot)

    return FitResult(
        p_th=popt[0],
        p_th_err=p_th_err,
        nu=popt[1],
        a=popt[2], b=popt[3], c=popt[4], 
        chi2_red=chi2_red,
        p_value=p_value,
        x_window=x_window,
        distances=distances,
        # D=popt[5],
        # mu=popt[6]
    )
