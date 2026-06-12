from dataclasses import dataclass
from typing import List, Tuple
import warnings
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd

@dataclass(frozen=True)
class FitResult:
    p_th: float                     # threshold estimate
    p_th_err: float                 # 1-sigma uncertainty on p_th (the more conservative of pcov / bootstrap)
    p_th_err_cov: float             # uncertainty from the fit covariance (sqrt(pcov[0,0]))
    p_th_err_boot: float            # uncertainty from Monte-Carlo bootstrap (std of resampled p_th)
    p_th_ci: Tuple[float, float]    # 16/84 percentile bootstrap interval
    nu: float                       # critical exponent
    nu_err: float                   # critical exponent uncertainity
    chi2_red: float                 # reduced chi-square of the windowed fit (~1 is good)
    n_points: int                   # number of data points inside the fit window
    window: float                   # half-width of the fit window in p
    popt: Tuple[float, ...]
    pcov: Tuple[Tuple[float, ...], ...]

def fss(X, p_th, nu, a, b, c):
    # Approximation of the scaling function by a second-order polynomial
    p, d = X
    x = (p - p_th) * d ** (1.0 / nu)
    return a + (b * x) + (c * x * x)

def fit(ps, ds, p_Ls, sigs, p0, absolute_sigma):
    return curve_fit(
        fss, (ps, ds), p_Ls,
        sigma=sigs, p0=p0,
        absolute_sigma=absolute_sigma, maxfev=30000,
    )

def estimate_threshold(
    points: pd.DataFrame,
    window_frac: float = 0.25,
    n_boot: int = 500,
    seed: int = 0,
) -> FitResult:
    ps_all = points["p"].to_numpy()
    ds_all = points["d"].to_numpy()
    p_Ls_all = points["pl"].to_numpy()
    sigs_all = points["sigma"].to_numpy()
    errs_all = points["errors"].to_numpy()
    shots_all = points["shots"].to_numpy()

    # (1) keep only points that actually observed logical errors
    measured = errs_all > 0
    ps, ds, p_Ls = ps_all[measured], ds_all[measured], p_Ls_all[measured]
    sigs, shots = sigs_all[measured], shots_all[measured]
    # guard against zero/NaN weights breaking the least squares
    sigs = np.where(sigs > 0, sigs, np.nanmin(sigs[sigs > 0]) if np.any(sigs > 0) else 1.0)

    # (2) coarse fit to locate the threshold
    p0 = [float(np.median(ps)), 1.5, float(np.median(p_Ls)), 0.0, 0.0]
    popt0, _ = fit(ps, ds, p_Ls, sigs, p0, absolute_sigma=False)
    p_th0 = popt0[0]

    # (3) window around p_th^(0) and refit
    window = window_frac * float(p_th0)
    in_win = np.abs(ps - p_th0) <= window
    n_points = int(np.count_nonzero(in_win))
    if n_points < 6:
        warnings.warn(
            f"FSS fit window holds only {n_points} points for 5 free parameters "
            f"(p_th0={p_th0:.4f}, window={window:.4f}); the threshold error bar will "
            f"be unreliable. Sample more densely near the threshold or add larger distances.",
            stacklevel=2,
        )
        in_win = np.ones_like(ps, dtype=bool)  # fall back to all measured points
        n_points = int(np.count_nonzero(in_win))

    pw, dw, plw, sw = ps[in_win], ds[in_win], p_Ls[in_win], sigs[in_win]
    p0w = [p_th0, popt0[1], float(np.median(plw)), 0.0, 0.0]
    popt, pcov = fit(pw, dw, plw, sw, p0w, absolute_sigma=True)

    # (4) reduced chi-square (guard the degrees of freedom)
    resid = (plw - fss((pw, dw), *popt)) / sw
    dof = max(n_points - len(popt), 1)
    chi2_red = float(np.sum(resid ** 2) / dof)

    p_th_err_cov = float(np.sqrt(pcov[0, 0])) if np.isfinite(pcov[0, 0]) else float("nan")
    nu_err = float(np.sqrt(pcov[1, 1])) if np.isfinite(pcov[1, 1]) else float("nan")

    # (5) parametric bootstrap on the windowed points
    rng = np.random.default_rng(seed)
    shots_w = shots[in_win]
    boot_p_th: List[float] = []
    for _ in range(n_boot):
        boot_errs = rng.binomial(shots_w, np.clip(plw, 0.0, 1.0))
        boot_pL = boot_errs / shots_w
        try:
            bopt, _ = fit(pw, dw, boot_pL, sw, popt, absolute_sigma=True)
            boot_p_th.append(float(bopt[0]))
        except (RuntimeError, ValueError):
            continue
    if boot_p_th:
        boot_arr = np.array(boot_p_th)
        p_th_err_boot = float(np.std(boot_arr))
        p_th_ci = (float(np.percentile(boot_arr, 16)), float(np.percentile(boot_arr, 84)))
    else:
        p_th_err_boot = float("nan")
        p_th_ci = (float("nan"), float("nan"))

    finite_errs = [e for e in (p_th_err_cov, p_th_err_boot) if np.isfinite(e)]
    p_th_err = max(finite_errs) if finite_errs else float("nan")

    return FitResult(
        p_th=float(popt[0]),
        p_th_err=p_th_err,
        p_th_err_cov=p_th_err_cov,
        p_th_err_boot=p_th_err_boot,
        p_th_ci=p_th_ci,
        nu=float(popt[1]),
        nu_err=nu_err,
        chi2_red=chi2_red,
        n_points=n_points,
        window=window,
        popt=tuple(float(v) for v in popt),
        pcov=tuple(tuple(float(v) for v in row) for row in pcov),
    )

def estimate_all_thresholds(df: pd.DataFrame) -> dict[tuple[float, str], FitResult]:
    thresholds: dict[tuple[float, str], FitResult] = dict()
    for (eta, code), points in df.groupby(["eta", "code"]):
        fit = estimate_threshold(points)
        thresholds[(eta, code)] = fit
    return thresholds
