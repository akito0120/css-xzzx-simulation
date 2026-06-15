from dataclasses import dataclass
from typing import List, Optional, Tuple
import warnings
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd

N_PARAMS = 7 # Number of free parameters in FSS (p_th, nu, a, b, c, D, mu)
WINDOW_MARGIN = 7 # A fit window must hold at least N_PARAMS + WINDOW_MARGIN points to be usable
DEFAULT_WINDOW_QUANTILES = (0.1, 0.12, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0) # Default |x| window cutoffs

@dataclass(frozen=True)
class FitResult:
    p_th: float                     # threshold estimate
    p_th_err: float                 # 1-sigma uncertainty on p_th (the more conservative of pcov / bootstrap)
    p_th_err_cov: float             # uncertainty from the fit covariance (sqrt(pcov[0,0]))
    p_th_err_boot: float            # uncertainty from Monte-Carlo bootstrap (std of resampled p_th)
    p_th_ci: Tuple[float, float]    # 16/84 percentile bootstrap interval
    nu: float                       # critical exponent
    nu_err: float                   # critical exponent uncertainity
    D: float                        # corrections-to-scaling amplitude
    D_err: float                    # 1-sigma uncertainty on D
    mu: float                       # corrections-to-scaling exponent (correction ~ d^(-1/mu))
    mu_err: float                   # 1-sigma uncertainty on mu
    chi2_red: float                 # reduced chi-square of the windowed fit (~1 is good)
    n_points: int                   # number of data points inside the fit window
    x_max: float                    # selected |x| half-width of the fit window (inf if fell back to all points)
    window: float                   # corresponding half-width in p ((p_max - p_min) / 2 over the window), for reference
    d_min: Optional[int]            # smallest code distance kept in the fit (None = all distances)
    popt: Tuple[float, ...]
    pcov: Tuple[Tuple[float, ...], ...]

def fss(X, p_th, nu, a, b, c, D, mu):
    # Quadratic approximation of the scaling function with a correction term
    p, d = X
    x = (p - p_th) * d ** (1.0 / nu)
    return a + (b * x) + (c * x * x) + D * d ** (-1.0 / mu)

def get_bounds(ps):
    # Keep p_th inside the swept p-range and the two exponents positive/finite so
    # the d^(1/nu) and d^(-1/mu) factors stay well behaved; a, b, c, D are free.
    p_hi = float(np.max(ps))
    lower = [0.0, 0.5, -np.inf, -np.inf, -np.inf, -np.inf, 0.3]
    upper = [p_hi, 5.0, np.inf, np.inf, np.inf, np.inf, 10.0]
    return (lower, upper)

def fit(ps, ds, p_Ls, sigs, p0, absolute_sigma, bounds):
    return curve_fit(
        fss, (ps, ds), p_Ls,
        sigma=sigs, p0=p0, bounds=bounds,
        absolute_sigma=absolute_sigma, maxfev=30000,
    )

def calc_chi2_red(pw, dw, plw, sw, popt) -> float:
    resid = (plw - fss((pw, dw), *popt)) / sw
    dof = max(len(plw) - N_PARAMS, 1)
    return float(np.sum(resid ** 2) / dof)

def select_window(ps, ds, p_Ls, sigs, p_th_c, nu_c, p0_seed, quantiles, bounds, min_points):
    x_abs = np.abs((ps - p_th_c) * ds ** (1.0 / nu_c))
    x_maxes = np.unique(np.quantile(x_abs, quantiles))
    best = None
    for x_max in x_maxes:
        in_win = x_abs <= x_max
        n_points = int(np.count_nonzero(in_win))
        if n_points < min_points:
            continue
        pw, dw, plw, sw = ps[in_win], ds[in_win], p_Ls[in_win], sigs[in_win]

        try:
            popt, pcov = fit(pw, dw, plw, sw, p0_seed, absolute_sigma=True, bounds=bounds)
        except (RuntimeError, ValueError):
            continue
        if not np.all(np.isfinite(pcov)):
            continue

        chi2_red = calc_chi2_red(pw, dw, plw, sw, popt)
        key = abs(chi2_red - 1.0)
        if best is None or key < best[0]:
            best = (key, float(x_max), in_win, popt, pcov, chi2_red, n_points)
    return best

def estimate_threshold(
    points: pd.DataFrame,
    d_min: Optional[int] = None,
    window_quantiles: Tuple[float, ...] = DEFAULT_WINDOW_QUANTILES,
    n_boot: int = 500,
    seed: int = 0,
) -> FitResult:
    ps_all = points["p"].to_numpy()
    ds_all = points["d"].to_numpy()
    p_Ls_all = points["pl"].to_numpy()
    sigs_all = points["sigma"].to_numpy()
    errs_all = points["errors"].to_numpy()
    shots_all = points["shots"].to_numpy()

    # (1) Keep only points that actually observed logical errors
    #     Optionally drop the smallest distances
    measured = errs_all > 0
    if d_min is not None:
        measured = measured & (ds_all >= d_min)
    ps, ds, p_Ls = ps_all[measured], ds_all[measured], p_Ls_all[measured]
    sigs, shots = sigs_all[measured], shots_all[measured]
    # guard against zero/NaN weights breaking the least squares
    sigs = np.where(sigs > 0, sigs, np.nanmin(sigs[sigs > 0]) if np.any(sigs > 0) else 1.0)

    bounds = get_bounds(ps)
    min_points = N_PARAMS + WINDOW_MARGIN

    # (2) Coarse fit to locate the threshold and a first nu
    nu0 = 1.5
    p0 = [float(np.median(ps)), nu0, float(np.median(p_Ls)), 0.0, 0.0, 0.0, nu0]
    popt0, _ = fit(ps, ds, p_Ls, sigs, p0, absolute_sigma=False, bounds=bounds)
    p_th0, nu0 = float(popt0[0]), float(popt0[1])

    # (3) Pick the |x| <= x_max window whose reduced chi-square is closest to 1
    #     Recompute x with the refitted nu once to relax the circularity
    best = select_window(ps, ds, p_Ls, sigs, p_th0, nu0, list(popt0), window_quantiles, bounds, min_points)
    if best is not None:
        popt_b = best[3]
        best2 = select_window(
            ps, ds, p_Ls, sigs, float(popt_b[0]), float(popt_b[1]),
            list(popt_b), window_quantiles, bounds, min_points,
        )
        if best2 is not None and best2[0] < best[0]:
            best = best2

    if best is not None:
        _, x_max, in_win, popt, pcov, chi2_red, n_points = best
    else:
        warnings.warn(
            f"No |x| window held at least {min_points} points for {N_PARAMS} free "
            f"parameters (p_th0={p_th0:.4f}, nu0={nu0:.2f}); the threshold error bar "
            f"will be unreliable. Sample more densely near the threshold, add larger "
            f"distances, or lower d_min.",
            stacklevel=2,
        )
        in_win = np.ones_like(ps, dtype=bool)
        n_points = int(np.count_nonzero(in_win))
        popt, pcov = fit(ps, ds, p_Ls, sigs, list(popt0), absolute_sigma=True, bounds=bounds)
        chi2_red = calc_chi2_red(ps, ds, p_Ls, sigs, popt)
        x_max = float("inf")

    pw, dw, plw, sw = ps[in_win], ds[in_win], p_Ls[in_win], sigs[in_win]
    window = float((pw.max() - pw.min()) / 2.0) if pw.size else float("nan")

    # (4) Uncertainties from the fit covariance
    def cov_err(i):
        return float(np.sqrt(pcov[i, i])) if np.isfinite(pcov[i, i]) else float("nan")
    p_th_err_cov = cov_err(0)
    nu_err = cov_err(1)
    D_err = cov_err(5)
    mu_err = cov_err(6)

    # (5) Parametric bootstrap on the windowed points
    rng = np.random.default_rng(seed)
    shots_w = shots[in_win]
    boot_p_th: List[float] = []
    for _ in range(n_boot):
        boot_errs = rng.binomial(shots_w, np.clip(plw, 0.0, 1.0))
        boot_pL = boot_errs / shots_w
        try:
            bopt, _ = fit(pw, dw, boot_pL, sw, popt, absolute_sigma=True, bounds=bounds)
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

    finite_errs = [err for err in (p_th_err_cov, p_th_err_boot) if np.isfinite(err)]
    p_th_err = max(finite_errs) if finite_errs else float("nan")

    return FitResult(
        p_th=float(popt[0]),
        p_th_err=p_th_err,
        p_th_err_cov=p_th_err_cov,
        p_th_err_boot=p_th_err_boot,
        p_th_ci=p_th_ci,
        nu=float(popt[1]),
        nu_err=nu_err,
        D=float(popt[5]),
        D_err=D_err,
        mu=float(popt[6]),
        mu_err=mu_err,
        chi2_red=chi2_red,
        n_points=n_points,
        x_max=float(x_max),
        window=window,
        d_min=d_min,
        popt=tuple(float(v) for v in popt),
        pcov=tuple(tuple(float(v) for v in row) for row in pcov),
    )

def estimate_all_thresholds(
    df: pd.DataFrame,
    d_min: Optional[int] = None,
    window_quantiles: Tuple[float, ...] = DEFAULT_WINDOW_QUANTILES,
) -> dict[tuple[float, str], FitResult]:
    thresholds: dict[tuple[float, str], FitResult] = dict()
    for (eta, code), points in df.groupby(["eta", "code"]):
        fit = estimate_threshold(points, d_min=d_min, window_quantiles=window_quantiles)
        thresholds[(eta, code)] = fit
    return thresholds
