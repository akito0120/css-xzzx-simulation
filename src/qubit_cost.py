import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from threshold import crossing_seed, estimate_threshold

def physical_qubits(d: int) -> int:
    return 2 * d * d - 1

KAPPA = 0.5

def supp_model(X, lnA):
    d, ln_ratio = X  # ln_ratio = ln(p / p_th)
    return lnA + KAPPA * d * ln_ratio

def suppression_fit(sub: pd.DataFrame, p_th: float):
    d = sub["d"].to_numpy(dtype=float)
    p = sub["p"].to_numpy(dtype=float)
    pl = sub["pl"].to_numpy(dtype=float)
    sig = sub["sigma"].to_numpy(dtype=float)
    ln_eps = np.log(pl / d)
    ln_ratio = np.log(p / p_th)
    sig_ln = sig / pl
    sig_ln = np.where(sig_ln > 0, sig_ln, np.nanmin(sig_ln[sig_ln > 0]))
    popt, pcov = curve_fit(
        supp_model, (d, ln_ratio), ln_eps, sigma=sig_ln,
        p0=[float(np.median(ln_eps))], absolute_sigma=True, maxfev=50000,
    )
    return float(popt[0]), KAPPA, pcov

def required_distance(lnA: float, kappa: float, p_th: float, p_common: float, eps_target: float) -> int:
    # Invert ln eps = lnA + kappa*d*ln(p/p_th) for smallest odd d with eps <= eps_target.
    denom = kappa * np.log(p_common / p_th)
    if not (denom < 0) or not np.isfinite(denom):
        return -1  # not suppressing (p >= p_th, or bad fit)
    d = (np.log(eps_target) - lnA) / denom
    if not np.isfinite(d) or d <= 0:
        return -1
    d = int(np.ceil(d))
    if d % 2 == 0:
        d += 1
    return max(d, 3)

def subthreshold_slice(g: pd.DataFrame, p_th: float, margin: float = 1.0, min_points: int = 6):
    sub = g[(g["errors"] > 0) & (g["p"] < margin * p_th)]
    if len(sub) < min_points or sub["d"].nunique() < 3:
        return None
    return sub

def qubit_cost_row(sub: pd.DataFrame, p_th: float, p_common: float, eps_target: float, n_boot: int, rng) -> dict:
    lnA, kappa, _ = suppression_fit(sub, p_th)
    d_star = required_distance(lnA, kappa, p_th, p_common, eps_target)
    n = physical_qubits(d_star) if d_star > 0 else -1

    d = sub["d"].to_numpy(dtype=float)
    p = sub["p"].to_numpy(dtype=float)
    pl = sub["pl"].to_numpy(dtype=float)
    sig = sub["sigma"].to_numpy(dtype=float)
    shots = sub["shots"].to_numpy(dtype=float)
    ln_ratio = np.log(p / p_th)
    sig_ln = np.where(sig / pl > 0, sig / pl, 1.0)
    ns: list[int] = []
    for _ in range(n_boot):
        boot_err = rng.binomial(shots.astype(int), np.clip(pl, 0.0, 1.0))
        boot_pl = np.maximum(boot_err, 0.5) / shots  # continuity floor avoids log(0)
        y = np.log(boot_pl / d)
        try:
            popt, _ = curve_fit(supp_model, (d, ln_ratio), y, sigma=sig_ln,
                                p0=[lnA], absolute_sigma=True, maxfev=50000)
        except (RuntimeError, ValueError):
            continue
        d_b = required_distance(popt[0], KAPPA, p_th, p_common, eps_target)
        if d_b > 0:
            ns.append(physical_qubits(d_b))
    if ns:
        arr = np.array(ns, dtype=float)
        n_err = float(np.std(arr))
        n_lo, n_hi = float(np.percentile(arr, 16)), float(np.percentile(arr, 84))
    else:
        n_err = n_lo = n_hi = float("nan")

    return {
        "p_th": p_th,
        "kappa": kappa,
        "d_star": d_star,
        "n": n,
        "n_err": n_err,
        "n_lo": n_lo,
        "n_hi": n_hi,
        "n_points": int(len(sub)),
    }

def qubit_cost_table(df: pd.DataFrame, p_common: float, eps_target: float, n_boot: int = 1000, seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    for code in df["code"].unique():
        for eta in sorted(df[df["code"] == code]["eta"].unique()):
            g = df[(df["code"] == code) & (df["eta"] == eta)]
            p_th = estimate_threshold(g, n_boot=1).p_th
            sub = subthreshold_slice(g, p_th)
            row = qubit_cost_row(sub, p_th, p_common, eps_target, n_boot, rng)
            row.update({"code": code, "eta": float(eta)})
            rows.append(row)
    cols = ["eta", "code", "p_th", "kappa", "d_star", "n", "n_err", "n_lo", "n_hi", "n_points"]
    return pd.DataFrame(rows)[cols].sort_values(["code", "eta"]).reset_index(drop=True)
