import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict

from config import DISTANCES, ETAS, PHYSICAL_ERROR_RATES
from code_builder import build_rotated_surface_code, build_xzzx_code
from circuit_builder import CodeCapacityCircuitBuilder
from simulation import wilson_interval
from threshold import estimate_threshold, SamplePoint

def fmt_eta(eta: float) -> str:
    if eta == float("inf"):
        return "inf"
    return str(int(eta)) if float(eta).is_integer() else str(eta)

def group_results(pairs: List[Tuple[str, SamplePoint]], code_type: str):
    # code_type + distance -> sample points (sorted by p), plus a flat list for the FSS fit.
    points = [sp for ct, sp in pairs if ct == code_type]
    results: Dict[str, List[SamplePoint]] = {}
    for distance in sorted({sp.distance for sp in points}):
        d_points = sorted(
            (sp for sp in points if sp.distance == distance),
            key=lambda sp: sp.physical_error_rate,
        )
        results[f"{code_type}_d{distance}"] = d_points
    return results, points

def draw_collapse(path, eta, css, xzzx):
    # Data-collapse plot: rescale each point to x = (p - p_th) * d^(1/nu) and plot p_L against it
    # Under the FSS ansatz all distances fall on a single curve
    fig, axes = plt.subplots(1, 2, figsize=(18, 8), dpi=600)
    for ax, (results, fit, colors, name) in zip(axes, (css, xzzx)):
        x_all, pL_all = [], []
        for i, (label, points) in enumerate(results.items()):
            ps_arr = np.array([sp.physical_error_rate for sp in points])
            ds_arr = np.array([sp.distance for sp in points])
            p_Ls = np.array([sp.logical_error_rate for sp in points])
            errs = np.array([sp.logical_errors for sp in points])
            measured = errs > 0
            x = (ps_arr - fit.p_th) * ds_arr ** (1.0 / fit.nu)
            ax.scatter(x[measured], p_Ls[measured], s=18, color=colors[i], label=label)
            x_all.append(x[measured]); pL_all.append(p_Ls[measured])
        
        # Overlay the fitted scaling function f(x) = a + b x + c x^2
        # drawn only over the span of the actual data
        a, b, c = fit.popt[2], fit.popt[3], fit.popt[4]
        x_all = np.concatenate(x_all) if x_all else np.array([0.0])
        pL_all = np.concatenate(pL_all) if pL_all else np.array([0.0])
        xs = np.linspace(x_all.min(), x_all.max(), 200)
        ax.plot(xs, a + b * xs + c * xs ** 2, color="black", lw=1.0,
                label="fitted f(x)")
        if pL_all.size:
            ax.set_ylim(0, 1.15 * float(pL_all.max()))
        ax.set_xlabel(r"$(p - p_{th})\, d^{1/\nu}$")
        ax.set_ylabel("Logical Error Rate")
        ax.set_title(
            f"{name}: p_th={fit.p_th:.4f}±{fit.p_th_err:.4f}, "
            f"ν={fit.nu:.2f}±{fit.nu_err:.2f}, χ²_red={fit.chi2_red:.2f}, "
            f"n={fit.n_points}"
        )
        ax.grid(True, alpha=0.4)
        ax.legend(fontsize=8)

    fig.suptitle(f"FSS data collapse (η = {eta})")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)

def render_eta(eta: float, pairs: List[Tuple[str, SamplePoint]], outdir: str) -> None:
    eta_label = fmt_eta(eta)
    css_results, css_sample_points = group_results(pairs, "css")
    xzzx_results, xzzx_sample_points = group_results(pairs, "xzzx")

    # FSS fitting
    css_fit = estimate_threshold(css_sample_points)
    xzzx_fit = estimate_threshold(xzzx_sample_points)

    n_dist = max(len(css_results), len(xzzx_results))
    css_colors = plt.cm.plasma(np.linspace(0.2, 0.8, n_dist + 1))
    xzzx_colors = plt.cm.viridis(np.linspace(0.2, 0.8, n_dist + 1))
    fig, ax = plt.subplots(figsize=(12, 8), dpi=600)

    # Draw logical error rates
    def draw_curve(label, points: List[SamplePoint], color):
        # Error bars are the 1-sigma Wilson interval (z=1.0)
        ps_arr = np.array([sp.physical_error_rate for sp in points])
        errs = np.array([sp.logical_errors for sp in points])
        lows = np.array([wilson_interval(sp.logical_errors, sp.shots, z=1.0)[0] for sp in points])
        highs = np.array([wilson_interval(sp.logical_errors, sp.shots, z=1.0)[1] for sp in points])
        p_Ls = np.array([sp.logical_error_rate for sp in points])

        measured = errs > 0
        zero = ~measured

        # Measured points: connect with a line and show Wilson error bars
        if np.any(measured):
            ax.errorbar(
                ps_arr[measured], p_Ls[measured],
                yerr=[p_Ls[measured] - lows[measured], highs[measured] - p_Ls[measured]],
                marker='o', linestyle='-', capsize=3,
                label=label, color=color,
            )

        # Zero-failure points: plot the Wilson upper bound as a downward arrow
        if np.any(zero):
            ax.errorbar(
                ps_arr[zero], highs[zero],
                yerr=highs[zero] * 0.5, uplims=True,
                marker='', linestyle='none', color=color,
                label=None if np.any(measured) else label,
            )

    for i, (label, points) in enumerate(css_results.items()):
        draw_curve(label, points, css_colors[i])
    for i, (label, points) in enumerate(xzzx_results.items()):
        draw_curve(label, points, xzzx_colors[i])

    # Draw thresholds (vertical line + 1-sigma uncertainty band)
    for fit, color, name in (
            (css_fit, css_colors[n_dist], "css"),
            (xzzx_fit, xzzx_colors[n_dist], "xzzx")
        ):
        ax.axvline(x=fit.p_th, color=color, linestyle="--",
                   label=f"{name} p_th = {fit.p_th:.4f} ± {fit.p_th_err:.4f}")
        ax.axvspan(fit.p_th - fit.p_th_err, fit.p_th + fit.p_th_err,
                   color=color, alpha=0.1)

    ax.set_xscale('linear')
    ax.set_yscale('log')

    ax.set_xlabel('Physical Error Rate')
    ax.set_ylabel('Logical Error Rate (error bars: 1σ Wilson)')
    ax.set_title(f'Rotated CSS Surface Code vs XZZX Code Simulation Results with η = {eta_label}')

    ax.grid(True, which="both", alpha=0.5)
    ax.legend()

    plt.savefig(f"{outdir}/result_{eta_label}.png")
    plt.close(fig)

    # Draw data-collapse figure: rescaled axis x = (p - p_th) d^(1/nu)
    draw_collapse(
        f"{outdir}/collapse_{eta_label}.png", eta_label,
        (css_results, css_fit, css_colors, "CSS"),
        (xzzx_results, xzzx_fit, xzzx_colors, "XZZX")
    )

def render_all(pairs: List[Tuple[str, SamplePoint]], outdir) -> None:
    os.makedirs(outdir, exist_ok=True)
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["font.size"] = 12
    mpl.rcParams["lines.linewidth"] = 1.5

    for eta in sorted({sp.eta for _, sp in pairs}):
        eta_pairs = [(ct, sp) for ct, sp in pairs if sp.eta == eta]
        render_eta(eta, eta_pairs, outdir)

def render_diagrams(outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    
    for distance in DISTANCES:
        eta = ETAS[0]
        p = PHYSICAL_ERROR_RATES[0]

        for code in [build_rotated_surface_code(distance), build_xzzx_code(distance)]:
            circuit = CodeCapacityCircuitBuilder(code, p, eta).build()

            detslice = circuit.diagram("detslice-svg")
            timeline = circuit.diagram("timeline-svg")

            with open(f"{outdir}/{code.name}_detslice.svg", "w") as f:
                f.write(str(detslice))
            with open(f"{outdir}/{code.name}_timeline.svg", "w") as f:
                f.write(str(timeline))
