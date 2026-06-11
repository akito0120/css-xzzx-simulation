import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict

from config import DISTANCES
from code_builder import build_rotated_surface_code, build_xzzx_code
from circuit_builder import CodeCapacityCircuitBuilder
from simulation import wilson_interval
from threshold import estimate_threshold, estimate_all_thresholds
from rich.status import Status
import pandas as pd
from rich.console import Console, Group
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

def fmt_eta(eta: float) -> str:
    if eta == float("inf"):
        return "inf"
    return str(int(eta)) if float(eta).is_integer() else str(eta)

def draw_collapse(path, eta, df: pd.DataFrame, colors):
    # Rescale each point to x = (p - p_th) * d^(1/nu) and plot p_L against it
    # Under the FSS ansatz all distances collapse onto a single curve
    fig, axes = plt.subplots(1, 2, figsize=(18, 8), dpi=600)
    for ax, (code, name) in zip(axes, [("css", "CSS"), ("xzzx", "XZZX")]):
        code_df = df[df["code"] == code]
        if code_df.empty:
            continue
        fit = estimate_threshold(code_df)
        x_all, pL_all = [], []
        for i, (d, d_df) in enumerate(code_df.groupby("d")):
            ps_arr = d_df["p"].to_numpy()
            ds_arr = d_df["d"].to_numpy()
            p_Ls = d_df["pl"].to_numpy()
            errs = d_df["errors"].to_numpy()
            measured = errs > 0
            x = (ps_arr - fit.p_th) * ds_arr ** (1.0 / fit.nu)
            ax.scatter(x[measured], p_Ls[measured], s=18, color=colors[code][i], label=f"{code}_d{d}")
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

def draw_pl(ax, label, points: pd.DataFrame, color):
    # Draw logical error rates
    # Error bars are the 1-sigma Wilson interval (z=1.0)
    ps = points["p"].to_numpy()
    pls = points["pl"].to_numpy()
    errs = points["errors"].to_numpy()
    shots = points["shots"].to_numpy()
    bounds = np.array([wilson_interval(e, s, z=1.0) for e, s in zip(errs, shots)]).reshape(-1, 2)
    lows, highs = bounds[:, 0], bounds[:, 1]

    measured = errs > 0
    zero = ~measured

    # Measured points: connect with a line and show Wilson error bars
    if np.any(measured):
        ax.errorbar(
            ps[measured], pls[measured],
            yerr=[pls[measured] - lows[measured], highs[measured] - pls[measured]],
            marker='o', linestyle='-', capsize=3,
            label=label, color=color,
        )

    # Zero-failure points: plot the Wilson upper bound as a downward arrow
    if np.any(zero):
        ax.errorbar(
            ps[zero], highs[zero],
            yerr=highs[zero] * 0.5, uplims=True,
            marker='', linestyle='none', color=color,
            label=None if np.any(measured) else label,
        )

def render_eta(eta: float, df: pd.DataFrame, outdir: str) -> None:
    eta_label = fmt_eta(eta)

    n_dist = df["d"].nunique()
    colors = {
        "css": plt.cm.plasma(np.linspace(0.2, 0.8, n_dist + 1)),
        "xzzx": plt.cm.viridis(np.linspace(0.2, 0.8, n_dist + 1))
    }
    fig, ax = plt.subplots(figsize=(12, 8), dpi=600)
    for code, code_df in df.groupby("code"):
        for i, (d, points) in enumerate(code_df.groupby("d")):
            draw_pl(ax, f"{code}_d{d}", points.sort_values("p", ascending=True), colors[code][i])

    # Draw thresholds (vertical line + 1-sigma uncertainty band)
    for code, points in df.groupby("code"):
        color = colors[code][points["d"].nunique()]
        fit = estimate_threshold(points)
        ax.axvline(x=fit.p_th, linestyle="--", label=f"{code} p_th = {fit.p_th:.6f} ± {fit.p_th_err:.6f}", color=color)
        ax.axvspan(fit.p_th - fit.p_th_err, fit.p_th + fit.p_th_err, alpha=0.1, color=color)

    ax.set_xscale('linear')
    ax.set_yscale('log')

    ax.set_xlabel('Physical Error Rate')
    ax.set_ylabel('Logical Error Rate (error bars: 1σ Wilson)')
    ax.set_title(f'Rotated CSS Surface Code vs XZZX Code Simulation Results with η = {eta_label}')

    ax.grid(True, which="both", alpha=0.5)
    ax.legend()

    plt.savefig(f"{outdir}/result_{eta_label}.png")
    plt.close(fig)

    draw_collapse(f"{outdir}/collapse_{eta_label}.png", eta_label, df, colors)

def render_threshold(points: pd.DataFrame, outdir: str):
    os.makedirs(outdir, exist_ok=True)
    etas = points["eta"].unique().tolist()

    # (eta, code) -> (eta, threshold, threshold error)
    thresholds: Dict[str, List[Tuple[float, float, float]]] = {"css": [], "xzzx": []}
    for (eta, code), fit_points in points.groupby(["eta", "code"]):
        if fit_points.empty:
            continue
        fit = estimate_threshold(fit_points)
        thresholds[code].append((float(eta), fit.p_th, fit.p_th_err))

    # Place inf one decade past the largest finite eta
    finite = [e for e in etas if e != float("inf")]
    inf_x = (max(finite) * 10.0) if finite else 1.0
    eta_to_x = lambda e: inf_x if e == float("inf") else e

    fig, ax = plt.subplots(figsize=(10, 7), dpi=600)
    styles = {"css": ("darkorange", "CSS"), "xzzx": ("limegreen", "XZZX")}

    for code_type, pts in thresholds.items():
        if not pts:
            continue
        color, label = styles[code_type]
        xs = [eta_to_x(e) for e, _, _ in pts]
        p_ths = [pth for _, pth, _ in pts]
        errs = [err for _, _, err in pts]
        ax.errorbar(xs, p_ths, yerr=errs, marker="o", linestyle="-", capsize=3, color=color, label=label)

    ax.set_xscale("log")
    ax.set_xticks([eta_to_x(e) for e in etas])
    ax.set_xticklabels([r"$\infty$" if e == float("inf") else fmt_eta(e) for e in etas])

    ax.set_xlabel(r"Bias $\eta = p_Z / (p_X + p_Y)$")
    ax.set_ylabel(r"Threshold $p_{th}$ (error bars: 1$\sigma$)")
    ax.set_title("Threshold vs bias: Rotated CSS Surface Code vs XZZX Code")
    ax.grid(True, which="both", alpha=0.4)
    ax.legend()

    fig.savefig(f"{outdir}/threshold.png")
    plt.close(fig)

def render_figures(df: pd.DataFrame, outdir):
    with Status("Rendering results", spinner="arc"):
        os.makedirs(outdir, exist_ok=True)
        mpl.rcParams["font.family"] = "serif"
        mpl.rcParams["font.size"] = 12
        mpl.rcParams["lines.linewidth"] = 1.5

        for eta, eta_df in df.groupby("eta"):
            render_eta(float(eta), eta_df, outdir)

        render_threshold(df, outdir)
        print(f"☑ Results saved to {outdir}")

def render_diagrams(outdir: str) -> None:
    os.makedirs(outdir, exist_ok=True)
    
    for distance in DISTANCES:
        for code in [build_rotated_surface_code(distance), build_xzzx_code(distance)]:
            circuit = CodeCapacityCircuitBuilder(code, 0.1, 0.5).build()

            detslice = circuit.diagram("detslice-svg")
            timeline = circuit.diagram("timeline-svg")

            with open(f"{outdir}/{code.name}_detslice.svg", "w") as f:
                f.write(str(detslice))
            with open(f"{outdir}/{code.name}_timeline.svg", "w") as f:
                f.write(str(timeline))

def print_summary(df: pd.DataFrame):
    console = Console()
    table = Table(expand=True)

    etas = df["eta"].unique().tolist()
    ds = df["d"].unique().tolist()
    codes = df["code"].unique().tolist()
    parameters = Text(f"Parameters\n - eta = {etas}\n - distances = {ds}\n - code types = {codes}")

    total_samples = Text(f"Total number of samples: {len(df)}")
    zero_error_count = Text(f"Number of zero-error samples: {str((df["errors"] == 0).sum())}")

    for col in ["eta", "code", "threshold", "error", "reduced chi-square"]:
        table.add_column(col)
    
    thresholds = estimate_all_thresholds(df)
    for (eta, code), fit in thresholds.items():
        table.add_row(str(eta), str(code), str(fit.p_th), str(fit.p_th_err), str(fit.chi2_red))

    content = Group(parameters, total_samples, zero_error_count, table)
    console.print(Panel(content, title="SUMMARY"))
