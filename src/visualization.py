import os
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

from config import DISTANCES
from code_builder import build_rotated_surface_code, build_xzzx_code
from circuit_builder import CodeCapacityCircuitBuilder
from simulation import wilson_interval
from threshold import estimate_threshold
import pandas as pd

def generate_pl_colors(n_dist: int):
    return {
        "css": plt.cm.plasma(np.linspace(0.2, 0.8, n_dist + 1)),
        "xzzx": plt.cm.viridis(np.linspace(0.2, 0.8, n_dist + 1))
    }

def generate_p_th_colors():
    return {
        "css": plt.cm.plasma(0.5),
        "xzzx": plt.cm.viridis(0.5)
    }

def fmt_eta(eta: float) -> str:
    if eta == float("inf"):
        return "inf"
    return str(int(eta)) if float(eta).is_integer() else str(eta)

def draw_collapse(path, eta, df: pd.DataFrame, colors):
    # Rescale each point to x = (p - p_th) * d^(1/nu) and plot p_L against it
    # Under the FSS ansatz all distances collapse onto a single curve
    fig, axes = plt.subplots(1, 2, figsize=(18, 10), dpi=600)
    for ax, (code, name) in zip(axes, [("css", "CSS"), ("xzzx", "XZZX")]):
        code_df = df[df["code"] == code]
        if code_df.empty:
            continue
        fit = estimate_threshold(code_df)
        x_all = []
        for i, (d, d_df) in enumerate(code_df.groupby("d")):
            if d not in fit.distances:
                continue

            ps = d_df["p"].to_numpy()
            ds = d_df["d"].to_numpy()
            # pls = d_df["pl"].to_numpy() - fit.D * ds ** (-1.0 / fit.mu)
            pls = d_df["pl"].to_numpy()

            x = (ps - fit.p_th) * ds ** (1.0 / fit.nu)
            x_abs = np.abs(x)
            keep = x_abs < fit.x_window

            ax.scatter(x[keep], pls[keep], s=18, color=colors[code][i], label=f"{code}_d{d}")
            x_all.append(x[keep]);

        # Overlay the fitted f(x)
        x_all = np.concatenate(x_all) if x_all else np.array([0.0])
        xs = np.linspace(x_all.min(), x_all.max(), 200)
        ax.plot(xs, fit.a + fit.b * xs + fit.c * xs ** 2, color="black", lw=1.0, label="fitted f(x)")
        ax.set_xlabel(r"$(p - p_{th})\, d^{1/\nu}$")
        ax.set_ylabel("Logical Error Rate")
        ax.set_title(
            f"{name}: " r"$p_{th}$" f"={fit.p_th:.6f}±{fit.p_th_err:.6f}, "
            r"$\nu$" f"={fit.nu:.2f}, "
            r"$\chi^2_{red}$" f"={fit.chi2_red:.6f}, "
            f"p-value={fit.p_value:.4f}"
        )
        ax.grid(True, alpha=0.4)
        ax.legend(fontsize=8)

    fig.suptitle(f"FSS data collapse (η = {eta})")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    print(f"- Rendered {path} ✔")

def draw_pl(ax, label, points: pd.DataFrame, color, linestyle="-", marker="o"):
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
            marker=marker, linestyle=linestyle, capsize=3,
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
    colors = generate_pl_colors(n_dist)
    p_th_colors = generate_p_th_colors()

    fig, axes = plt.subplots(1, 2, figsize=(22, 8), dpi=600, sharey=True)
    for ax, (code, name) in zip(axes, [("css", "CSS"), ("xzzx", "XZZX")]):
        code_df = df[df["code"] == code]
        if code_df.empty:
            continue
        for i, (d, points) in enumerate(code_df.groupby("d")):
            draw_pl(ax, f"{code}_d{d}", points.sort_values("p", ascending=True), colors[code][i])

        # Draw threshold (vertical line + 1-sigma uncertainty band)
        color = p_th_colors[code]
        fit = estimate_threshold(code_df)
        ax.axvline(x=fit.p_th, linestyle="--", label=f"{code} p_th = {fit.p_th:.6f} ± {fit.p_th_err:.6f}", color=color)
        ax.axvspan(fit.p_th - fit.p_th_err, fit.p_th + fit.p_th_err, alpha=0.1, color=color)

        ax.set_xscale('linear')
        ax.set_yscale('log')
        ax.set_xlabel('Physical Error Rate')
        ax.set_ylabel('Logical Error Rate (error bars: 1σ Wilson)')
        ax.set_title(name)
        ax.grid(True, which="both", alpha=0.5)
        ax.legend()

    fig.suptitle(f'Rotated CSS Surface Code vs XZZX Code Simulation Results with η = {eta_label}')
    fig.tight_layout()

    path = f"{outdir}/result_{eta_label}.pdf"
    fig.savefig(path)
    plt.close(fig)
    print(f"- Rendered {path} ✔")

    draw_collapse(f"{outdir}/collapse_{eta_label}.pdf", eta_label, df, colors)

def render_eta_basis(eta: float, df: pd.DataFrame, outdir: str) -> None:
    # Intra-code comparison of X- vs Z-memory logical error rate
    # Evidence that X-memory is the worst-case basis (its pl sits above Z-memory's)
    eta_label = fmt_eta(eta)

    n_dist = df["d"].nunique()
    colors = generate_pl_colors(n_dist)
    fig, axes = plt.subplots(1, 2, figsize=(22, 8), dpi=600)
    basis_style = {"x": ("-", "o"), "z": ("--", "s")}

    for ax, (code, name) in zip(axes, [("css", "CSS"), ("xzzx", "XZZX")]):
        code_df = df[df["code"] == code]
        for i, (d, d_df) in enumerate(code_df.groupby("d")):
            for basis, (ls, marker) in basis_style.items():
                points = d_df[d_df["basis"] == basis].sort_values("p", ascending=True)
                if points.empty:
                    continue
                draw_pl(ax, f"d{d} {basis.upper()}", points, colors[code][i],
                        linestyle=ls, marker=marker)
        ax.set_yscale("log")
        ax.set_xlabel("Physical Error Rate")
        ax.set_ylabel("Logical Error Rate (error bars: 1σ Wilson)")
        ax.set_title(name)
        ax.grid(True, which="both", alpha=0.5)
        ax.legend(fontsize=8)

    fig.suptitle(f"Memory-basis comparison (solid: X, dashed: Z) — η = {eta_label}")
    fig.tight_layout()
    fig.savefig(f"{outdir}/basis_{eta_label}.pdf")
    plt.close(fig)

def render_threshold(points: pd.DataFrame, outdir: str):
    os.makedirs(outdir, exist_ok=True)
    etas = points["eta"].unique().tolist()

    # (eta, code) -> (eta, threshold, threshold error)
    thresholds: dict[str, list[tuple[float, float, float]]] = {"css": [], "xzzx": []}
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
    colors = generate_p_th_colors()

    for code, pts in thresholds.items():
        if not pts:
            continue
        xs = [eta_to_x(e) for e, _, _ in pts]
        p_ths = [pth for _, pth, _ in pts]
        errs = [err for _, _, err in pts]
        ax.errorbar(xs, p_ths, yerr=errs, marker="o", linestyle="-", capsize=3, color=colors[code], label=code)

    ax.set_xscale("log")
    ax.set_xticks([eta_to_x(e) for e in etas])
    ax.set_xticklabels([r"$\infty$" if e == float("inf") else fmt_eta(e) for e in etas])

    ax.set_xlabel(r"Bias $\eta = p_Z / (p_X + p_Y)$")
    ax.set_ylabel(r"Threshold $p_{th}$ (error bars: 1$\sigma$)")
    ax.set_title("Threshold vs bias: Rotated CSS Surface Code vs XZZX Code")
    ax.grid(True, which="both", alpha=0.4)
    ax.legend()

    path = f"{outdir}/threshold.pdf"
    fig.savefig(path)
    plt.close(fig)
    print(f"- Rendered {path} ✔")

def render_figures(df: pd.DataFrame, outdir):
    print("Rendering results...")
    os.makedirs(outdir, exist_ok=True)
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["font.size"] = 12
    mpl.rcParams["lines.linewidth"] = 1.5

    # Draw inter-code comparison of pl and the thresholds for X-memory results
    for eta, eta_df in df.groupby("eta"):
        render_eta(float(eta), eta_df, outdir)
    render_threshold(df, outdir)

    print(f"Results saved to {outdir} ✔\n")

def render_diagrams(outdir: str) -> None:
    print("Rendering diagrams...")
    os.makedirs(outdir, exist_ok=True)
    
    for distance in DISTANCES:
        for code in [build_rotated_surface_code(distance), build_xzzx_code(distance)]:
            circuit = CodeCapacityCircuitBuilder(code, 0.1, 0.5).build()

            detslice = circuit.diagram("detslice-svg")
            timeline = circuit.diagram("timeline-svg")

            detslice_path = f"{outdir}/{code.name}_detslice.svg"
            timeline_path = f"{outdir}/{code.name}_timeline.svg"

            with open(detslice_path, "w") as f:
                f.write(str(detslice))
                print(f"- Rendered {detslice_path} ✔")
            with open(timeline_path, "w") as f:
                f.write(str(timeline))
                print(f"- Rendered {timeline_path} ✔")
    
    print(f"Diagrams saved to {outdir} ✔\n")
