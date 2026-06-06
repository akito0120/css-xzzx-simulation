import os
import numpy as np
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
from typing import Dict, List
from rich.progress import Progress, BarColumn, TextColumn

from code_builder import QecCode, build_rotated_surface_code, build_xzzx_code
from circuit_builder import CodeCapacityCircuitBuilder, CircuitLevelCircuitBuilder
from simulation import estimate_logical_error_rate, wilson_interval
from threshold import SamplePoint, estimate_threshold

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

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--max-shots", type=int, default=2_000_000)
    ap.add_argument("--target-errors", type=int, default=200)
    ap.add_argument("--batch-size", type=int, default=100_000)
    ap.add_argument("--diagram-only", action="store_true", default=False)
    args = ap.parse_args()

    etas = [30, 100, float("inf")]
    code_types = ["css", "xzzx"]
    distances = [3, 5, 7]
    physical_error_rates = ps = list(np.linspace(0.002, 0.04, 20))

    os.makedirs(f"{args.outdir}", exist_ok=True)

    if args.diagram_only == False:
        with Progress(
            TextColumn("{task.description}"),
            BarColumn(bar_width=30),
            TextColumn("{task.completed}/{task.total}")) as progress:

            eta_task = progress.add_task("Simulation", total=len(etas))

            for eta in etas:
                code_type_task = progress.add_task(f"Simulating with eta = {eta}", total=len(code_types))

                # Sample points for FSS fitting
                css_sample_points: List[SamplePoint] = []
                xzzx_sample_points: List[SamplePoint] = []

                # Simulated logical errors: code type + distance -> sample points
                css_results: Dict[str, List[SamplePoint]] = {}
                xzzx_results: Dict[str, List[SamplePoint]] = {}

                # Simulation
                for code_type in code_types:
                    distance_task = progress.add_task(f"Simulating {code_type}", total=len(distances))

                    for distance in distances:
                        p_error_task = progress.add_task(f"Simulating with d = {distance}", total=len(physical_error_rates))
                        result: List[SamplePoint] = []

                        for physicall_error_rate in physical_error_rates:

                            code = build_rotated_surface_code(distance) if code_type == "css" else build_xzzx_code(distance)
                            circuit = CircuitLevelCircuitBuilder(code, physicall_error_rate, eta).build()
                            p_L, sigma, errors, shots = estimate_logical_error_rate(
                                circuit,
                                max_shots=args.max_shots,
                                target_errors=args.target_errors,
                                batch_size=args.batch_size,
                            )

                            sample_point = SamplePoint(
                                eta=eta,
                                distance=distance,
                                physical_error_rate=physicall_error_rate,
                                logical_error_rate=p_L,
                                standard_deviation=sigma,
                                logical_errors=errors,
                                shots=shots
                            )

                            result.append(sample_point)

                            if(code_type == "css"):
                                css_sample_points.append(sample_point)
                            else:
                                xzzx_sample_points.append(sample_point)

                            progress.update(p_error_task, advance=1)

                        if(code_type == "css"):
                            css_results[f"{code_type}_d{distance}"] = result
                        else:
                            xzzx_results[f"{code_type}_d{distance}"] = result
                        
                        progress.update(distance_task, advance=1)
                        progress.remove_task(p_error_task)
                    
                    progress.update(code_type_task, advance=1)
                    progress.remove_task(distance_task)

                # FSS fitting
                css_fit = estimate_threshold(css_sample_points)
                xzzx_fit = estimate_threshold(xzzx_sample_points)

                # Draw results
                mpl.rcParams["font.family"] = "serif"
                mpl.rcParams["font.size"] = 12
                mpl.rcParams["lines.linewidth"] = 1.5

                css_colors = plt.cm.plasma(np.linspace(0.2, 0.8, len(distances) + 1))
                xzzx_colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(distances) + 1))
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
                        (css_fit, css_colors[len(distances)], "css"),
                        (xzzx_fit, xzzx_colors[len(distances)], "xzzx")
                    ):
                    ax.axvline(x=fit.p_th, color=color, linestyle="--",
                               label=f"{name} p_th = {fit.p_th:.4f} ± {fit.p_th_err:.4f}")
                    ax.axvspan(fit.p_th - fit.p_th_err, fit.p_th + fit.p_th_err,
                               color=color, alpha=0.1)

                ax.set_xscale('linear')
                ax.set_yscale('log')

                ax.set_xlabel('Physical Error Rate')
                ax.set_ylabel('Logical Error Rate (error bars: 1σ Wilson)')
                ax.set_title(f'Rotated CSS Surface Code vs XZZX Code Simulation Results with η = {eta}')

                ax.grid(True, which="both", alpha=0.5)
                ax.legend()

                plt.savefig(f"{args.outdir}/result_{eta}.png")
                plt.close(fig)

                # Draw data-collapse figure: rescaled axis x = (p - p_th) d^(1/nu)
                draw_collapse(
                    f"{args.outdir}/collapse_{eta}.png", eta,
                    (css_results, css_fit, css_colors, "CSS"),
                    (xzzx_results, xzzx_fit, xzzx_colors, "XZZX")
                )
            
                progress.update(eta_task, advance=1)
                progress.remove_task(code_type_task)
    
    # Generate circuit diagrams
    os.makedirs(f"{args.outdir}/diagrams", exist_ok=True)
    for distance in distances:
        eta = etas[0]
        p = physical_error_rates[0]

        codes: List[QecCode] = [build_rotated_surface_code(distance), build_xzzx_code(distance)]

        for code in codes:
            circuit = CodeCapacityCircuitBuilder(code, p, eta).build()

            detslice = circuit.diagram("detslice-svg")
            timeline = circuit.diagram("timeline-svg")

            with open(f"{args.outdir}/diagrams/{code.name}_detslice.svg", "w") as f:
                f.write(str(detslice))
            with open(f"{args.outdir}/diagrams/{code.name}_timeline.svg", "w") as f:
                f.write(str(timeline))
