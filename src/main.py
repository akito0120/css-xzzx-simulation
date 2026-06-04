import os
import numpy as np
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
from typing import Dict, List
from rich.progress import Progress, BarColumn, TextColumn

from code_builder import QecCode, build_rotated_surface_code, build_xzzx_code
from circuit_builder import CodeCapacityCircuitBuilder, PhenomenologicalCircuitBuilder
from simulation import estimate_logical_error_rate
from threshold import SamplePoint, estimate_threshold

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--threshold", action="store_true", default=False)
    ap.add_argument("--shots", default=100_000)
    ap.add_argument("--diagramonly", action="store_true", default=False)
    args = ap.parse_args()

    etas = [3]
    code_types = ["css", "xzzx"]
    distances = [3, 5, 7]
    physical_error_rates = ps = list(np.linspace(0.01, 0.20, 12))

    os.makedirs(f"{args.outdir}", exist_ok=True)

    if args.diagramonly == False:
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

                # Simulated logical errors: code type + distance -> logical errors
                css_results: Dict[str, np.ndarray] = {}
                xzzx_results: Dict[str, np.ndarray] = {}

                # Simulation
                for code_type in code_types:
                    distance_task = progress.add_task(f"Simulating {code_type}", total=len(distances))

                    for distance in distances:
                        p_error_task = progress.add_task(f"Simulating with d = {distance}", total=len(physical_error_rates))
                        result: List[float] = []

                        for physicall_error_rate in physical_error_rates:

                            code = build_rotated_surface_code(distance) if code_type == "css" else build_xzzx_code(distance)
                            circuit = PhenomenologicalCircuitBuilder(code, physicall_error_rate, eta).build()
                            p_L, sigma = estimate_logical_error_rate(circuit, shots=int(args.shots))

                            result.append(p_L)

                            sample_point = SamplePoint(
                                eta=eta, 
                                distance=distance, 
                                physical_error_rate=physicall_error_rate,
                                logical_error_rate=p_L,
                                standard_deviation=sigma
                            )

                            if(code_type == "css"):
                                css_sample_points.append(sample_point)
                            else:
                                xzzx_sample_points.append(sample_point)
                            
                            progress.update(p_error_task, advance=1)

                        if(code_type == "css"):
                            css_results[f"{code_type}_d{distance}"] = np.array(result)
                        else:
                            xzzx_results[f"{code_type}_d{distance}"] = np.array(result)
                        
                        progress.update(distance_task, advance=1)
                        progress.remove_task(p_error_task)
                    
                    progress.update(code_type_task, advance=1)
                    progress.remove_task(distance_task)

                # FSS fitting
                css_p_th = estimate_threshold(css_sample_points)
                xzzx_p_th = estimate_threshold(xzzx_sample_points)

                # Draw results
                mpl.rcParams["font.family"] = "serif"
                mpl.rcParams["font.size"] = 12
                mpl.rcParams["lines.linewidth"] = 1.5

                css_colors = plt.cm.plasma(np.linspace(0.2, 0.8, len(distances) + 1))
                xzzx_colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(distances) + 1))
                fig, ax = plt.subplots(figsize=(12, 8), dpi=600)

                # Draw logical error rates
                for i, (label, y_data) in enumerate(css_results.items()):
                    ax.plot(physical_error_rates, y_data, marker='o', label=label, color=css_colors[i])
                for i, (label, y_data) in enumerate(xzzx_results.items()):
                    ax.plot(physical_error_rates, y_data, marker='o', label=label, color=xzzx_colors[i])

                # Draw thresholds
                if(args.threshold):
                    ax.axvline(x = css_p_th, label=f"CSS p_th = {css_p_th:.3f}", color=css_colors[len(distances)], linestyle="--")
                    ax.axvline(x = xzzx_p_th, label=f"XZZX p_th = {xzzx_p_th:.3f}", color=xzzx_colors[len(distances)], linestyle="--")

                ax.set_xscale('linear')
                ax.set_yscale('log')

                ax.set_xlabel('Physical Error Rate')
                ax.set_ylabel('Logical Error Rate')
                ax.set_title(f'Rotated CSS Surface Code vs XZZX Code Simulation Results with η = {eta}')

                ax.grid(True, which="both", alpha=0.5)
                ax.legend()

                plt.savefig(f"{args.outdir}/result_{eta}.png")
            
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
