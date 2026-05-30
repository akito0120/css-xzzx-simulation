import os
from code_builder import build_rotated_surface_code, build_xzzx_code
from circuit_builder import build_circuit
from simulation import estimate_logical_error_rate
import numpy as np
from typing import Dict, List
import matplotlib.pyplot as plt
import matplotlib as mpl
from threshold import SamplePoint, estimate_p_th

if __name__ == "__main__":
    etas = [0.5, 3, 10]
    distances = [3, 5, 7, 9, 11]
    physical_error_rates = ps = list(np.linspace(0.05, 0.50, 12))

    for eta in etas:
        # Sample points for FSS fitting
        css_sample_points: List[SamplePoint] = []
        xzzx_sample_points: List[SamplePoint] = []

        # Simulated logical errors: code type + distance -> logical errors
        css_results: Dict[str, np.ndarray] = {}
        xzzx_results: Dict[str, np.ndarray] = {}

        # Simulation
        for code_type in ["css", "xzzx"]:
            for distance in distances:
                result: List[float] = []

                for physicall_error_rate in physical_error_rates:
                    code = build_rotated_surface_code(distance) if code_type == "css" else build_xzzx_code(distance)
                    circuit = build_circuit(code, physicall_error_rate, eta)

                    print(f"Simulating {code.name} with p = {physicall_error_rate}, eta = {eta}")

                    svg =  circuit.diagram("detslice-svg")

                    outdir = "results"
                    os.makedirs(outdir, exist_ok=True)
                    with open(f"{outdir}/{code.name}.svg", "w") as f:
                        f.write(str(svg))

                    p_L, sigma = estimate_logical_error_rate(circuit, shots=100_000)
                    print(f"Result: p_L = {p_L}, sigma = {sigma}")

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

                if(code_type == "css"):
                    css_results[f"{code_type}_d{distance}"] = np.array(result)
                else:
                    xzzx_results[f"{code_type}_d{distance}"] = np.array(result)

        # FSS fitting
        css_p_th = estimate_p_th(css_sample_points)
        xzzx_p_th = estimate_p_th(xzzx_sample_points)

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
        ax.axvline(x = css_p_th, label=f"CSS p_th = {css_p_th:.3f}", color=css_colors[len(distances)], linestyle="--")
        ax.axvline(x = xzzx_p_th, label=f"XZZX p_th = {xzzx_p_th:.3f}", color=xzzx_colors[len(distances)], linestyle="--")

        ax.set_xscale('linear')
        ax.set_yscale('log')

        ax.set_xlabel('Physical Error Rate')
        ax.set_ylabel('Logical Error Rate')
        ax.set_title(f'Rotated CSS Surface Code vs XZZX Code Simulation Results with η = {eta}')

        ax.grid(True, which="both", alpha=0.5)
        ax.legend()

        plt.savefig(f"{outdir}/result_{eta}.png")
