import os
from code_builder import build_rotated_surface_code, build_xzzx_code
from circuit_builder import build_circuit
from simulation import estimate_logical_error_rate
import numpy as np
from typing import Dict, List
import matplotlib.pyplot as plt
import matplotlib as mpl

if __name__ == "__main__":
    etas = [0.5, 3, 10]
    distances = [3, 5, 7, 9, 11]
    physical_error_rates = ps = list(np.linspace(0.05, 0.30, 9))

    for eta in etas:
        results: Dict[str, np.ndarray] = {}

        # Simulation
        for code_type in ["surface", "xzzx"]:
            for distance in distances:
                result: List[float] = []

                for physicall_error_rate in physical_error_rates:
                    code = build_rotated_surface_code(distance) if code_type == "surface" else build_xzzx_code(distance)
                    circuit = build_circuit(code, physicall_error_rate, eta)

                    print(f"Simulating {code.name} with physical error rate = {physicall_error_rate}")

                    svg =  circuit.diagram("detslice-svg")

                    outdir = "results"
                    os.makedirs(outdir, exist_ok=True)
                    with open(f"{outdir}/{code.name}.svg", "w") as f:
                        f.write(str(svg))

                    logical_error_rate = estimate_logical_error_rate(circuit, shots=100_000)
                    print(f"Result: logical error rate = {logical_error_rate}")

                    result.append(logical_error_rate)

                results[f"{code_type}_d{distance}"] = np.array(result)

        # Draw results
        mpl.rcParams["font.family"] = "serif"
        mpl.rcParams["font.size"] = 12
        mpl.rcParams["lines.linewidth"] = 1.5

        colors = plt.cm.plasma(
            np.linspace(0.2, 0.8, 2 * len(distances))
        )

        fig, ax = plt.subplots(figsize=(12, 8), dpi=600)
        for i, (label, y_data) in enumerate(results.items()):
            ax.plot(physical_error_rates, y_data, marker='o', label=label, color=colors[i])

        ax.set_xscale('log')
        ax.set_yscale('log')

        ax.set_xlabel('Physical Error Rate')
        ax.set_ylabel('Logical Error Rate')
        ax.set_title(f'Rotated Surface Code vs XZZX Code Simulation Results with η = {eta}')

        ax.grid(True, which="both", alpha=0.5)
        ax.legend()

        plt.savefig(f"{outdir}/result_{eta}.png")
