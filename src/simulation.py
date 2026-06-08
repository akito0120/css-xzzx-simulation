import stim
import pymatching
import numpy as np
import zlib
from typing import Optional, Tuple, List
from rich.progress import Progress, BarColumn, TextColumn
from code_builder import build_code
from circuit_builder import CircuitLevelCircuitBuilder
from config import ETAS, CODE_TYPES, DISTANCES, PHYSICAL_ERROR_RATES

def wilson_interval(errors: int, shots: int, z: float = 1.96) -> Tuple[float, float]:
    # Wilson interval for a binomial proportion
    if shots <= 0:
        return 0.0, 0.0
    
    p_hat = errors / shots
    denom = 1.0 + z * z / shots
    nom = (p_hat + (z * z / (2.0 * shots)))

    center = nom / denom
    half = (z / denom) * np.sqrt(p_hat * (1.0 - p_hat) / shots + z * z / (4.0 * shots * shots))

    low = max(0.0, center - half)
    high = min(1.0, center + half)

    return low, high

def estimate_logical_error_rate(
    circuit: stim.Circuit,
    max_shots: int = 2_000_000,
    target_errors: int = 200,
    batch_size: int = 100_000,
    seed: Optional[int] = None,
):
    dem = circuit.detector_error_model(decompose_errors=True, approximate_disjoint_errors=True)
    matching = pymatching.Matching.from_detector_error_model(dem)
    sampler = circuit.compile_detector_sampler(seed=seed)

    total_errors = 0
    total_shots = 0

    # Adaptive Monte Carlo sampling
    while total_shots < max_shots and total_errors < target_errors:
        batch = min(batch_size, max_shots - total_shots)
        det, obs = sampler.sample(batch, separate_observables=True)
        pred = matching.decode_batch(det)
        mism = np.any(pred != obs, axis=1)
        total_errors += int(np.count_nonzero(mism))
        total_shots += batch

    # Estimated logical error rate
    p_L = total_errors / total_shots
    low, high = wilson_interval(total_errors, total_shots, z=1.0)
    sigma = (high - low) / 2.0

    return p_L, sigma, total_errors, total_shots

def point_seed(base_seed: int, code_type: str, distance: int, eta: float, p: float) -> int:
    key = f"{base_seed}|{code_type}|{distance}|{eta}|{p!r}".encode()
    return zlib.crc32(key) & 0xFFFFFFFF

def sweep(max_shots: int, target_errors: int, batch_size: int, seed: int) -> List[dict]:
    # Run the full sweep and return one CSV-ready row per sample point
    rows: List[dict] = []
    with Progress(
            TextColumn("{task.description}"),
            BarColumn(bar_width=30),
            TextColumn("{task.completed}/{task.total}")
        ) as progress:

        eta_task = progress.add_task("Simulation", total=len(ETAS))

        for eta in ETAS:
            code_type_task = progress.add_task(f"Simulating with eta = {eta}", total=len(CODE_TYPES))

            for code_type in CODE_TYPES:
                distance_task = progress.add_task(f"Simulating {code_type}", total=len(DISTANCES))

                for distance in DISTANCES:
                    p_error_task = progress.add_task(f"Simulating with d = {distance}", total=len(PHYSICAL_ERROR_RATES))

                    for physical_error_rate in PHYSICAL_ERROR_RATES:
                        code = build_code(code_type, distance)
                        circuit = CircuitLevelCircuitBuilder(code, physical_error_rate, eta).build()
                        p_seed = point_seed(seed, code_type, distance, eta, physical_error_rate)

                        p_L, sigma, errors, shots = estimate_logical_error_rate(
                            circuit,
                            max_shots=max_shots,
                            target_errors=target_errors,
                            batch_size=batch_size,
                            seed=p_seed,
                        )
                        rows.append({
                            "eta": eta,
                            "distance": distance,
                            "physical_error_rate": physical_error_rate,
                            "logical_error_rate": p_L,
                            "standard_deviation": sigma,
                            "logical_errors": errors,
                            "shots": shots,
                            "code_type": code_type,
                            "seed": p_seed
                        })
                        progress.update(p_error_task, advance=1)

                    progress.update(distance_task, advance=1)
                    progress.remove_task(p_error_task)

                progress.update(code_type_task, advance=1)
                progress.remove_task(distance_task)

            progress.update(eta_task, advance=1)
            progress.remove_task(code_type_task)

    return rows

def verify_distance_preservation():
    with Progress(
            TextColumn("{task.description}"),
            BarColumn(bar_width=30),
            TextColumn("{task.completed}/{task.total}")
        ) as progress:

        task = progress.add_task("Verifying distance preservation", total=len(ETAS) * len(CODE_TYPES) * len(DISTANCES))
        for eta in ETAS:
            for code_type in CODE_TYPES:
                for distance in DISTANCES:
                    code = build_code(code_type, distance)
                    circuit = CircuitLevelCircuitBuilder(code, 0.1, eta).build()
                    err = circuit.shortest_graphlike_error()
                    if len(err) != distance:
                        raise AssertionError(f"Graphlike fault distance {len(err)} != code distance {distance}")
                    progress.update(task, advance=1)
