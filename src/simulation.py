import sinter
import numpy as np
from typing import Tuple
from code_builder import build_code
from circuit_builder import CircuitLevelCircuitBuilder
from config import ETAS, CODE_TYPES, DISTANCES, P_STEP, P_WINDOWS
from rich.status import Status
from joblib import Parallel, delayed
import pandas as pd

def physical_error_rates(eta: float, code_type: str) -> list[float]:
    if (eta, code_type) not in P_WINDOWS:
        raise KeyError(f"No p-window configured for (eta={eta!r}, code_type={code_type!r})")
    p_min, p_max = P_WINDOWS[(eta, code_type)]
    return list(np.arange(p_min, p_max + P_STEP * 1e-9, P_STEP))

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

def generate_params():
    for eta in ETAS:
        for code_type in CODE_TYPES:
            ps = physical_error_rates(eta, code_type)
            for distance in DISTANCES:
                for p in ps:
                    yield (eta, code_type, distance, p)

def generate_tasks():
    def build_task(param: tuple[float, str, int, float]) -> sinter.Task:
        eta, code_type, distance, p = param
        code = build_code(code_type, distance)
        circuit = CircuitLevelCircuitBuilder(code, p, eta).build()
        dem = circuit.detector_error_model(
            decompose_errors=True,
            approximate_disjoint_errors=True
        )
        return sinter.Task(
            circuit=circuit,
            detector_error_model=dem,
            json_metadata={
                "code": code_type,
                "eta": str(eta),
                "d": distance,
                "p": p,
            },
        )

    for param in generate_params():
        yield build_task(param)

def sweep(max_shots: int, target_errors: int, num_workers: int) -> pd.DataFrame:
    with Status("Sweeping", spinner="arc"):
        stats = sinter.collect(
            num_workers=num_workers,
            tasks=generate_tasks(),
            decoders=["pymatching"],
            max_shots=max_shots,
            max_errors=target_errors,
        )

        rows: list[dict] = list()
        for stat in stats:
            meta, shots, errors = stat.json_metadata, stat.shots, stat.errors
            pl = errors / shots if shots > 0 else 0.0
            low, high = wilson_interval(errors, shots, z=1.0)
            rows.append({
                "code": meta["code"],
                "eta": float(meta["eta"]),
                "d": meta["d"],
                "p": meta["p"],
                "pl": pl,
                "sigma": (high - low) / 2.0,
                "errors": errors,
                "shots": shots,
            })

        print("☑ Sweeping completed")
        return pd.DataFrame(rows)

def verify_distance_preservation():
    def verify(eta: float, code_type: str, distance: int):
        code = build_code(code_type, distance)
        circuit = CircuitLevelCircuitBuilder(code, 0.1, eta).build()
        err = circuit.shortest_graphlike_error()
        if len(err) != distance:
            raise AssertionError(f"Graphlike fault distance {len(err)} != code distance {distance}")

    with Status("Verifying circuits", spinner="arc"):
        Parallel(n_jobs=-1)(
            delayed(verify)(eta, code_type, distance)
            for eta in ETAS
            for code_type in CODE_TYPES
            for distance in DISTANCES
        )
        print("☑ Circuit verification completed")
