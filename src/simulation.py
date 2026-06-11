import sinter
import numpy as np
from typing import Tuple
from code_builder import build_code
from circuit_builder import CircuitLevelCircuitBuilder
from config import ETAS, CODE_TYPES, DISTANCES, P_STEP, P_WINDOWS
from rich.status import Status
from joblib import Parallel, delayed

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

def build_tasks() -> list[sinter.Task]:
    task_params: list[tuple[float, str, int, float]] = list()
    for eta in ETAS:
        for code_type in CODE_TYPES:
            ps = physical_error_rates(eta, code_type)
            for distance in DISTANCES:
                for p in ps:
                    task_params.append((eta, code_type, distance, p))

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
                "code_type": code_type,
                "eta": str(eta),
                "distance": distance,
                "physical_error_rate": p,
            },
        )
    
    with Status("Building tasks", spinner="arc"):
        tasks = Parallel(n_jobs=-1)(delayed(build_task)(param) for param in task_params)
        print("☑ Tasks built successfully")
        return tasks

def sweep(max_shots: int, target_errors: int, num_workers: int) -> list[dict]:
    tasks = build_tasks()

    with Status("Sweeping", spinner="arc"):
        rows: list[dict] = list()
        stats = sinter.collect(
            num_workers=num_workers,
            tasks=tasks,
            decoders=["pymatching"],
            max_shots=max_shots,
            max_errors=target_errors
        )
    
        for stat in stats:
            meta = stat.json_metadata
            shots = stat.shots
            errors = stat.errors
            p_L = errors / shots if shots > 0 else 0.0
            low, high = wilson_interval(errors, shots, z=1.0)
            sigma = (high - low) / 2.0
            rows.append({
                "code_type": meta["code_type"],
                "eta": float(meta["eta"]),
                "distance": meta["distance"],
                "physical_error_rate": meta["physical_error_rate"],
                "logical_error_rate": p_L,
                "standard_deviation": sigma,
                "logical_errors": errors,
                "shots": shots,
            })

        print("☑ Sweeping completed")
        return rows

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
