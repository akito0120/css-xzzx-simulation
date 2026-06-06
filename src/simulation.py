import stim
import pymatching
import numpy as np
from typing import Tuple

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
):
    dem = circuit.detector_error_model(decompose_errors=True, approximate_disjoint_errors=True)
    matching = pymatching.Matching.from_detector_error_model(dem)
    sampler = circuit.compile_detector_sampler()

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
