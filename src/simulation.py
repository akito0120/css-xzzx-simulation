import stim
import pymatching
import numpy as np

def estimate_logical_error_rate(circuit: stim.Circuit, shots: int = 100_000):
    dem = circuit.detector_error_model(decompose_errors=True)
    matching = pymatching.Matching.from_detector_error_model(dem)
    sampler = circuit.compile_detector_sampler()

    det, obs = sampler.sample(shots, separate_observables=True)
    pred = matching.decode_batch(det)
    mism = np.any(pred != obs, axis=1)
    total_errors = int(np.count_nonzero(mism))

    p_L = total_errors / shots
    return p_L

