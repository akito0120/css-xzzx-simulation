from dataclasses import dataclass
from typing import List
import numpy as np
from scipy.optimize import curve_fit

@dataclass(frozen=True)
class SamplePoint:
    eta: float
    distance: int
    physical_error_rate: float
    logical_error_rate: float
    standard_deviation: float
    logical_errors: int
    shots: int

def fss(X, p_th, nu, a, b, c):
    # Approximation of the scaling function by a second-order polynomial
    # TODO: take the error bias into account
    p, d = X
    x = (p - p_th) * d ** (1.0 / nu)
    return a + (b * x) + (c * x * x)

def estimate_threshold(sample_points: List[SamplePoint]):
    # FSS fitting with least squares
    ps = np.array([sample.physical_error_rate for sample in sample_points])
    ds = np.array([sample.distance for sample in sample_points])
    p_Ls = np.array([sample.logical_error_rate for sample in sample_points])
    sigs = np.array([sample.standard_deviation for sample in sample_points])

    # Initial guess (p_th, nu, a, b, c)
    p0 = [float(np.median(ps)), 1.5, float(np.median(p_Ls)), 0.0, 0.0]

    popt, pcov = curve_fit(fss, (ps, ds), p_Ls, sigma=sigs, p0=p0, maxfev=20000)
    p_th = popt[0]
    return p_th
