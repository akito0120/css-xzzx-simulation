from typing import Tuple, List

CGATE = {"X": "CX", "Y": "CY", "Z": "CZ"}

def biased_pauli_rates(p: float, eta: float) -> Tuple[float, float, float]:
    if eta == float("inf"):
        return 0.0, 0.0, p
    px = py = p / (2.0 * (1.0 + eta))
    pz = p * eta / (1.0 + eta)
    return px, py, pz


def biased_two_qubit_rates(p: float, eta: float) -> List[float]:
    # Biased correlated 2-qubit Pauli channel for the two-qubit-gate error
    # Order must match stim PAULI_CHANNEL_2
    order = ["IX", "IY", "IZ", "XI", "XX", "XY", "XZ", "YI", "YX", "YY", "YZ", "ZI", "ZX", "ZY", "ZZ"]
    high = {"IZ", "ZI", "ZZ"}
    if eta == float("inf"):
        return [p / 3.0 if ab in high else 0.0 for ab in order]
    p_high = p * eta / (3.0 * (1.0 + eta))
    p_low = p / (12.0 * (1.0 + eta))
    return [p_high if ab in high else p_low for ab in order]
