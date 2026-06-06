from typing import Tuple, List, Dict, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from code_builder import QecCode, Coord

CGATE = {"X": "CX", "Y": "CY", "Z": "CZ"}

# Distance-preserving 4-step CNOT schedule
# Maps a leg's relative offset (dx, dy) -> time step
# Adding a Z-memory experiment would require revisiting this (the standard X/Z-transposed schedule).
STEP_OF_OFFSET = {(-1, -1): 0, (+1, -1): 1, (-1, +1): 2, (+1, +1): 3}

def build_cnot_schedule(code: "QecCode") -> Dict["Coord", List[Optional["Coord"]]]:
    # Assign each stabilizer's legs to one of 4 parallel time steps
    schedule: Dict["Coord", List[Optional["Coord"]]] = {}
    for (x, y), legs in code.stabilizers.items():
        steps: List[Optional["Coord"]] = [None, None, None, None]
        for (dx_coord, dy_coord) in legs:
            steps[STEP_OF_OFFSET[(dx_coord - x, dy_coord - y)]] = (dx_coord, dy_coord)
        schedule[(x, y)] = steps
    return schedule

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
