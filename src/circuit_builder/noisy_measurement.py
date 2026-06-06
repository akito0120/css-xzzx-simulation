from typing import Dict, List, Optional
from code_builder import QecCode, Coord
from .base import BaseCircuitBuilder

class NoisyMeasurementCircuitBuilder(BaseCircuitBuilder):
    # Shared base for the noisy-measurement models (phenomenological, circuit-level).
    # These run multiple rounds because a single measurement cannot be trusted
    # and close the time boundary with detectors reconstructed from the final data readout.

    rounds: int
    p_meas: float

    def __init__(self, code: QecCode, p: float, eta: float,
                 rounds: Optional[int] = None, p_meas: Optional[float] = None):
        super().__init__(code, p, eta)
        self.rounds = code.distance if rounds is None else rounds
        self.p_meas = p if p_meas is None else p_meas

    def initial_boundary_detectors(self):
        # Mirror of final_boundary_detectors at the bottom time boundary
        # Use if data prep is noisy
        prep_basis: Dict[Coord, str] = {
            q: ("Z" if q in self.code.hset else "X") for q in self.code.data_qubits
        }
        eligible: List[Coord] = [
            anc for anc, legs in self.code.stabilizers.items()
            if all(pauli == prep_basis[q] for q, pauli in legs.items())
        ]
        for ancilla in eligible:
            d0 = self.ancilla_record[(ancilla, 0)]
            self.circuit.append("DETECTOR", [self.rel(d0)], [ancilla[0], ancilla[1], 0])

    def final_boundary_detectors(self, data_record: Dict[Coord, int]):
        # Per-qubit preparation/measurement basis for X-memory
        # |+> (basis X) off hset, |0> (basis Z) on hset
        prep_basis: Dict[Coord, str] = {
            q: ("Z" if q in self.code.hset else "X") for q in self.code.data_qubits
        }

        # "eligible" stabilizers = checks deterministic under that preparation
        # Every leg matches the prep basis
        # Checks reconstructable from the final data readout
        eligible: List[Coord] = [
            anc for anc, legs in self.code.stabilizers.items()
            if all(pauli == prep_basis[q] for q, pauli in legs.items())
        ]

        # Close the time boundary: reconstructed check vs. last round
        last = self.current_round
        for ancilla in eligible:
            targets = [self.rel(data_record[dc])
                       for dc in self.code.stabilizers[ancilla]]
            targets.append(self.rel(self.ancilla_record[(ancilla, last)]))
            self.circuit.append("DETECTOR", targets,
                                [ancilla[0], ancilla[1], last + 1])
