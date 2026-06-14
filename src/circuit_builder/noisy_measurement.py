import stim
from typing import Callable, Dict, List, Optional
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

    def repeat_rounds(self, n: int, round_body: Callable[[], stim.Circuit]) -> stim.Circuit:
        # Emit n identical rounds as a stim REPEAT block.
        if n <= 0:
            return stim.Circuit()
        n_ancilla = len(self.ancilla_order)

        # Save values before repeating rounds
        rc0 = self.record_counter
        cr0 = self.current_round
        cb0 = self.coord_base

        # Build one representative iteration.
        # Relative record offsets rel() and the SHIFT_COORDS-relative time coordinate are identical for every iteration.
        body = stim.Circuit()
        body.append("SHIFT_COORDS", [], [0, 0, 1])
        self.coord_base += 1
        self.current_round += 1
        body += round_body()
        body.append("TICK")

        # Fix up bookkeeping to reflect the full n-round unroll
        # so post-rounds helpers resolve the correct relative offsets
        self.current_round = cr0 + n
        self.record_counter = rc0 + n * n_ancilla
        self.coord_base = cb0 + n
        for i, ancilla in enumerate(self.ancilla_order):
            self.ancilla_record[(ancilla, self.current_round)] = self.record_counter - n_ancilla + i

        return body * n

    def initial_boundary_detectors(self) -> stim.Circuit:
        # Mirror of final_boundary_detectors at the bottom time boundary
        # Use if data prep is noisy
        prep_basis: Dict[Coord, str] = {
            q: ("Z" if q in self.code.hset else "X") for q in self.code.data_qubits
        }
        eligible: List[Coord] = [
            anc for anc, legs in self.code.stabilizers.items()
            if all(pauli == prep_basis[q] for q, pauli in legs.items())
        ]
        circuit = stim.Circuit()
        for ancilla in eligible:
            d0 = self.ancilla_record[(ancilla, 0)]
            circuit.append("DETECTOR", [self.rel(d0)], [ancilla[0], ancilla[1], 0 - self.coord_base])
        return circuit

    def final_boundary_detectors(self, data_record: Dict[Coord, int]) -> stim.Circuit:
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
        circuit = stim.Circuit()
        last = self.current_round
        for ancilla in eligible:
            targets = [self.rel(data_record[dc])
                       for dc in self.code.stabilizers[ancilla]]
            targets.append(self.rel(self.ancilla_record[(ancilla, last)]))
            circuit.append("DETECTOR", targets, [ancilla[0], ancilla[1], last + 1])
            
        return circuit
