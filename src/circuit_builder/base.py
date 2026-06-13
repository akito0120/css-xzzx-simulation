import stim
from typing import Tuple, Dict, List
from code_builder import QecCode, Coord
from .shared import CGATE


class BaseCircuitBuilder:
    # Common machinery shared by every noise-model circuit builder.
    # Memory basis: X
    # Subclasses implement build() and reuse the helpers below

    code: QecCode
    p: float
    eta: float

    circuit: stim.Circuit

    record_counter: int
    current_round: int
    coord_base: int
    ancilla_record: Dict[Tuple[Coord, int], int]
    ancilla_order: List[Coord]

    def __init__(self, code: QecCode, p: float, eta: float):
        self.code = code
        self.p = p
        self.eta = eta

        self.record_counter = 0
        self.current_round = 0
        self.coord_base = 0
        self.ancilla_record = {}

    def rel(self, abs_idx: int) -> stim.GateTarget:
        # Return the relative index of a measurement result
        return stim.target_rec(abs_idx - self.record_counter)

    def data_basis_partition(self):
        # Split data qubits by memory basis: X basis off hset, Z basis on hset.
        x_basis = [(coord, idx) for coord, idx in self.code.data_qubits.items()
                   if coord not in self.code.hset]
        z_basis = [(coord, idx) for coord, idx in self.code.data_qubits.items()
                   if coord in self.code.hset]
        return x_basis, z_basis

    def prep_data(self):
        # Prepare each data qubit in its memory basis
        x_basis, z_basis = self.data_basis_partition()
        if x_basis:
            self.circuit.append("RX", [idx for _, idx in x_basis])
        if z_basis:
            self.circuit.append("R", [idx for _, idx in z_basis])

    def init_qubit_coords(self):
        for coord, idx in self.code.data_qubits.items():
            self.circuit.append("QUBIT_COORDS", [idx], [coord[0], coord[1]])
        for coord, idx in self.code.ancilla_qubits.items():
            self.circuit.append("QUBIT_COORDS", [idx], [coord[0], coord[1]])

    def syndrome_meas(self, flip: float = 0.0):
        # 1 round of syndrome extraction.
        # flip > 0 injects a measurement readout error with that probability.
        for ancilla in self.ancilla_order:
            self.circuit.append("RX", [self.code.ancilla_qubits[ancilla]])

        for ancilla in self.ancilla_order:
            ancilla_idx = self.code.ancilla_qubits[ancilla]
            for dcoord, pauli in self.code.stabilizers[ancilla].items():
                self.circuit.append(CGATE[pauli],
                                    [ancilla_idx, self.code.data_qubits[dcoord]])

        for ancilla in self.ancilla_order:
            ancilla_idx = self.code.ancilla_qubits[ancilla]
            if flip > 0.0:
                self.circuit.append("MX", [ancilla_idx], flip)
            else:
                self.circuit.append("MX", [ancilla_idx])
            self.ancilla_record[(ancilla, self.current_round)] = self.record_counter
            self.record_counter += 1

    def consecutive_round_detectors(self):
        # Compare every ancilla's syndrome in the current round against the previous round.
        # A change between consecutive rounds flags an error that occurred in between.
        for ancilla in self.ancilla_order:
            d_now = self.ancilla_record[(ancilla, self.current_round)]
            d_prev = self.ancilla_record[(ancilla, self.current_round - 1)]
            self.circuit.append(
                "DETECTOR",
                [self.rel(d_now), self.rel(d_prev)],
                [ancilla[0], ancilla[1], self.current_round - self.coord_base],
            )

    def data_readout(self, flip: float = 0.0) -> Dict[Coord, int]:
        # Final data readout: measure each data qubit in its memory basis
        # Returns coord -> absolute measurement index.
        x_basis, z_basis = self.data_basis_partition()
        data_record: Dict[Coord, int] = {}
        for gate, group in (("MX", x_basis), ("M", z_basis)):
            if not group:
                continue
            idxs = [idx for _, idx in group]
            if flip > 0.0:
                self.circuit.append(gate, idxs, flip)
            else:
                self.circuit.append(gate, idxs)
            for coord, _ in group:
                data_record[coord] = self.record_counter
                self.record_counter += 1
        return data_record

    def define_observable(self, data_record: Dict[Coord, int]):
        # The logical-X observable: parity of the final data measurements on logical_x.
        observable_targets = [self.rel(data_record[dc]) for dc in self.code.logical_x]
        self.circuit.append("OBSERVABLE_INCLUDE", observable_targets, 0)
