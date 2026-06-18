import stim
from typing import Tuple, Dict, List
from code_builder import QecCode, Coord
from .shared import CGATE


class BaseCircuitBuilder:
    # Common machinery shared by every noise-model circuit builder.
    # Subclasses implement build() and reuse the helpers below

    code: QecCode
    p: float
    eta: float
    basis: str

    record_counter: int
    current_round: int
    coord_base: int
    ancilla_record: Dict[Tuple[Coord, int], int]
    ancilla_order: List[Coord]

    def __init__(self, code: QecCode, p: float, eta: float, basis: str = "X"):
        if basis not in ("X", "Z"):
            raise ValueError(f"basis must be 'X' or 'Z', got {basis!r}")
        self.code = code
        self.p = p
        self.eta = eta
        self.basis = basis

        self.record_counter = 0
        self.current_round = 0
        self.coord_base = 0
        self.ancilla_record = {}

    def rel(self, abs_idx: int) -> stim.GateTarget:
        # Return the relative index of a measurement result
        return stim.target_rec(abs_idx - self.record_counter)

    def prep_basis_of(self, q: Coord) -> str:
        # Per-qubit prep/measurement basis
        conj = {"X": "Z", "Z": "X"}
        return conj[self.basis] if q in self.code.hset else self.basis

    def data_basis_partition(self):
        # Split data qubits by their prep/measurement basis (X-prepared vs Z-prepared)
        x_basis = [(coord, idx) for coord, idx in self.code.data_qubits.items()
                   if self.prep_basis_of(coord) == "X"]
        z_basis = [(coord, idx) for coord, idx in self.code.data_qubits.items()
                   if self.prep_basis_of(coord) == "Z"]
        return x_basis, z_basis

    def prep_data(self) -> stim.Circuit:
        # Prepare each data qubit in its memory basis
        circuit = stim.Circuit()
        x_basis, z_basis = self.data_basis_partition()
        if x_basis:
            circuit.append("RX", [idx for _, idx in x_basis])
        if z_basis:
            circuit.append("R", [idx for _, idx in z_basis])
        return circuit

    def init_qubit_coords(self) -> stim.Circuit:
        circuit = stim.Circuit()
        for coord, idx in self.code.data_qubits.items():
            circuit.append("QUBIT_COORDS", [idx], [coord[0], coord[1]])
        for coord, idx in self.code.ancilla_qubits.items():
            circuit.append("QUBIT_COORDS", [idx], [coord[0], coord[1]])
        return circuit

    def syndrome_meas(self, flip: float = 0.0) -> stim.Circuit:
        # 1 round of syndrome extraction.
        # flip > 0 injects a measurement readout error with that probability.
        circuit = stim.Circuit()
        for ancilla in self.ancilla_order:
            circuit.append("RX", [self.code.ancilla_qubits[ancilla]])

        for ancilla in self.ancilla_order:
            ancilla_idx = self.code.ancilla_qubits[ancilla]
            for dcoord, pauli in self.code.stabilizers[ancilla].items():
                circuit.append(CGATE[pauli],
                               [ancilla_idx, self.code.data_qubits[dcoord]])

        for ancilla in self.ancilla_order:
            ancilla_idx = self.code.ancilla_qubits[ancilla]
            if flip > 0.0:
                circuit.append("MX", [ancilla_idx], flip)
            else:
                circuit.append("MX", [ancilla_idx])
            self.ancilla_record[(ancilla, self.current_round)] = self.record_counter
            self.record_counter += 1
        return circuit

    def consecutive_round_detectors(self) -> stim.Circuit:
        # Compare every ancilla's syndrome in the current round against the previous round.
        # A change between consecutive rounds flags an error that occurred in between.
        circuit = stim.Circuit()
        for ancilla in self.ancilla_order:
            d_now = self.ancilla_record[(ancilla, self.current_round)]
            d_prev = self.ancilla_record[(ancilla, self.current_round - 1)]
            circuit.append(
                "DETECTOR",
                [self.rel(d_now), self.rel(d_prev)],
                [ancilla[0], ancilla[1], self.current_round - self.coord_base],
            )
        return circuit

    def data_readout(self, flip: float = 0.0) -> Tuple[stim.Circuit, Dict[Coord, int]]:
        # Final data readout: measure each data qubit in its memory basis.
        # Returns (circuit, coord -> absolute measurement index).
        circuit = stim.Circuit()
        x_basis, z_basis = self.data_basis_partition()
        data_record: Dict[Coord, int] = {}
        for gate, group in (("MX", x_basis), ("M", z_basis)):
            if not group:
                continue
            idxs = [idx for _, idx in group]
            if flip > 0.0:
                circuit.append(gate, idxs, flip)
            else:
                circuit.append(gate, idxs)
            for coord, _ in group:
                data_record[coord] = self.record_counter
                self.record_counter += 1
        return circuit, data_record

    def define_observable(self, data_record: Dict[Coord, int]) -> stim.Circuit:
        # The logical observable
        circuit = stim.Circuit()
        logical = self.code.logical_x if self.basis == "X" else self.code.logical_z
        observable_targets = [self.rel(data_record[dc]) for dc in logical]
        circuit.append("OBSERVABLE_INCLUDE", observable_targets, 0)
        return circuit
