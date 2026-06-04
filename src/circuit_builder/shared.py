import stim
from typing import Tuple, Dict, List
from code_builder import QecCode, Coord

CGATE = {"X": "CX", "Y": "CY", "Z": "CZ"}

def biased_pauli_rates(p: float, eta: float) -> Tuple[float, float, float]:
    if eta == float("inf"):
        return 0.0, 0.0, p
    px = py = p / (2.0 * (1.0 + eta))
    pz = p * eta / (1.0 + eta)
    return px, py, pz


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
    ancilla_record: Dict[Tuple[Coord, int], int]
    ancilla_order: List[Coord]

    def __init__(self, code: QecCode, p: float, eta: float):
        self.code = code
        self.p = p
        self.eta = eta

        self.record_counter = 0
        self.current_round = 0
        self.ancilla_record = {}

    def rel(self, abs_idx: int) -> stim.GateTarget:
        # Return the relative index of a measurement result
        return stim.target_rec(abs_idx - self.record_counter)

    def deform_x_basis_data(self):
        # Apply H to data qubits prepared/measured in the X basis (those off hset)
        x_basis_data_list = [idx for coord, idx in self.code.data_qubits.items()
                             if coord not in self.code.hset]
        if x_basis_data_list:
            self.circuit.append("H", x_basis_data_list)

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

    def final_boundary_detectors(self, data_record: Dict[Coord, int]):
        # Close the time boundary by reconstructing checks from the final data readout.
        # Used wherever the bulk rounds have noisy measurement (phenomenological, circuit-level).
        # Overridden as a no-op in code capacity, where the final readout is noiseless.

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

    def data_readout_and_observable(self) -> Dict[Coord, int]:
        # Final perfect data readout in the X-memory basis, then the logical-X observable.
        # Returns coord -> absolute measurement index.
        self.deform_x_basis_data()
        self.circuit.append("M", list(self.code.data_qubits.values()))

        data_record: Dict[Coord, int] = {}
        for dcoord, idx in self.code.data_qubits.items():
            data_record[dcoord] = self.record_counter
            self.record_counter += 1

        self.final_boundary_detectors(data_record)

        observable_targets = [self.rel(data_record[dc]) for dc in self.code.logical_x]
        self.circuit.append("OBSERVABLE_INCLUDE", observable_targets, 0)
        return data_record
