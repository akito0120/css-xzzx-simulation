import stim
from qec_code import QecCode, Coord
from typing import Tuple, Dict, List

CGATE = {"X": "CX", "Y": "CY", "Z": "CZ"}

def biased_pauli_rates(p: float, eta: float) -> Tuple[float, float, float]:
    if eta == float("inf"):
        return 0.0, 0.0, p
    px = py = p / (2.0 * (1.0 + eta))
    pz = p * eta / (1.0 + eta)
    return px, py, pz

class CircuitBuilder:
    # Noise model: code capacity
    # Memory basis: X

    code: QecCode
    p: float
    eta: float

    circuit: stim.Circuit

    record_counter: int
    current_round: int
    ancilla_record: Dict[Tuple[Coord, int], int]
    ancilla_order: List[Coord]

    def __init__(self, code: QecCode, p: float, eta: float):
        self.record_counter = 0
        self.current_round = 0
        self.ancilla_record = {}

        self.code = code
        self.p = p
        self.eta = eta

    def syndrome_meas(self):
        # Measure ancillas and record the (coordinate, round) -> absolute index information
        for ancilla in self.ancilla_order:
            ancilla_idx = self.code.ancilla_qubits[ancilla]
            self.circuit.append("RX", [ancilla_idx])
        
        for ancilla in self.ancilla_order:
            ancilla_idx = self.code.ancilla_qubits[ancilla]
            for dcoord, pauli in self.code.stabilizers[ancilla].items():
                data_idx = self.code.data_qubits[dcoord]
                self.circuit.append(CGATE[pauli], [ancilla_idx, data_idx])
        
        for ancilla in self.ancilla_order:
            ancilla_idx = self.code.ancilla_qubits[ancilla]
            self.circuit.append("MX", [ancilla_idx])
            self.ancilla_record[(ancilla, self.current_round)] = self.record_counter
            self.record_counter += 1

    def rel(self, abs_idx):
        # Return the relative index of measurement result
        return stim.target_rec(abs_idx - self.record_counter)

    def deform_x_basis_data(self):
        # Apply H gate to prepare data qubits in X basis
        x_basis_data_list = [idx for coord, idx in self.code.data_qubits.items() if coord not in self.code.hset]
        if x_basis_data_list:
            self.circuit.append("H", x_basis_data_list)
    
    def build(self) -> stim.Circuit:
        self.circuit = stim.Circuit()
        self.ancilla_order = list(self.code.stabilizers.keys())

        px, py, pz = biased_pauli_rates(self.p, self.eta)
        data_list = list(self.code.data_qubits.values())

        # Initialize qubits
        for coord, idx in self.code.data_qubits.items():
            self.circuit.append("QUBIT_COORDS", [idx], [coord[0], coord[1]])
        for coord, idx in self.code.ancilla_qubits.items():
            self.circuit.append("QUBIT_COORDS", [idx], [coord[0], coord[1]])

        self.circuit.append("R", data_list)
        self.deform_x_basis_data()

        # Create the reference point of detector
        # Map the state to the code space
        self.syndrome_meas()
        self.circuit.append("TICK")

        # Noisy round
        self.current_round += 1
        self.circuit.append("PAULI_CHANNEL_1", data_list, [px, py, pz])

        # Syndrome extraction
        self.syndrome_meas()
        for ancilla in self.ancilla_order:
            d_now = self.ancilla_record[(ancilla, self.current_round)]
            d_prev = self.ancilla_record[(ancilla, self.current_round - 1)]
            self.circuit.append(
                "DETECTOR",
                [self.rel(d_now), self.rel(d_prev)],
                [ancilla[0], ancilla[1], self.current_round],
            )
        
        # Data readout
        self.deform_x_basis_data()
        self.circuit.append("M", data_list)

        data_record = {}
        for dcoord, idx in self.code.data_qubits.items():
            data_record[dcoord] = self.record_counter
            self.record_counter += 1

        # Final-round detector
        # TODO: uncomment for phenomenological and circuit-level noise model
        # for ancilla in self.ancilla_order:
        #     legs = self.code.stabilizers[ancilla]
        #     if(all(pauli == ("Z" if dc in self.code.hset else "X") for dc, pauli in legs.items())):
        #         targets = [self.rel(data_record[dc]) for dc in legs]
        #         targets.append(self.rel(self.ancilla_record[(ancilla, 1)]))
        #         self.circuit.append("DETECTOR", targets, [ancilla[0], ancilla[1], 2])

        # Define observable
        observable_targets = [self.rel(data_record[dc]) for dc in self.code.logical_x]
        self.circuit.append("OBSERVABLE_INCLUDE", observable_targets, 0)

        return self.circuit
