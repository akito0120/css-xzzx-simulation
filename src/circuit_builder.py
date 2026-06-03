import stim
from typing import Tuple, Dict, List, Optional
from code_builder import QecCode, Coord

CGATE = {"X": "CX", "Y": "CY", "Z": "CZ"}

def biased_pauli_rates(p: float, eta: float) -> Tuple[float, float, float]:
    if eta == float("inf"):
        return 0.0, 0.0, p
    px = py = p / (2.0 * (1.0 + eta))
    pz = p * eta / (1.0 + eta)
    return px, py, pz

class CodeCapacityCircuitBuilder:
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

        # Define observable
        observable_targets = [self.rel(data_record[dc]) for dc in self.code.logical_x]
        self.circuit.append("OBSERVABLE_INCLUDE", observable_targets, 0)

        return self.circuit

class PhenomenologicalCircuitBuilder:
    # Noise model: phenomenological
    # Memory basis: X

    code: QecCode
    p: float
    eta: float
    rounds: int
    p_meas: float
 
    circuit: stim.Circuit

    record_counter: int
    current_round: int
    ancilla_record: Dict[Tuple[Coord, int], int]
    ancilla_order: List[Coord]
 
    def __init__(self, code: QecCode, p: float, eta: float, 
                 rounds: Optional[int] = None, p_meas: Optional[float] = None):
        self.code = code
        self.p = p
        self.eta = eta
        self.rounds = code.distance if rounds is None else rounds
        self.p_meas = p if p_meas is None else p_meas
 
        self.record_counter = 0
        self.current_round = 0
        self.ancilla_record = {}
 
    def rel(self, abs_idx: int) -> stim.GateTarget:
        return stim.target_rec(abs_idx - self.record_counter)
 
    def deform_x_basis_data(self):
        x_basis_data_list = [idx for coord, idx in self.code.data_qubits.items() if coord not in self.code.hset]
        if x_basis_data_list:
            self.circuit.append("H", x_basis_data_list)
 
    def syndrome_meas(self, noisy: bool):
        # 1 round of syndrome extraction
        # noisy=True injects a readout flip
        flip = self.p_meas if noisy else 0.0
 
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
 
    def build(self) -> stim.Circuit:
        self.circuit = stim.Circuit()
        self.ancilla_order = list(self.code.stabilizers.keys())
 
        px, py, pz = biased_pauli_rates(self.p, self.eta)
        data_list = list(self.code.data_qubits.values())
 
        for coord, idx in self.code.data_qubits.items():
            self.circuit.append("QUBIT_COORDS", [idx], [coord[0], coord[1]])
        for coord, idx in self.code.ancilla_qubits.items():
            self.circuit.append("QUBIT_COORDS", [idx], [coord[0], coord[1]])
 
        # Prepare |+> / |0> 
        self.circuit.append("R", data_list)
        self.deform_x_basis_data()
 
        # Round 0: perfect reference (projects into the code space) 
        self.syndrome_meas(noisy=False)
        self.circuit.append("TICK")
 
        # Noisy rounds
        for _ in range(self.rounds):
            self.current_round += 1
            self.circuit.append("PAULI_CHANNEL_1", data_list, [px, py, pz])
            self.syndrome_meas(noisy=True)
 
            for ancilla in self.ancilla_order:
                d_now = self.ancilla_record[(ancilla, self.current_round)]
                d_prev = self.ancilla_record[(ancilla, self.current_round - 1)]
                self.circuit.append(
                    "DETECTOR",
                    [self.rel(d_now), self.rel(d_prev)],
                    [ancilla[0], ancilla[1], self.current_round],
                )
            self.circuit.append("TICK")
 
        # Final perfect data readout
        self.deform_x_basis_data()
        self.circuit.append("M", data_list)
 
        data_record: Dict[Coord, int] = {}
        for dcoord, idx in self.code.data_qubits.items():
            data_record[dcoord] = self.record_counter
            self.record_counter += 1
 
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
 
        # Logical X observable
        observable_targets = [self.rel(data_record[dc]) for dc in self.code.logical_x]
        self.circuit.append("OBSERVABLE_INCLUDE", observable_targets, 0)
 
        return self.circuit
