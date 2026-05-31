import stim
from qec_code import QecCode, Coord
from typing import Tuple, Dict

CGATE = {"X": "CX", "Y": "CY", "Z": "CZ"}

def biased_pauli_rates(p: float, eta: float) -> Tuple[float, float, float]:
    if eta == float("inf"):
        return 0.0, 0.0, p
    px = py = p / (2.0 * (1.0 + eta))
    pz = p * eta / (1.0 + eta)
    return px, py, pz

def build_circuit(code: QecCode, p: float, eta: float) -> stim.Circuit:
    # Noise model: code capacity
    # Memory: X
    
    px, py, pz = biased_pauli_rates(p, eta)
    circuit = stim.Circuit()

    # Initialize qubits
    for coord, idx in code.data_qubits.items():
        circuit.append("QUBIT_COORDS", [idx], [coord[0], coord[1]])
    for coord, idx in code.ancilla_qubits.items():
        circuit.append("QUBIT_COORDS", [idx], [coord[0], coord[1]])

    data_list = list(code.data_qubits.values())
    x_basis_data_list = [idx for coord, idx in code.data_qubits.items() if coord not in code.hset]
    circuit.append("R", data_list)
    if x_basis_data_list:
        circuit.append("H", x_basis_data_list)

    circuit.append("TICK")

    record_counter = 0
    current_round = 0
    ancilla_record: Dict[Tuple[Coord, int], int] = {} # (coordinate, round) -> absolute index
    ancilla_order = list(code.stabilizers.keys()) # To fix the order of measurement
    
    def syndrome_meas():
        # Measure ancillas and record the (coordinate, round) -> absolute index information
        nonlocal record_counter, current_round
        for ancilla in ancilla_order:
            ancilla_idx = code.ancilla_qubits[ancilla]
            circuit.append("RX", [ancilla_idx])
        
        for ancilla in ancilla_order:
            ancilla_idx = code.ancilla_qubits[ancilla]
            for dcoord, pauli in code.stabilizers[ancilla].items():
                data_idx = code.data_qubits[dcoord]
                circuit.append(CGATE[pauli], [ancilla_idx, data_idx])
        
        for ancilla in ancilla_order:
            ancilla_idx = code.ancilla_qubits[ancilla]
            circuit.append("MX", [ancilla_idx])
            ancilla_record[(ancilla, current_round)] = record_counter
            record_counter += 1

    # Create the reference point of detector
    # Map the state to the code space
    syndrome_meas()

    circuit.append("TICK")

    def rel(abs_idx):
        # Return the relative index of measurement result
        nonlocal record_counter
        return stim.target_rec(abs_idx - record_counter)

    # Noisy round (1 round)
    current_round += 1
    if p > 0:
        circuit.append("PAULI_CHANNEL_1", data_list, [px, py, pz])

    # Syndrome extraction after noisy round
    syndrome_meas()
    for ancilla in ancilla_order:
        d_now = ancilla_record[(ancilla, current_round)]
        d_prev = ancilla_record[(ancilla, current_round - 1)]
        circuit.append(
            "DETECTOR",
            [rel(d_now), rel(d_prev)],
            [ancilla[0], ancilla[1], current_round],
        )
    circuit.append("TICK")

    # Data readout
    if x_basis_data_list:
        circuit.append("H", x_basis_data_list)
    circuit.append("M", data_list)

    data_record = {}
    for dcoord, idx in code.data_qubits.items():
        data_record[dcoord] = record_counter
        record_counter += 1

    # Final-round detector
    for ancilla in ancilla_order:
        legs = code.stabilizers[ancilla]
        if(all(pauli == ("Z" if dc in code.hset else "X") for dc, pauli in legs.items())):
            targets = [rel(data_record[dc]) for dc in legs]
            targets.append(rel(ancilla_record[(ancilla, 1)]))
            circuit.append("DETECTOR", targets, [ancilla[0], ancilla[1], 2])
    
    # Define observable
    observable_targets = [rel(data_record[dc]) for dc in code.logical_x]
    circuit.append("OBSERVABLE_INCLUDE", observable_targets, 0)

    return circuit
