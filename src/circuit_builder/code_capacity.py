import stim
from .base import BaseCircuitBuilder
from .shared import biased_pauli_rates

class CodeCapacityCircuitBuilder(BaseCircuitBuilder):
    # Noise model: code capacity
    # Memory basis: X

    def build(self) -> stim.Circuit:
        self.circuit = stim.Circuit()
        self.ancilla_order = list(self.code.stabilizers.keys())

        px, py, pz = biased_pauli_rates(self.p, self.eta)
        data_list = list(self.code.data_qubits.values())

        # Initialize qubits in their memory basis
        self.init_qubit_coords()
        self.prep_data()

        # Create the reference point of detector
        # Map the state to the code space
        self.syndrome_meas()
        self.circuit.append("TICK")

        # Noisy round
        self.current_round += 1
        self.circuit.append("PAULI_CHANNEL_1", data_list, [px, py, pz])

        # Syndrome extraction
        self.syndrome_meas()
        self.consecutive_round_detectors()

        # Data readout + logical observable (no final-round detectors)
        data_record = self.data_readout()
        self.define_observable(data_record)

        return self.circuit
