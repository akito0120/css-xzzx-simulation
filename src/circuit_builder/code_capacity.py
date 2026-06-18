import stim
from .base import BaseCircuitBuilder
from .shared import biased_pauli_rates

class CodeCapacityCircuitBuilder(BaseCircuitBuilder):
    # Noise model: code capacity

    def build(self) -> stim.Circuit:
        circuit = stim.Circuit()
        self.ancilla_order = list(self.code.stabilizers.keys())

        px, py, pz = biased_pauli_rates(self.p, self.eta)
        data_list = list(self.code.data_qubits.values())

        # Initialize qubits in their memory basis
        circuit += self.init_qubit_coords()
        circuit += self.prep_data()

        # Create the reference point of detector
        # Map the state to the code space
        circuit += self.syndrome_meas()
        circuit.append("TICK")

        # Noisy round
        self.current_round += 1
        circuit.append("PAULI_CHANNEL_1", data_list, [px, py, pz])

        # Syndrome extraction
        circuit += self.syndrome_meas()
        circuit += self.consecutive_round_detectors()

        # Data readout + logical observable (no final-round detectors)
        readout, data_record = self.data_readout()
        circuit += readout
        circuit += self.define_observable(data_record)

        return circuit
