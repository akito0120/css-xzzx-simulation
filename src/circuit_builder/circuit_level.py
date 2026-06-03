from .shared import biased_pauli_rates, CGATE
from .phenomenological import PhenomenologicalCircuitBuilder

class CircuitLevelCircuitBuilder(PhenomenologicalCircuitBuilder):
    # Noise model: circuit-level
    # Memory basis: X

    def syndrome_meas(self, flip: float = 0.0):
        # current_round == 0 is the perfect reference round; everything later is noisy.
        noisy = self.current_round > 0
        px, py, pz = biased_pauli_rates(self.p, self.eta)

        ancilla_idxs = [self.code.ancilla_qubits[a] for a in self.ancilla_order]
        self.circuit.append("RX", ancilla_idxs)
        if noisy:
            # Reset preparation error
            self.circuit.append("PAULI_CHANNEL_1", ancilla_idxs, [px, py, pz])

        for ancilla in self.ancilla_order:
            ancilla_idx = self.code.ancilla_qubits[ancilla]
            for dcoord, pauli in self.code.stabilizers[ancilla].items():
                data_idx = self.code.data_qubits[dcoord]
                self.circuit.append(CGATE[pauli], [ancilla_idx, data_idx])
                if noisy:
                    # Two-qubit gate error: independent biased 1q channel on each qubit
                    self.circuit.append("PAULI_CHANNEL_1", [ancilla_idx, data_idx], [px, py, pz])

        for ancilla in self.ancilla_order:
            ancilla_idx = self.code.ancilla_qubits[ancilla]
            if flip > 0.0:
                self.circuit.append("MX", [ancilla_idx], flip)
            else:
                self.circuit.append("MX", [ancilla_idx])
            self.ancilla_record[(ancilla, self.current_round)] = self.record_counter
            self.record_counter += 1
