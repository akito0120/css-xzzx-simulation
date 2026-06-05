import stim
from .noisy_measurement import NoisyMeasurementCircuitBuilder
from .shared import biased_pauli_rates, biased_two_qubit_rates, CGATE

class CircuitLevelCircuitBuilder(NoisyMeasurementCircuitBuilder):
    # Noise model: circuit-level
    # Memory basis: X

    def syndrome_meas(self, flip: float = 0.0):
        # current_round == 0 is the perfect reference round; everything later is noisy.
        noisy = self.current_round > 0
        px, py, pz = biased_pauli_rates(self.p, self.eta)
        pc2 = biased_two_qubit_rates(self.p, self.eta)

        ancilla_idxs = [self.code.ancilla_qubits[a] for a in self.ancilla_order]
        data_list = list(self.code.data_qubits.values())
        self.circuit.append("RX", ancilla_idxs)
        if noisy:
            # Reset preparation error
            self.circuit.append("PAULI_CHANNEL_1", ancilla_idxs, [px, py, pz])
            # Idle error: data qubits wait while the ancillas are reset
            self.circuit.append("PAULI_CHANNEL_1", data_list, [px, py, pz])

        for ancilla in self.ancilla_order:
            ancilla_idx = self.code.ancilla_qubits[ancilla]
            for dcoord, pauli in self.code.stabilizers[ancilla].items():
                data_idx = self.code.data_qubits[dcoord]
                self.circuit.append(CGATE[pauli], [ancilla_idx, data_idx])
                if noisy:
                    # Two-qubit gate error: biased correlated 2-qubit Pauli channel
                    self.circuit.append("PAULI_CHANNEL_2", [ancilla_idx, data_idx], pc2)

        if noisy:
            # Idle error: data qubits wait while the ancillas are measured
            self.circuit.append("PAULI_CHANNEL_1", data_list, [px, py, pz])

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

        data_list = list(self.code.data_qubits.values())

        self.init_qubit_coords()

        # Prepare |+> / |0>
        self.circuit.append("R", data_list)
        self.deform_x_basis_data()

        # Round 0: perfect reference (projects into the code space)
        self.syndrome_meas()
        self.circuit.append("TICK")

        # Noisy rounds. Unlike phenomenological, there is no bulk per-round data channel.
        # The data noise comes from the gate and idle errors inside syndrome_meas.
        for _ in range(self.rounds):
            self.current_round += 1
            self.syndrome_meas(flip=self.p_meas)
            self.consecutive_round_detectors()
            self.circuit.append("TICK")

        # Final perfect data readout, time-boundary detectors and logical observable
        data_record = self.data_readout()
        self.final_boundary_detectors(data_record)
        self.define_observable(data_record)

        return self.circuit
