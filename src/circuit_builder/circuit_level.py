import stim
from .noisy_measurement import NoisyMeasurementCircuitBuilder
from .shared import biased_pauli_rates, biased_two_qubit_rates, build_cnot_schedule, CGATE

class CircuitLevelCircuitBuilder(NoisyMeasurementCircuitBuilder):
    # Noise model: circuit-level (biased-SD6)
    # Memory basis: X

    def syndrome_meas(self):
        # One full noisy round of syndrome extraction
        px, py, pz = biased_pauli_rates(self.p, self.eta)
        pc2 = biased_two_qubit_rates(self.p, self.eta)

        ancilla_idxs = [self.code.ancilla_qubits[a] for a in self.ancilla_order]
        data_list = list(self.code.data_qubits.values())
        all_qubits = data_list + ancilla_idxs

        # Phase 1: Reset
        self.circuit.append("RX", ancilla_idxs)
        # Ancilla reset preparation error
        self.circuit.append("Z_ERROR", ancilla_idxs, self.p)
        # Data qubit idle error in reset window
        self.circuit.append("PAULI_CHANNEL_1", data_list, [px, py, pz])

        # Phase 2: 4-step CX/CZ
        schedule = build_cnot_schedule(self.code)
        for step in range(4):
            busy_qubits: set = set()
            for ancilla in self.ancilla_order:
                dcoord = schedule[ancilla][step]
                if dcoord is None:
                    continue
                ancilla_idx = self.code.ancilla_qubits[ancilla]
                data_idx = self.code.data_qubits[dcoord]
                pauli = self.code.stabilizers[ancilla][dcoord]
                self.circuit.append(CGATE[pauli], [ancilla_idx, data_idx])
                # Two-qubit gate error (biased SD6)
                self.circuit.append("PAULI_CHANNEL_2", [ancilla_idx, data_idx], pc2)
                busy_qubits.add(ancilla_idx)
                busy_qubits.add(data_idx)
            # Idle error: qubits resting this step (boundary data + idle ancillas)
            idle_qubits = [q for q in all_qubits if q not in busy_qubits]
            if idle_qubits:
                self.circuit.append("PAULI_CHANNEL_1", idle_qubits, [px, py, pz])
            self.circuit.append("TICK")

        # Phase 3: Measurement
        # Data qubit idle error in measurement window
        self.circuit.append("PAULI_CHANNEL_1", data_list, [px, py, pz])
        for ancilla in self.ancilla_order:
            ancilla_idx = self.code.ancilla_qubits[ancilla]
            if self.p_meas > 0.0:
                self.circuit.append("MX", [ancilla_idx], self.p_meas)
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
        # The following H deformation turns this into the per-basis prep error (Z on |+>, X on |0>)
        self.circuit.append("R", data_list)
        self.circuit.append("X_ERROR", data_list, self.p)
        self.deform_x_basis_data()

        # Round 0 (noisy)
        self.syndrome_meas()
        self.initial_boundary_detectors()
        self.circuit.append("TICK")

        # Remaining noisy rounds
        # Unlike phenomenological, data noise comes from gate and idle errors in syndrome_meas.
        for _ in range(self.rounds - 1):
            self.current_round += 1
            self.syndrome_meas()
            self.consecutive_round_detectors()
            self.circuit.append("TICK")

        # Final data readout (with readout error), time-boundary detectors and logical observable
        data_record = self.data_readout(flip=self.p_meas)
        self.final_boundary_detectors(data_record)
        self.define_observable(data_record)

        return self.circuit
