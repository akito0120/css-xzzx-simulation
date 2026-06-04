import stim
from .noisy_measurement import NoisyMeasurementCircuitBuilder
from .shared import biased_pauli_rates

class PhenomenologicalCircuitBuilder(NoisyMeasurementCircuitBuilder):
    # Noise model: phenomenological
    # Memory basis: X

    def data_round_noise(self):
        # Per-round bulk data noise (the phenomenological data error).
        px, py, pz = biased_pauli_rates(self.p, self.eta)
        self.circuit.append("PAULI_CHANNEL_1",
                            list(self.code.data_qubits.values()), [px, py, pz])

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

        # Noisy rounds
        for _ in range(self.rounds):
            self.current_round += 1
            self.data_round_noise()
            self.syndrome_meas(flip=self.p_meas)

            for ancilla in self.ancilla_order:
                d_now = self.ancilla_record[(ancilla, self.current_round)]
                d_prev = self.ancilla_record[(ancilla, self.current_round - 1)]
                self.circuit.append(
                    "DETECTOR",
                    [self.rel(d_now), self.rel(d_prev)],
                    [ancilla[0], ancilla[1], self.current_round],
                )
            self.circuit.append("TICK")

        # Final perfect data readout, time-boundary detectors and logical observable
        data_record = self.data_readout()
        self.final_boundary_detectors(data_record)
        self.define_observable(data_record)

        return self.circuit
