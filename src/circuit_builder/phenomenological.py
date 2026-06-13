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

        # Initialize qubits in their memory basis
        self.init_qubit_coords()
        self.prep_data()

        # Round 0: perfect reference (projects into the code space)
        self.syndrome_meas()
        self.circuit.append("TICK")

        # Noisy rounds
        def round_body():
            self.data_round_noise()
            self.syndrome_meas(flip=self.p_meas)
            self.consecutive_round_detectors()
        self.repeat_rounds(self.rounds, round_body)

        # Final perfect data readout, time-boundary detectors and logical observable
        data_record = self.data_readout()
        self.final_boundary_detectors(data_record)
        self.define_observable(data_record)

        return self.circuit
