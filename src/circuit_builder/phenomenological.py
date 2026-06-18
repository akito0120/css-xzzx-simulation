import stim
from .noisy_measurement import NoisyMeasurementCircuitBuilder
from .shared import biased_pauli_rates

class PhenomenologicalCircuitBuilder(NoisyMeasurementCircuitBuilder):
    # Noise model: phenomenological

    def data_round_noise(self) -> stim.Circuit:
        # Per-round bulk data noise (the phenomenological data error).
        px, py, pz = biased_pauli_rates(self.p, self.eta)
        circuit = stim.Circuit()
        circuit.append("PAULI_CHANNEL_1",
                       list(self.code.data_qubits.values()), [px, py, pz])
        return circuit

    def build(self) -> stim.Circuit:
        circuit = stim.Circuit()
        self.ancilla_order = list(self.code.stabilizers.keys())

        # Initialize qubits in their memory basis
        circuit += self.init_qubit_coords()
        circuit += self.prep_data()

        # Round 0: perfect reference (projects into the code space)
        circuit += self.syndrome_meas()
        circuit.append("TICK")

        # Noisy rounds
        def round_body() -> stim.Circuit:
            c = stim.Circuit()
            c += self.data_round_noise()
            c += self.syndrome_meas(flip=self.p_meas)
            c += self.consecutive_round_detectors()
            return c
        circuit += self.repeat_rounds(self.rounds, round_body)

        # Final perfect data readout, time-boundary detectors and logical observable
        readout, data_record = self.data_readout()
        circuit += readout
        circuit += self.final_boundary_detectors(data_record)
        circuit += self.define_observable(data_record)

        return circuit
