import stim
from typing import Dict
from code_builder import Coord
from .shared import BaseCircuitBuilder, biased_pauli_rates

class CodeCapacityCircuitBuilder(BaseCircuitBuilder):
    # Noise model: code capacity (single noisy data round, perfect measurements)
    # Memory basis: X

    def final_boundary_detectors(self, data_record: Dict[Coord, int]):
        # No final-round detectors: the readout is noiseless, so the round-1
        # syndrome already captures every data error. Opt out of the base impl.
        pass

    def build(self) -> stim.Circuit:
        self.circuit = stim.Circuit()
        self.ancilla_order = list(self.code.stabilizers.keys())

        px, py, pz = biased_pauli_rates(self.p, self.eta)
        data_list = list(self.code.data_qubits.values())

        # Initialize qubits
        self.init_qubit_coords()
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

        # Data readout + logical observable (no final-round detectors; see above)
        self.data_readout_and_observable()

        return self.circuit
