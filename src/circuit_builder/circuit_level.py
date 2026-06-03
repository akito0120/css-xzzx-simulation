import stim
from typing import Dict, List, Optional
from code_builder import QecCode, Coord
from .shared import BaseCircuitBuilder, biased_pauli_rates, CGATE

class CircuitLevelCircuitBuilder(BaseCircuitBuilder):
    # Noise model: circuit-level
    # Memory basis: X

    # Self-contained on purpose: it owns build / __init__ / time-boundary detectors
    # and only reuses the model-independent plumbing from BaseCircuitBuilder
    # (rel, deform_x_basis_data, init_qubit_coords, data_readout_and_observable).

    rounds: int
    p_meas: float

    def __init__(self, code: QecCode, p: float, eta: float,
                 rounds: Optional[int] = None, p_meas: Optional[float] = None):
        super().__init__(code, p, eta)
        self.rounds = code.distance if rounds is None else rounds
        self.p_meas = p if p_meas is None else p_meas

    def syndrome_meas(self, flip: float = 0.0):
        # current_round == 0 is the perfect reference round; everything later is noisy.
        noisy = self.current_round > 0
        px, py, pz = biased_pauli_rates(self.p, self.eta)

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
                    # Two-qubit gate error: independent biased 1q channel on each qubit
                    self.circuit.append("PAULI_CHANNEL_1", [ancilla_idx, data_idx], [px, py, pz])

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

    def final_boundary_detectors(self, data_record: Dict[Coord, int]):
        # Per-qubit preparation/measurement basis for X-memory
        # |+> (basis X) off hset, |0> (basis Z) on hset
        prep_basis: Dict[Coord, str] = {
            q: ("Z" if q in self.code.hset else "X") for q in self.code.data_qubits
        }

        # "eligible" stabilizers = checks deterministic under that preparation
        # Every leg matches the prep basis
        # Checks reconstructable from the final data readout
        eligible: List[Coord] = [
            anc for anc, legs in self.code.stabilizers.items()
            if all(pauli == prep_basis[q] for q, pauli in legs.items())
        ]

        # Close the time boundary: reconstructed check vs. last round
        last = self.current_round
        for ancilla in eligible:
            targets = [self.rel(data_record[dc])
                       for dc in self.code.stabilizers[ancilla]]
            targets.append(self.rel(self.ancilla_record[(ancilla, last)]))
            self.circuit.append("DETECTOR", targets,
                                [ancilla[0], ancilla[1], last + 1])

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
        self.data_readout_and_observable()

        return self.circuit
