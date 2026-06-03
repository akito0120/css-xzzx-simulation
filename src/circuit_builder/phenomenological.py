import stim
from typing import Dict, List, Optional
from code_builder import QecCode, Coord
from .shared import BaseCircuitBuilder, biased_pauli_rates

class PhenomenologicalCircuitBuilder(BaseCircuitBuilder):
    # Noise model: phenomenological (per-round data noise + measurement flips)
    # Memory basis: X

    rounds: int
    p_meas: float

    def __init__(self, code: QecCode, p: float, eta: float,
                 rounds: Optional[int] = None, p_meas: Optional[float] = None):
        super().__init__(code, p, eta)
        self.rounds = code.distance if rounds is None else rounds
        self.p_meas = p if p_meas is None else p_meas

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

        px, py, pz = biased_pauli_rates(self.p, self.eta)
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
            self.circuit.append("PAULI_CHANNEL_1", data_list, [px, py, pz])
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
