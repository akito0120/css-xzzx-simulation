from dataclasses import dataclass, field
from typing import Dict, Tuple

Coord = Tuple[int, int]

@dataclass(frozen=True)
class QecCode:
    distance: int # Code distance
    name: str # Code name
    data_qubits: Dict[Coord, int] # Data qubit coordinate on lattice -> qubit index
    ancilla_qubits: Dict[Coord, int] # Ancilla qubit coordinate on lattice -> qubit index
    stabilizers: Dict[Coord, Dict[Coord, str]] # Ancilla coordinate -> (data qubit coordinate -> X/Z)
    logical_x: Dict[Coord, str] # Logical X operator
    logical_z: Dict[Coord, str] # Logical Z operator
    hset: set = field(default_factory=set)
