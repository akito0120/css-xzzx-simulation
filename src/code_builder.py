from qec_code import QecCode, Coord
from typing import Dict, List, Tuple

def build_rotated_surface_code(distance: int) -> QecCode:
    # Data qubits sit at (2c+1, 2r+1) for r,c in 0..d-1.
    # Plaquette centres sit at even-even coordinates (2c, 2r).

    # Building data qubits
    data: Dict[Coord, int] = {}
    for r in range(distance):
        for c in range(distance):
            data[(2 * c + 1, 2 * r + 1)] = len(data)

    # Building ancilla qubits and stabilizers
    ancillas: Dict[Coord, int] = {}
    stabilizers: Dict[Coord, Dict[Coord, str]] = {}

    next_ancilla = len(data)
    for r in range(0, distance + 1):
        for c in range(0, distance + 1):
            x, y = 2 * c, 2 * r
            # data neighbours = the 4 diagonal data sites
            neighbours = [
                (x - 1, y - 1), (x + 1, y - 1),
                (x - 1, y + 1), (x + 1, y + 1),
            ]
            present = [q for q in neighbours if q in data]
            if len(present) <= 1:
                continue
            # Checkerboard type by (r+c) parity
            stabilizer_type = 'Z' if (r + c) % 2 == 0 else 'X'
            # Boundary handling: keep weight-2 checks only on the right boundaries.
            if len(present) == 2:
                # determine if this is a top/bottom (horizontal) or
                # left/right (vertical) boundary plaquette
                on_vertical_boundary = (c == 0 or c == distance)   # left/right edge
                on_horizontal_boundary = (r == 0 or r == distance)  # top/bottom edge
                # Standard rotated code: Z checks on top/bottom, X on left/right
                # (consistent with the (r+c) colouring chosen above).
                if on_horizontal_boundary and stabilizer_type != 'Z':
                    continue
                if on_vertical_boundary and stabilizer_type != 'X':
                    continue
                # reject corners with only the wrong-type 2-body check
                if on_horizontal_boundary and on_vertical_boundary:
                    continue
            ancillas[(x, y)] = next_ancilla
            next_ancilla += 1
            stabilizers[(x, y)] = {}
            
            for q in present:
                stabilizers[(x, y)][q] = stabilizer_type
    
    logical_x = {(2 * c + 1, 1): 'X' for c in range(distance)}
    logical_z = {(1, 2 * r + 1): 'Z' for r in range(distance)}

    return QecCode(
        distance=distance,
        name=f"rotated_surface_d{distance}",
        data_qubits=data,
        ancilla_qubits=ancillas,
        stabilizers=stabilizers,
        logical_x=logical_x,
        logical_z=logical_z
    )

H_CONJ = {"X": "Z", "Z": "X", "Y": "Y"}

def deform(pauli_dict: Dict[Coord, str], hset: set) -> Dict[Coord, str]:
    return {q: (H_CONJ[p] if q in hset else p) for q, p in pauli_dict.items()}

def build_xzzx_code(distance: int) -> QecCode:
    # Build from rotated surface code
    rotated_surface_code = build_rotated_surface_code(distance)

    hset = set()
    for (cx, cy) in rotated_surface_code.data_qubits:
        gc, gr = (cx - 1) // 2, (cy - 1) // 2
        if (gc + gr) % 2 == 1:
            hset.add((cx, cy))
    
    stabilizers = {anc: deform(legs, hset)
                   for anc, legs in rotated_surface_code.stabilizers.items()}

    logical_x = deform(rotated_surface_code.logical_x, hset)
    logical_z = deform(rotated_surface_code.logical_z, hset)

    return QecCode(
        distance=distance,
        name=f"xzzx_d{distance}",
        data_qubits=rotated_surface_code.data_qubits,
        ancilla_qubits=rotated_surface_code.ancilla_qubits,
        stabilizers=stabilizers,
        logical_x=logical_x,
        logical_z=logical_z,
        hset=hset
    )
