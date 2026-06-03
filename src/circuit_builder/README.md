# Circuit Builders — Explained

How the CSS rotated surface code / XZZX code are assembled into Stim circuits under three
noise models (code capacity / phenomenological / circuit-level), and the design rationale behind each.

Files covered:
- [shared.py](css_xzzx_simulation/src/circuit_builder/shared.py) — common base class `BaseCircuitBuilder`
- [code_capacity.py](css_xzzx_simulation/src/circuit_builder/code_capacity.py) — code capacity model
- [phenomenological.py](css_xzzx_simulation/src/circuit_builder/phenomenological.py) — phenomenological model
- [circuit_level.py](css_xzzx_simulation/src/circuit_builder/circuit_level.py) — circuit-level model

The codes themselves are defined in [code_builder.py](css_xzzx_simulation/src/code_builder.py).

## 0. The Big Picture

All three builders share the same skeleton — an **X-memory experiment** (hold a logical state in the X
basis and measure how well it survives) — and differ **only in how noise is injected**. The inheritance is:

```
BaseCircuitBuilder                 … shared parts: qubit init, syndrome extraction, readout
├── CodeCapacityCircuitBuilder     … data noise once, perfect measurements
└── PhenomenologicalCircuitBuilder … per-round data noise + measurement flips
    └── CircuitLevelCircuitBuilder  … the above plus reset/gate/idle noise
```

`CircuitLevel` inherits from `Phenomenological` for a reason, not by accident:
**circuit-level is a superset of phenomenological** (data + measurement noise, plus gate-induced noise).
Details below.

The "why" shared by all three models:

| Design choice | Reason |
|---|---|
| **X-memory basis** | The goal is to show XZZX's advantage under Z-biased noise. We hold a logical X and read out in the X basis at the end |
| **Round 0 is perfect (noiseless)** | To make the detector "reference point" deterministic. It projects the syndrome into the code space |
| **Detector = comparison of consecutive rounds' syndromes** | If nothing happened, the syndrome is unchanged. A change = an error occurred in between — a difference detector |
| **Single ancilla + controlled-Pauli** | Prepare one ancilla in \|+⟩, apply controlled Paulis to the data, and measure it in X — this measures any Pauli stabilizer |


## 1. Shared Building Blocks (BaseCircuitBuilder)

The common methods in [shared.py](css_xzzx_simulation/src/circuit_builder/shared.py).

### 1.1 Measurement-record relative indexing — `rel`

Stim's detectors/observables refer to "the N-th most recent measurement" via `target_rec(-k)`,
i.e. a **relative** reference. The builder, on the other hand, finds it easier to track
(ancilla coord, round) → **absolute** measurement number in `ancilla_record`.
`rel` bridges that gap:

```python
def rel(self, abs_idx: int) -> stim.GateTarget:
    return stim.target_rec(abs_idx - self.record_counter)
```

`record_counter` is "the number of measurements emitted so far." `abs_idx - record_counter` is always
negative and yields the correct relative offset at that moment. The key is to call `rel`
**at the exact point where the detector is appended**.

### 1.2 Basis deformation (Hadamard) — `deform_x_basis_data`

```python
x_basis_data_list = [idx for coord, idx in self.code.data_qubits.items()
                     if coord not in self.code.hset]
self.circuit.append("H", x_basis_data_list)
```

- **CSS code**: `hset` is empty. H is applied to every data qubit, preparing them all in \|+⟩ (pure X memory).
- **XZZX code**: the qubits in `hset` (a checkerboard subset of data qubits, see `build_xzzx_code` in
  [code_builder.py](css_xzzx_simulation/src/code_builder.py)) are **not** given an H → they stay in
  \|0⟩ (the Z basis).

Why: the XZZX code is the **Hadamard-deformed** version of the CSS code, obtained by applying H to the
data qubits in `hset`. The stabilizers and logical operators undergo the same deformation (the `deform`
function). On the circuit side, we choose each data qubit's preparation basis (\|+⟩ or \|0⟩) to match
that deformation, so that **the logical X in the deformed frame is deterministic**.
The same function is also called **right before the final readout** (H is self-inverse, so the
preparation basis and the measurement basis coincide).

### 1.3 Syndrome measurement — `syndrome_meas(flip=0.0)`

One round of syndrome extraction, in three stages:

```python
# (1) reset every ancilla to |+>
for ancilla in self.ancilla_order:
    self.circuit.append("RX", [self.code.ancilla_qubits[ancilla]])

# (2) controlled Pauli on each leg of the stabilizer
for ancilla in ...:
    for dcoord, pauli in self.code.stabilizers[ancilla].items():
        self.circuit.append(CGATE[pauli], [ancilla_idx, data_idx])  # CX/CY/CZ

# (3) measure the ancilla in the X basis; inject a flip if flip > 0
self.circuit.append("MX", [ancilla_idx], flip)   # or MX without noise
self.ancilla_record[(ancilla, current_round)] = record_counter
```

Why this measures the stabilizer:
put the ancilla in \|+⟩ (an X eigenstate) and apply controlled Paulis targeting the data; the ancilla's
phase flips according to the eigenvalue of the stabilizer operator $S$. Measuring the ancilla in the X
basis (`MX`) then yields the value of $S$. `CGATE = {"X":"CX","Y":"CY","Z":"CZ"}` selects the controlled
gate matching each leg's Pauli type ([shared.py](css_xzzx_simulation/src/circuit_builder/shared.py)).

The `flip` argument is the **injection point for measurement noise**. Switching it on/off is all it takes
to toggle between "perfect" and "noisy" measurements, which is why all three models can build on the same
`syndrome_meas`.

### 1.4 Final readout and logical observable — `data_readout_and_observable`

```python
self.deform_x_basis_data()                 # return to the same basis as preparation
self.circuit.append("M", data_list)        # measure all data in Z (effectively the X-memory readout)
# record coord -> absolute measurement number
self.final_boundary_detectors(data_record) # hook (see below)
observable_targets = [self.rel(data_record[dc]) for dc in self.code.logical_x]
self.circuit.append("OBSERVABLE_INCLUDE", observable_targets, 0)
```

- The logical-X observable is defined by the parity of the final measurements of the data qubits listed
  in `code.logical_x`.
- `final_boundary_detectors` is a **hook method**. In the base class it does nothing (no-op).
  It is unnecessary for code capacity (perfect measurements), but for models with noisy measurement it
  adds **time-boundary detectors** that cross-check "the stabilizer reconstructed from the final data
  readout" against "the last round's syndrome" (implemented via an override). Since this single point is
  the only model-specific difference, it is factored out as a hook.


## 2. Code Capacity Model

[code_capacity.py](css_xzzx_simulation/src/circuit_builder/code_capacity.py)

### Structure

```
QUBIT_COORDS …                     declare qubit coordinates
R(data); H(deform)                 prepare |+> / |0>
syndrome_meas()                    round 0: perfect reference (project into the code space)
TICK
PAULI_CHANNEL_1(data, [px,py,pz])  ★the only noise: one biased Pauli channel on data
syndrome_meas()                    round 1: perfect measurement
DETECTOR(round1 vs round0) × all ancillas
data_readout_and_observable()      perfect readout + logical X
```

### Why it looks like this

**Code capacity is the simplest model: errors land only on data qubits, while measurements and gates are
perfect.** Therefore:

- Noise is a **single** `PAULI_CHANNEL_1`. `biased_pauli_rates(p, eta)` provides
  $p_x=p_y=\dfrac{p}{2(1+\eta)},\ p_z=\dfrac{p\eta}{1+\eta}$ (summing to $p$, biased toward Z by a factor $\eta$).
- Measurements are perfect (`flip=0`).
- There is effectively a single round: "reference → noise → measurement."
- **No time-boundary detectors.** Because measurement is error-free, round 1's syndrome already
  represents the final state completely and is consistent with the final readout (the detectors would only
  be redundant). This decision is also noted in the TODO comment of the legacy implementation,
  [circuit_builder_legacy.py](css_xzzx_simulation/src/circuit_builder/circuit_builder_legacy.py).

The detectors compare round 1 against round 0 because (as in §0) we look at the difference of a syndrome
that should be unchanged in the absence of noise. Since round 0 establishes a perfect baseline, the
detectors that fire in round 1 are exactly the traces of the errors introduced by `PAULI_CHANNEL_1`.


## 3. Phenomenological Model

[phenomenological.py](css_xzzx_simulation/src/circuit_builder/phenomenological.py)

### Structure

```
QUBIT_COORDS … ; R(data); H(deform)
syndrome_meas()                         round 0: perfect reference
TICK
for _ in range(rounds):                 ★multiple rounds (default: code distance d)
    PAULI_CHANNEL_1(data, [px,py,pz])   data noise
    syndrome_meas(flip=p_meas)          ★with measurement flip
    DETECTOR(round_t vs round_{t-1}) × all ancillas
    TICK
data_readout_and_observable()           perfect readout + time-boundary detectors + logical X
```

### Why it looks like this

**Phenomenological adds "measurement errors" on top of "data errors."** Since measurements can no longer
be trusted:

- We must run **multiple rounds** (a single measurement result cannot distinguish a data error from a
  measurement error). The default round count is the code distance `code.distance`
  (the `__init__` in [phenomenological.py](css_xzzx_simulation/src/circuit_builder/phenomenological.py)).
  This is the standard choice of "give the time direction the same distance-d redundancy."
- `syndrome_meas(flip=self.p_meas)` injects the **measurement flip**. `p_meas` defaults to `p`.
- Each round adds a detector between consecutive rounds. A measurement error shows up as a
  "detector that fires in two consecutive rounds" (a time-like edge), while a data error shows up as a
  space-like edge ⇒ both are correctable by matching.

### Time-boundary detectors (`final_boundary_detectors` override)

Immediately after the last round, we "close" the noisy measurements using the **perfect data readout**:

```python
prep_basis = {q: ("Z" if q in hset else "X") for q in data_qubits}   # each data qubit's prep/meas basis
eligible = [anc for anc, legs in stabilizers.items()
            if all(pauli == prep_basis[q] for q, pauli in legs.items())]  # reconstructable checks
for ancilla in eligible:
    targets = [rel(data_record[dc]) for dc in stabilizers[ancilla]]   # check reconstructed from final data
    targets.append(rel(ancilla_record[(ancilla, last)]))             # the same check in the last round
    DETECTOR(targets, ...)
```

Why restrict to `eligible`:
a stabilizer can be **reconstructed from a product of single-qubit measurements** only when the Pauli
type of every one of its legs matches that qubit's measurement basis (`prep_basis`). It is the obvious
constraint that you can only read X-type legs off data measured in the X basis. We rebuild exactly those
checks from the final data readout and compare them against the last round's syndrome, so that
**errors occurring near the final round can also be detected**, closing the time boundary.
(The reason code capacity omits these detectors is as in §2.)


## 4. Circuit-Level Model

[circuit_level.py](css_xzzx_simulation/src/circuit_builder/circuit_level.py)

### Structure (inherits from phenomenological)

`CircuitLevelCircuitBuilder` inherits from `PhenomenologicalCircuitBuilder` and **overrides only
`syndrome_meas`**. The `build` / `__init__` / time-boundary detectors are reused as-is.

Why inheritance suffices:
phenomenological's `build` does, each round,
"`PAULI_CHANNEL_1` on data (= idle/data noise)" → "`syndrome_meas(flip=p_meas)` (= measurement noise)."
Circuit-level is **just that plus reset and gate errors**. So we repurpose the per-round data noise in
`build` as the **idle noise**, and inject only what is missing (reset and two-qubit gate errors) inside
`syndrome_meas`.

### The overridden `syndrome_meas`

```python
def syndrome_meas(self, flip=0.0):
    noisy = self.current_round > 0          # round 0 is the perfect reference
    px, py, pz = biased_pauli_rates(self.p, self.eta)

    self.circuit.append("RX", ancilla_idxs)
    if noisy:
        self.circuit.append("PAULI_CHANNEL_1", ancilla_idxs, [px,py,pz])      # (1) reset error

    for ancilla in ...:
        for dcoord, pauli in legs.items():
            self.circuit.append(CGATE[pauli], [ancilla_idx, data_idx])
            if noisy:
                self.circuit.append("PAULI_CHANNEL_1", [ancilla_idx, data_idx], [px,py,pz])  # (2) 2q gate error

    for ancilla in ...:
        self.circuit.append("MX", [ancilla_idx], flip)   # (3) measurement flip
        ...
```

### Why it looks like this

**Circuit-level is the most realistic model: errors land on every operation (preparation, gates,
measurement, idling).** Here, every noise location uses an **independent biased single-qubit Pauli channel
at rate $(p,\eta)$**:

| Location | Implementation | Intent |
|---|---|---|
| (1) reset error | `PAULI_CHANNEL_1` on ancillas after `RX` | imperfection of the prepared \|+⟩ |
| (2) 2q gate error | `PAULI_CHANNEL_1` on **both qubits independently** after each controlled-Pauli | imperfection of the entangling gate |
| (3) idle/data error | `PAULI_CHANNEL_1` on data each round (**from the inherited build**) | decoherence of data during the round |
| (4) measurement flip | `MX(p_meas)` | readout error |

The two-qubit error is modeled as an **independent biased 1q channel on both qubits** rather than a
**correlated 2q channel (`PAULI_CHANNEL_2`)**, because it lets us reuse the existing `biased_pauli_rates`
directly, keeping the implementation simple and consistent (a deliberate design choice).

- `noisy = self.current_round > 0` keeps **only round 0 perfect**. This avoids mis-deciding whether to
  inject gate errors even for edge values like `p_meas=0` (it does not depend on the value of `flip`).
- The final data readout stays perfect (inherited). The time-boundary detectors are reused as-is from
  phenomenological.



## 5. Model Comparison

| Item | code capacity | phenomenological | circuit-level |
|---|---|---|---|
| Data noise | once | every round | every round (as idle) |
| Measurement noise | none | yes (`p_meas`) | yes (`p_meas`) |
| Reset noise | none | none | yes |
| Gate noise | none | none | yes (per 2q gate) |
| Rounds | effectively 1 | `rounds` (default d) | `rounds` (default d) |
| Time-boundary detectors | none | yes | yes |
| Inherits from | `BaseCircuitBuilder` | `BaseCircuitBuilder` | `PhenomenologicalCircuitBuilder` |

All of them share the skeleton: **X memory**, **round 0 as a perfect reference**,
**consecutive-round comparison detectors**, and **Z-biased noise via `biased_pauli_rates(p, eta)`**.



