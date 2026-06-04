# Circuit Builders — Explained

How the CSS rotated surface code / XZZX code are assembled into Stim circuits under three
noise models (code capacity / phenomenological / circuit-level), and the design rationale behind each.

Files covered:
- [base.py](./base.py) — common base class `BaseCircuitBuilder`
- [shared.py](./shared.py) — model-independent helpers (`biased_pauli_rates`, `CGATE`)
- [noisy_measurement.py](./noisy_measurement.py) — intermediate class `NoisyMeasurementCircuitBuilder` (phenomenological + circuit-level)
- [code_capacity.py](./code_capacity.py) — code capacity model
- [phenomenological.py](./phenomenological.py) — phenomenological model
- [circuit_level.py](./circuit_level.py) — circuit-level model

The codes themselves are defined in [code_builder.py](../code_builder.py).

## 0. The Big Picture

All three builders share the same skeleton — an **X-memory experiment** (hold a logical state in the X
basis and measure how well it survives) — and differ **only in how noise is injected**. The hierarchy is
**two-level** and the key invariant is that **no model inherits from a sibling model** — only from an
abstract base. Code capacity extends `BaseCircuitBuilder` directly; the two noisy-measurement models
(phenomenological, circuit-level) extend an intermediate `NoisyMeasurementCircuitBuilder`.

```
BaseCircuitBuilder                       # model-independent plumbing:
│                                          rel, deform, qubit init, syndrome_meas (code-cap/pheno only),
│                                          consecutive_round_detectors,
│                                          data_readout + define_observable (final-readout primitives)
├── CodeCapacityCircuitBuilder           # data noise once, perfect measurements
└── NoisyMeasurementCircuitBuilder       # shared by the noisy-measurement models:
      │                                    __init__(rounds, p_meas) + real time-boundary detectors
      ├── PhenomenologicalCircuitBuilder # per-round bulk data noise + measurement flips
      └── CircuitLevelCircuitBuilder     # operation-attached noise: reset/gate/idle + measurement
```

Why the intermediate class (and not `CircuitLevel ← Phenomenological`):
`CircuitLevel` used to inherit from `Phenomenological` to reuse its round/detector loop, but that made
circuit-level's correctness depend on a *sibling* model's internals (a *fragile base class*). The fix is an
**abstract intermediate**, not sibling inheritance: `NoisyMeasurementCircuitBuilder` holds only what the two
noisy-measurement models genuinely share and code capacity does not need — the `__init__(rounds, p_meas)`
constructor and the time-boundary `final_boundary_detectors`. Since **code-capacity and phenomenological are
finished/frozen and only circuit-level keeps evolving**, circuit-level still **owns its own `build` and
`syndrome_meas`** (the model-specific parts) and merely *inherits* the stable shared constructor and
boundary detectors. It can change its gate scheduling, idle model, etc. with zero risk to the frozen models.

Where the data noise lives differs between the two multi-round models:

- **Phenomenological** lumps each round's data decoherence into a single **bulk** `PAULI_CHANNEL_1` on
  all data qubits (the `data_round_noise()` helper).
- **Circuit-level** has **no bulk channel at all**; its data decoherence comes from the reset, two-qubit
  gate, and **idle** errors injected at the actual operations inside `syndrome_meas`.

This is why circuit-level is a faithful "every operation is noisy" model and not just phenomenological
with extra terms. Adding a bulk channel *on top of* the idle noise would **double-count** the data noise,
which is exactly why circuit-level's own `build` omits it. Details in §3–§4.

The "why" shared by all three models:

| Design choice | Reason |
|---|---|
| **X-memory basis** | The goal is to show XZZX's advantage under Z-biased noise. We hold a logical X and read out in the X basis at the end |
| **Round 0 is perfect (noiseless)** | To make the detector "reference point" deterministic. It projects the syndrome into the code space |
| **Detector = comparison of consecutive rounds' syndromes** | If nothing happened, the syndrome is unchanged. A change = an error occurred in between — a difference detector |
| **Single ancilla + controlled-Pauli** | Prepare one ancilla in $\ket{+}$, apply controlled Paulis to the data, and measure it in X — this measures any Pauli stabilizer |


## 1. Shared Building Blocks (BaseCircuitBuilder)

The common methods in [base.py](./base.py). (The model-independent helpers `biased_pauli_rates` and
`CGATE` live in [shared.py](./shared.py).)

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

- **CSS code**: `hset` is empty. H is applied to every data qubit, preparing them all in $\ket{+}$ (pure X memory).
- **XZZX code**: the qubits in `hset` (a checkerboard subset of data qubits, see `build_xzzx_code` in
  [code_builder.py](../code_builder.py)) are **not** given an H → they stay in
  $\ket{0}$ (the Z basis).

Why: the XZZX code is the **Hadamard-deformed** version of the CSS code, obtained by applying H to the
data qubits in `hset`. The stabilizers and logical operators undergo the same deformation (the `deform`
function). On the circuit side, we choose each data qubit's preparation basis ($\ket{+}$ or $\ket{0}$) to match
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
put the ancilla in $\ket{+}$ (an X eigenstate) and apply controlled Paulis targeting the data; the ancilla's
phase flips according to the eigenvalue of the stabilizer operator $S$. Measuring the ancilla in the X
basis (`MX`) then yields the value of $S$. `CGATE = {"X":"CX","Y":"CY","Z":"CZ"}` selects the controlled
gate matching each leg's Pauli type ([shared.py](./shared.py)).

The `flip` argument is the **injection point for measurement noise**. Switching it on/off is all it takes
to toggle between "perfect" and "noisy" measurements, which is why **code-capacity and phenomenological**
build on this same `syndrome_meas`. **Circuit-level is the exception**: it overrides `syndrome_meas` to
attach per-operation reset/gate/idle noise that a `flip` toggle cannot express (see §4).

### 1.4 Consecutive-round detectors — `consecutive_round_detectors`

```python
for ancilla in self.ancilla_order:
    d_now  = self.ancilla_record[(ancilla, self.current_round)]
    d_prev = self.ancilla_record[(ancilla, self.current_round - 1)]
    self.circuit.append("DETECTOR", [self.rel(d_now), self.rel(d_prev)], [*ancilla, self.current_round])
```

Emits one detector per ancilla comparing **this round's syndrome against the previous round's**. If no
error occurred in between, the two measurements agree and the detector stays silent; a flip signals an
error in that interval (the difference-detector idea from §0). It reads `self.current_round` at call time,
so the **same method serves every model**: code capacity calls it once (round 1 vs round 0), the two
multi-round models call it each round inside their loop. Being model-independent, it lives on `Base`
alongside `rel` / `syndrome_meas`.

### 1.5 Final readout and logical observable — `data_readout` + `define_observable`

`Base` provides two single-responsibility primitives; there is **no composite method**. Each `build()`
composes them explicitly (and the noisy models slot `final_boundary_detectors` in between):

```python
# data_readout(): final perfect readout, returns coord -> absolute measurement index
self.deform_x_basis_data()                 # return to the same basis as preparation
self.circuit.append("M", data_list)        # measure all data in Z (effectively the X-memory readout)
# ... record coord -> absolute measurement number, return it

# define_observable(data_record): the logical-X observable
observable_targets = [self.rel(data_record[dc]) for dc in self.code.logical_x]
self.circuit.append("OBSERVABLE_INCLUDE", observable_targets, 0)
```

A `build()` therefore ends with:

```python
data_record = self.data_readout()
self.final_boundary_detectors(data_record)   # noisy models only — omitted by code capacity
self.define_observable(data_record)
```

- The logical-X observable is defined by the parity of the final measurements of the data qubits listed
  in `code.logical_x`.
- `final_boundary_detectors` adds **time-boundary detectors** that cross-check "the stabilizer
  reconstructed from the final data readout" against "the last round's syndrome." It lives **only** in
  `NoisyMeasurementCircuitBuilder` (so it is **shared by phenomenological and circuit-level**) and `Base`
  knows nothing about it. The two noisy models call it explicitly in their `build()`; code capacity, whose
  final readout is noiseless, simply omits the call — those detectors would only be redundant (see §2).


## 2. Code Capacity Model

[code_capacity.py](./code_capacity.py)

### Structure

```
QUBIT_COORDS …                     declare qubit coordinates
R(data); H(deform)                 prepare |+> / |0>
syndrome_meas()                    round 0: perfect reference (project into the code space)
TICK
PAULI_CHANNEL_1(data, [px,py,pz])  ▶ the only noise: one biased Pauli channel on data
syndrome_meas()                    round 1: perfect measurement
DETECTOR(round1 vs round0) × all ancillas
data_readout(); define_observable() perfect readout + logical X (no boundary detectors)
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
  be redundant). Code capacity's `build()` therefore calls just `data_readout` + `define_observable` and
  **omits the `final_boundary_detectors` call** — that logic exists only in `NoisyMeasurementCircuitBuilder`,
  which code capacity does not extend.

The detectors compare round 1 against round 0 because (as in §0) we look at the difference of a syndrome
that should be unchanged in the absence of noise. Since round 0 establishes a perfect baseline, the
detectors that fire in round 1 are exactly the traces of the errors introduced by `PAULI_CHANNEL_1`.


## 3. Phenomenological Model

[phenomenological.py](./phenomenological.py)

### Structure

```
QUBIT_COORDS … ; R(data); H(deform)
syndrome_meas()                         round 0: perfect reference
TICK
for _ in range(rounds):                 ▶ multiple rounds (default: code distance d)
    data_round_noise()                  bulk data noise (PAULI_CHANNEL_1 on all data)
    syndrome_meas(flip=p_meas)          ▶ with measurement flip
    DETECTOR(round_t vs round_{t-1}) × all ancillas
    TICK
data_readout()                          perfect readout
final_boundary_detectors()              time-boundary detectors
define_observable()                     logical X
```

### Why it looks like this

**Phenomenological adds "measurement errors" on top of "data errors."** Since measurements can no longer
be trusted:

- We must run **multiple rounds** (a single measurement result cannot distinguish a data error from a
  measurement error). The default round count is the code distance `code.distance`
  (the `__init__` in [phenomenological.py](./phenomenological.py)).
  This is the standard choice of "give the time direction the same distance-d redundancy."
- `data_round_noise()` injects the **bulk data noise**: one biased `PAULI_CHANNEL_1` over *all* data
  qubits per round. At the phenomenological level we do not model individual gates, so every source of
  data decoherence during a round is **lumped into this single channel** placed just before the
  measurements. Circuit-level does not use a bulk channel at all — it models that decoherence through
  per-operation gate and idle errors instead (§4).
- `syndrome_meas(flip=self.p_meas)` injects the **measurement flip**. `p_meas` defaults to `p`.
- Each round adds a detector between consecutive rounds. A measurement error shows up as a
  "detector that fires in two consecutive rounds" (a time-like edge), while a data error shows up as a
  space-like edge ⇒ both are correctable by matching.

### Time-boundary detectors (`final_boundary_detectors`, defined in `NoisyMeasurementCircuitBuilder`)

Immediately after the last round, we "close" the noisy measurements using the **perfect data readout**.
This logic is common to both noisy-measurement models, so it lives in the intermediate
`NoisyMeasurementCircuitBuilder`; each model calls it explicitly in its `build()`, between `data_readout`
and `define_observable` (shared verbatim with circuit-level, §4):

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
(Code capacity does not extend this class and simply omits the call in its `build()`, skipping these
detectors — see §2.)


## 4. Circuit-Level Model

[circuit_level.py](./circuit_level.py)

### Structure (owns its model-specific parts, extends `NoisyMeasurementCircuitBuilder`)

`CircuitLevelCircuitBuilder` extends `NoisyMeasurementCircuitBuilder` and owns its **model-specific** parts
— `build` and `syndrome_meas` — while *inheriting* the constructor `__init__(rounds, p_meas)` and the
time-boundary `final_boundary_detectors` from that intermediate class, plus the model-independent plumbing
from `Base` (`rel`, `deform_x_basis_data`, `init_qubit_coords`, `consecutive_round_detectors`, and the
`data_readout` / `define_observable` final-readout primitives). It does **not** inherit from
`PhenomenologicalCircuitBuilder` (see §0 for why).

Its `build` mirrors phenomenological's round/detector loop — round 0 reference, then `rounds` noisy
rounds each emitting a consecutive-round detector, then a perfect final readout — with **one difference**:
there is **no bulk per-round data channel**. Circuit-level's data noise comes entirely from the **actual
operations** inside `syndrome_meas`: reset error, two-qubit gate error, and **idle noise on data during
the reset and measurement windows**.

```python
# circuit-level build's noisy round — no bulk data channel, just syndrome_meas:
for _ in range(self.rounds):
    self.current_round += 1
    self.syndrome_meas(flip=self.p_meas)   # reset + gate + idle + measurement noise live here
    ... emit consecutive-round detectors ...
```

If `build` *also* applied a phenomenological bulk channel here, it would fire alongside the idle channels
each round and **double-count** the data decoherence — so it deliberately omits it.

### The overridden `syndrome_meas`

```python
def syndrome_meas(self, flip=0.0):
    noisy = self.current_round > 0          # round 0 is the perfect reference
    px, py, pz = biased_pauli_rates(self.p, self.eta)

    self.circuit.append("RX", ancilla_idxs)
    if noisy:
        self.circuit.append("PAULI_CHANNEL_1", ancilla_idxs, [px,py,pz])      # (1) reset error
        self.circuit.append("PAULI_CHANNEL_1", data_list,    [px,py,pz])      # (3a) idle: data waits during reset

    for ancilla in ...:
        for dcoord, pauli in legs.items():
            self.circuit.append(CGATE[pauli], [ancilla_idx, data_idx])
            if noisy:
                self.circuit.append("PAULI_CHANNEL_1", [ancilla_idx, data_idx], [px,py,pz])  # (2) 2q gate error

    if noisy:
        self.circuit.append("PAULI_CHANNEL_1", data_list, [px,py,pz])         # (3b) idle: data waits during measure

    for ancilla in ...:
        self.circuit.append("MX", [ancilla_idx], flip)   # (4) measurement flip
        ...
```

### Why it looks like this

**Circuit-level is the most realistic model: errors land on every operation (preparation, gates,
measurement, idling).** Here, every noise location uses an **independent biased single-qubit Pauli channel
at rate $(p,\eta)$**:

| Location | Implementation | Intent |
|---|---|---|
| (1) reset error | `PAULI_CHANNEL_1` on ancillas after `RX` | imperfection of the prepared $\ket{+}$ |
| (2) 2q gate error | `PAULI_CHANNEL_1` on **both qubits independently** after each controlled-Pauli | imperfection of the entangling gate |
| (3a) idle error (reset window) | `PAULI_CHANNEL_1` on **all data** right after the ancillas are reset | data decoheres while it waits for the ancillas to be reset |
| (3b) idle error (measure window) | `PAULI_CHANNEL_1` on **all data** right before the ancillas are measured | data decoheres while it waits for the ancillas to be measured |
| (4) measurement flip | `MX(p_meas)` | readout error |

The two-qubit error is modeled as an **independent biased 1q channel on both qubits** rather than a
**correlated 2q channel (`PAULI_CHANNEL_2`)**, because it lets us reuse the existing `biased_pauli_rates`
directly, keeping the implementation simple and consistent (a deliberate design choice).

**Why idle noise replaces the bulk channel.**
In a real device the data qubits are not actually idle during a round: they are repeatedly entangled by
the gates (2), and in between they *wait* — first while every ancilla is reset, then while every ancilla
is measured. Those two waiting windows are the genuine "idling" periods, so the idle noise is injected
exactly there (3a, 3b), once on every data qubit per window. This is what makes circuit-level
"operation-attached": each channel (1)–(4) maps to a concrete hardware step, with no leftover abstract
bulk term. (Adjacent identical `PAULI_CHANNEL_1` instructions are fused by Stim — e.g. reset error + 3a
appear as one combined channel — but that is purely cosmetic.)

- `noisy = self.current_round > 0` keeps **only round 0 perfect** (the reference round gets no reset,
  gate, or idle noise). This decision does not depend on the value of `flip`, so it stays correct even
  for edge cases like `p_meas=0`.
- The final data readout stays perfect (via `Base.data_readout` + `Base.define_observable`, composed
  directly in `build`). Between them, `build` calls the time-boundary `final_boundary_detectors` from
  `NoisyMeasurementCircuitBuilder` — the same prep-basis reconstruction described in §3. Because the logic
  is common to both noisy-measurement models, phenomenological and circuit-level **call that one
  intermediate-class implementation** rather than each carrying a copy.



## 5. Model Comparison

| Item | code capacity | phenomenological | circuit-level |
|---|---|---|---|
| Data noise | once | bulk, every round (`data_round_noise`) | none (bulk disabled) |
| Idle noise | none | none | every round, reset + measure windows |
| Measurement noise | none | yes (`p_meas`) | yes (`p_meas`) |
| Reset noise | none | none | yes |
| Gate noise | none | none | yes (per 2q gate) |
| Rounds | effectively 1 | `rounds` (default d) | `rounds` (default d) |
| Time-boundary detectors | none | yes | yes |
| Inherits from | `BaseCircuitBuilder` | `NoisyMeasurementCircuitBuilder` | `NoisyMeasurementCircuitBuilder` |

All of them share the skeleton: **X memory**, **round 0 as a perfect reference**,
**consecutive-round comparison detectors**, and **Z-biased noise via `biased_pauli_rates(p, eta)`**.

### On the meaning of `p` across models (important)

`p` is **not** a single universal physical error rate shared by the three models — each model defines it
in its own operational terms, so **the same `(p, eta)` produces different effective error rates**, and
this is intentional, not a bug:

- **Code capacity**: `p` is the probability that a data qubit takes a Pauli error in the **single**
  noisy channel applied once (rates summing to `p`); gates and measurements are perfect. With no
  measurement noise and effectively one round, `p` is just the "per-data-qubit error of a single shot."
- **Phenomenological**: `p` is the probability that a data qubit takes a Pauli error **per round**
  (one lumped `data_round_noise()` channel, rates summing to `p`); the **same `p` also sets the
  measurement flip** (`p_meas` defaults to `p`).
- **Circuit-level**: `p` is the **per-operation** error rate. Reset, every two-qubit gate, each idle
  window (reset / measure), and the measurement *each* fail independently at rate ≈ `p`. A data qubit
  therefore passes through several fault locations per round, so its **effective per-round error is
  larger than `p`** — exactly the standard circuit-level convention (and the reason circuit-level
  thresholds in `p` are numerically far below phenomenological ones).

Consequences:

1. **Do not compare the models at the same numeric `p`.** Report each threshold on its own model's
   `p`-axis; a smaller circuit-level threshold is expected, not a regression.
2. **Within one model, the `(p, eta)` sweep is self-consistent**, so the CSS-vs-XZZX comparison — the
   actual goal — is valid because it is always done at the same `p` *and* the same model.
3. Applying the idle channel in **both** the reset window (3a) and the measure window (3b) is **not**
   double-counting: they are two distinct physical waiting periods. (The earlier bulk-vs-idle
   double-count, fixed by disabling `data_round_noise` in circuit-level, was a separate issue.)
4. **Known simplification**: idle noise is injected only in the reset/measure windows, not between the
   CX layers (no per-layer idling). This slightly *over*-estimates the threshold but affects CSS and
   XZZX equally, so it does not bias the comparison.



