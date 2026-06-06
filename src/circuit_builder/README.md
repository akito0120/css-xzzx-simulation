# Circuit Builders — Explained

How the CSS rotated surface code / XZZX code are assembled into Stim circuits under three
noise models (code capacity / phenomenological / circuit-level), and the design rationale behind each.

Files covered:
- [base.py](./base.py) — common base class `BaseCircuitBuilder`
- [shared.py](./shared.py) — model-independent helpers (`biased_pauli_rates`, `biased_two_qubit_rates`, `build_cnot_schedule`, `CGATE`)
- [noisy_measurement.py](./noisy_measurement.py) — intermediate class `NoisyMeasurementCircuitBuilder` (phenomenological + circuit-level)
- [code_capacity.py](./code_capacity.py) — code capacity model
- [phenomenological.py](./phenomenological.py) — phenomenological model
- [circuit_level.py](./circuit_level.py) — circuit-level model

The codes themselves are defined in [code_builder.py](../code_builder.py).

## 0. The Big Picture

All three builders share the same skeleton — an **X-memory experiment** (hold a logical state in the X
basis and measure how well it survives) — and differ **only in how noise is injected**. The hierarchy is
**two-level**; code capacity extends `BaseCircuitBuilder` directly, and the two noisy-measurement models
(phenomenological, circuit-level) extend an intermediate `NoisyMeasurementCircuitBuilder`.

```
BaseCircuitBuilder                       # model-independent plumbing:
│                                          rel, prep_data, qubit init, syndrome_meas (code-cap/pheno only),
│                                          consecutive_round_detectors,
│                                          data_readout + define_observable (final-readout primitives)
├── CodeCapacityCircuitBuilder           # data noise once, perfect measurements
└── NoisyMeasurementCircuitBuilder       # shared by the noisy-measurement models:
      │                                    __init__(rounds, p_meas) + real time-boundary detectors
      ├── PhenomenologicalCircuitBuilder # per-round bulk data noise + measurement flips
      └── CircuitLevelCircuitBuilder     # operation-attached noise: reset/gate/idle + measurement
```

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
| **Round 0 is perfect (noiseless)** | To make the detector "reference point" deterministic. It projects the syndrome into the code space. (Code-capacity / phenomenological only; **circuit-level instead uses a noisy round 0 + `initial_boundary_detectors`** — see §4) |
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

### 1.2 Per-basis preparation — `prep_data` (and `data_basis_partition`)

Each data qubit is prepared **directly in its memory basis** (no Hadamard deformation), matching the
standard surface-code circuit convention (cf. Stim's generated circuits):

```python
x_basis, z_basis = self.data_basis_partition()   # split data qubits by memory basis
self.circuit.append("RX", [idx for _, idx in x_basis])   # |+> off hset (X basis)
self.circuit.append("R",  [idx for _, idx in z_basis])   # |0> on hset (Z basis)
```

- **CSS code**: `hset` is empty → every data qubit is reset with `RX` into $\ket{+}$ (pure X memory).
- **XZZX code**: the qubits in `hset` (a checkerboard subset, see `build_xzzx_code` in
  [code_builder.py](../code_builder.py)) are reset with `R` into $\ket{0}$ (Z basis); the rest with `RX`.

Why: the XZZX code is the **Hadamard-deformed** version of the CSS code (the `deform` function deforms its
stabilizers and logical operators). Rather than realize that deformation with explicit H gates, we simply
**reset each data qubit directly in the basis its deformed frame requires**, so that **the logical X in the
deformed frame is deterministic**. The readout mirrors this: `data_readout` measures the X-basis qubits with
`MX` and the Z-basis qubits with `M` (§1.5). Resetting/measuring in-basis avoids prep/readout H gates
entirely — so there is no separate single-qubit-gate noise location to model there.

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
# data_readout(flip=0.0): measure each data qubit directly in its memory basis
x_basis, z_basis = self.data_basis_partition()
self.circuit.append("MX", x_idxs, flip)    # X-basis qubits (off hset); flip>0 = readout error
self.circuit.append("M",  z_idxs, flip)    # Z-basis qubits (on hset)
# ... record coord -> absolute measurement number (in append order), return it
# flip defaults to 0 (perfect) for code-capacity / phenomenological; circuit-level passes p_meas

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
prep_data()                        RX/R: prepare |+> off hset, |0> on hset
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
QUBIT_COORDS … ; prep_data()  (RX off hset / R on hset)
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
from `Base` (`rel`, `prep_data`, `init_qubit_coords`, `consecutive_round_detectors`, and the
`data_readout` / `define_observable` final-readout primitives). It does **not** inherit from
`PhenomenologicalCircuitBuilder` (see §0 for why).

Its `build` mirrors phenomenological's round/detector loop — but with **both time boundaries noisy**: a
noisy round 0 closed by `initial_boundary_detectors` (instead of a perfect reference round), then the
remaining noisy rounds each emitting a consecutive-round detector, then a final data readout **carrying a
readout error** (`data_readout(flip=p_meas)`, unlike phenomenological's perfect readout) — and with the
usual **one difference**: there is **no bulk per-round data channel**. Circuit-level's data noise comes entirely from the **actual
operations** inside `syndrome_meas`: reset error, two-qubit gate error, and **idle noise during the reset
window, the measurement window, and each of the 4 parallel CNOT steps** (on whichever qubits rest that step).

```python
# circuit-level build: noisy bottom boundary (no perfect round 0), no bulk data channel
self.prep_data()                                     # RX off hset, R on hset (direct basis reset)
# SD6 reset error: |+> (RX) fails by Z, |0> (R) fails by X
self.circuit.append("Z_ERROR", x_idxs, self.p)
self.circuit.append("X_ERROR", z_idxs, self.p)
self.syndrome_meas()                 # round 0 is the first *noisy* round
self.initial_boundary_detectors()                    # deterministic-basis checks vs known +1 (bottom boundary)
for _ in range(self.rounds - 1):                     # remaining rounds; total = self.rounds
    self.current_round += 1
    self.syndrome_meas()             # reset + gate + idle + measurement noise live here
    ... emit consecutive-round detectors ...
```

If `build` *also* applied a phenomenological bulk channel here, it would fire alongside the idle channels
each round and **double-count** the data decoherence — so it deliberately omits it.

Unlike phenomenological, circuit-level has **no perfect reference round**: round 0 is itself noisy, and the
bottom time boundary is closed by `initial_boundary_detectors` — the mirror image of
`final_boundary_detectors` (§3). The deterministic-basis checks (the same `eligible` set, here the X-type
checks for X-memory) have a known $+1$ eigenvalue in round 0, so a single-measurement detector on their
round-0 syndrome catches the reset error (the per-basis prep error: `Z_ERROR(p)` on the `RX`/$\ket{+}$
qubits, `X_ERROR(p)` on the `R`/$\ket{0}$ qubits) and any round-0 measurement error. Naively adding reset noise *before* a perfect
round 0 would instead make it **undetectable** (absorbed into the reference, flipping only the observable)
and collapse the effective distance to $\sim1$ — which is why round 0 is left noisy rather than perfect.

### The overridden `syndrome_meas`

```python
def syndrome_meas(self):
    # every round (round 0 included) is noisy; all rates read from self (no flip arg);
    # the bottom boundary is closed by initial_boundary_detectors
    px, py, pz = biased_pauli_rates(self.p, self.eta)
    pc2 = biased_two_qubit_rates(self.p, self.eta)   # 15 probs for PAULI_CHANNEL_2

    self.circuit.append("RX", ancilla_idxs)
    self.circuit.append("Z_ERROR", ancilla_idxs, self.p)                  # (1) reset error: Z flip of |+>
    self.circuit.append("PAULI_CHANNEL_1", data_list,    [px,py,pz])      # (3a) idle: data waits during reset

    schedule = build_cnot_schedule(self.code)   # ancilla -> [leg per step], 4 parallel steps
    for step in range(4):
        busy = set()
        for ancilla in self.ancilla_order:
            dcoord = schedule[ancilla][step]
            if dcoord is None:                  # this ancilla rests this step (boundary check)
                continue
            pauli = self.code.stabilizers[ancilla][dcoord]
            self.circuit.append(CGATE[pauli], [ancilla_idx, data_idx])
            # (2) 2q gate error: same biased correlated channel for both CX and CZ
            self.circuit.append("PAULI_CHANNEL_2", [ancilla_idx, data_idx], pc2)
            busy |= {ancilla_idx, data_idx}
        idle = [q for q in all_qubits if q not in busy]   # (3c) per-step idle: whoever rests this step
        self.circuit.append("PAULI_CHANNEL_1", idle, [px,py,pz])
        self.circuit.append("TICK")

    self.circuit.append("PAULI_CHANNEL_1", data_list, [px,py,pz])         # (3b) idle: data waits during measure

    for ancilla in ...:
        self.circuit.append("MX", [ancilla_idx], self.p_meas)   # (4) measurement flip
        ...
```

The controlled gates are emitted **step by step** rather than ancilla by ancilla. Each of the 4
steps applies one leg from every stabilizer that still has work to do, all in parallel, so a round
has a well-defined depth and a well-defined notion of "which qubits are idle right now."

### Why it looks like this

**Circuit-level is the most realistic model: errors land on every operation (preparation, gates,
measurement, idling).** The **idle** locations use a **biased single-qubit Pauli channel at rate
$(p,\eta)$**, the **ancilla reset** uses a **pure $Z$ flip** at rate $p$ (the only fault that affects an
$\ket{+}$ preparation), and the two-qubit gate uses a **biased correlated two-qubit Pauli channel**:

| Location | Implementation | Intent |
|---|---|---|
| (1) reset error | `Z_ERROR(p)` on ancillas after `RX` | imperfection of the prepared $\ket{+}$: a $Z$ flip sends $\ket{+}\to\ket{-}$. A full Pauli channel would be wasteful here since $X\ket{+}=\ket{+}$ |
| (2) 2q gate error | after **both CZ and CX**: `PAULI_CHANNEL_2` (the same biased correlated 2-qubit channel) | imperfection of the entangling gate — uniformly Z-biased, assuming bias-preserving gates (biased SD6) |
| (3a) idle error (reset window) | `PAULI_CHANNEL_1` on **all data** right after the ancillas are reset | data decoheres while it waits for the ancillas to be reset |
| (3b) idle error (measure window) | `PAULI_CHANNEL_1` on **all data** right before the ancillas are measured | data decoheres while it waits for the ancillas to be measured |
| (3c) idle error (per CNOT step) | `PAULI_CHANNEL_1` on the qubits **not gated** in that step (boundary data + resting ancillas) | data/ancillas decohere while they sit out a gate step they have no leg in |
| (4) measurement flip | `MX(p_meas)` | readout error |

The two-qubit gate error is a genuine **correlated** channel (each gate fault produces a weight-2 error at
$O(p)$, as in standard circuit-level models), built by `biased_two_qubit_rates` in
[shared.py](./shared.py) next to `biased_pauli_rates`. It follows the **bias-preserving-gate** convention of
the XZZX biased-noise literature (Darmawan *et al.*, [arXiv:2104.09539](https://arxiv.org/abs/2104.09539),
XZZX + Kerr-cat qubits): the 15 non-identity two-qubit Paulis are
**partitioned** into a high-rate Z-subgroup $H=\{IZ, ZI, ZZ\}$ and the 12 remaining low-rate errors, with
bias $\eta = P(H)/P(L)$ — exactly the two-qubit generalization of the single-qubit convention
$\eta = p_Z/(p_X+p_Y)$ (high-rate $\{Z\}$ vs low-rate $\{X,Y\}$). Probabilities are uniform within each set
and sum to $p$: each high-rate component $= p\eta/(3(1+\eta))$, each low-rate component $= p/(12(1+\eta))$.
As $\eta\to\infty$ the channel concentrates **uniformly on the Z-subgroup** $\{IZ, ZI, ZZ\}$ (each $p/3$ — a
single $Z$ on *either* qubit survives, the physically correct high-bias limit). Note the standard two-qubit
depolarizing point (all 15 $=p/15$) is at $\eta=1/4$, **not** $1/2$, because $|H|=3$ and $|L|=12$ (unlike the
single-qubit $1$-vs-$2$ split, which makes $\eta=1/2$ depolarizing); prior research defines bias as the
high/low probability ratio rather than anchoring on a depolarizing point. The distribution is symmetric in
the two qubits, so the control/target order is irrelevant.

This is a **uniform biased SD6** model: the same biased channel is applied after **both** the CZ (Z-type
checks) and CX (X-type checks) gates, i.e. **all** two-qubit gates are assumed **bias-preserving**. This is
physically the **bias-preserving-hardware** regime (e.g. cat qubits, where the CX can be made
bias-preserving), *not* a generic two-level-qubit device — on two-level qubits a no-go theorem forbids a
bias-preserving CNOT, so a faithful transmon model would instead give the CX plain depolarizing noise (the
"HBD hybrid" of [arXiv:2505.17718](https://arxiv.org/abs/2505.17718)). The choice matters: under the uniform
biasing here, XZZX's advantage **keeps growing with $\eta$**, whereas an unbiased/depolarizing CX would
reintroduce $X/Y$ noise every round and **cap** the advantage (the threshold would saturate around $\sim1\%$
as $\eta\to\infty$). To model the hybrid instead, branch on the leg's Pauli and emit `DEPOLARIZE2(p)` for the
CX legs. These codes only ever emit CX and CZ (stabilizers are X/Z only), so CY never arises.
(Building the decoder DEM from `PAULI_CHANNEL_2` requires `approximate_disjoint_errors=True`, which
[simulation.py](../simulation.py) passes.)

**Why idle noise, not a bulk channel.**
In a real device the data qubits are not actually idle during a round: they are repeatedly entangled by
the gates (2), and in between they *wait*. There are three kinds of waiting period, each getting its own
idle channel: while every ancilla is reset (3a), while every ancilla is measured (3b), and — because the
entangling gates run in **4 parallel steps** — while a qubit *sits out a gate step it has no leg in* (3c).
The per-step idle (3c) only ever lands on **boundary** qubits: a weight-4 bulk plaquette engages all four
of its data qubits across the four steps, but weight-2/3 boundary checks leave their ancilla (and the
unused data sites) resting in some steps. This is what makes circuit-level "operation-attached": each
channel (1)–(4) maps to a concrete hardware step, with no leftover abstract bulk term. (Where Stim sees
adjacent identical noise instructions on the same targets it may fuse them in the printed circuit, but
that is purely cosmetic.)

**The CNOT schedule and distance preservation.**
The four steps come from `build_cnot_schedule` ([shared.py](./shared.py)), which assigns every
stabilizer's legs to time steps by their lattice offset (`STEP_OF_OFFSET`). Two properties matter:

- **No double-booking.** Within a step, no data qubit is targeted by two ancillas — a physical
  requirement, since a qubit can only take part in one two-qubit gate at a time. The offset-based
  assignment guarantees this for the rotated-surface-code geometry (a data qubit's four neighbouring
  ancillas sit at four distinct offsets, hence four distinct steps).
- **Distance preservation.** A careless gate order lets a single mid-round ancilla fault propagate into a
  weight-2 *hook error* aligned with a logical operator, which would halve the effective distance
  ($d \to \lceil (d{+}1)/2 \rceil$). The chosen row-major order keeps hook errors off the logical-X
  direction, so the effective fault distance stays exactly $d$.

A single common order works here because the experiment measures **only the logical-X observable**; a
two-sided memory (X *and* Z) would need the standard X/Z-transposed schedule instead (noted in the
[shared.py](./shared.py) comment).

- **Every round, including round 0, is noisy** — `syndrome_meas` injects its noise unconditionally.　The bottom boundary is instead closed by
  `initial_boundary_detectors` together with the `X_ERROR(p)` reset error, so reset / round-0 faults are
  detectable (see above). This keeps both time boundaries faithful.
- The final data readout **carries a readout error** at rate `p_meas` (via `Base.data_readout(flip=p_meas)`
  + `Base.define_observable`, composed directly in `build`) — so the upper time boundary is noisy too, like
  a standard circuit-level memory experiment. The `final_boundary_detectors` pick these readout flips up as
  detection events, and they also feed the logical-X observable. Between `data_readout` and
  `define_observable`, `build` calls the time-boundary `final_boundary_detectors` from
  `NoisyMeasurementCircuitBuilder` — the same prep-basis reconstruction described in §3. Because the logic
  is common to both noisy-measurement models, phenomenological and circuit-level **call that one
  intermediate-class implementation** rather than each carrying a copy.



## 5. Model Comparison

| Item | code capacity | phenomenological | circuit-level |
|---|---|---|---|
| Data noise | once | bulk, every round (`data_round_noise`) | none (bulk disabled) |
| Idle noise | none | none | every round, reset + measure windows + per CNOT step |
| Measurement noise | none | yes (`p_meas`) | yes (`p_meas`) |
| Reset noise | none | none | yes |
| Gate noise | none | none | yes (per 2q gate) |
| Rounds | effectively 1 | `rounds` (default d) | `rounds` (default d) |
| Time-boundary detectors | none | yes | yes |
| Inherits from | `BaseCircuitBuilder` | `NoisyMeasurementCircuitBuilder` | `NoisyMeasurementCircuitBuilder` |

All of them share the skeleton: **X memory**, **consecutive-round comparison detectors**, and
**Z-biased noise via `biased_pauli_rates(p, eta)`**. Code-capacity and phenomenological use a **perfect
reference round 0**; circuit-level instead makes **both boundaries noisy** (noisy round 0 +
`initial_boundary_detectors`, noisy final readout), the faithful circuit-level convention.

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
- **Circuit-level**: `p` is the **per-operation** error rate. Reset, each idle window (reset / measure),
  and the measurement each fail at rate ≈ `p` (single-qubit channels); every two-qubit gate fails at total
  rate `p` via the **correlated** `PAULI_CHANNEL_2` (one weight-≤2 Pauli fault per gate). A data qubit
  therefore passes through several fault locations per round, so its **effective per-round error is
  larger than `p`** — exactly the standard circuit-level convention (and the reason circuit-level
  thresholds in `p` are numerically far below phenomenological ones).

Consequences:

1. **Do not compare the models at the same numeric `p`.** Report each threshold on its own model's
   `p`-axis; a smaller circuit-level threshold is expected, not a regression.
2. **Within one model, the `(p, eta)` sweep is self-consistent**, so the CSS-vs-XZZX comparison — the
   actual goal — is valid because it is always done at the same `p` *and* the same model.
3. Applying the idle channel in the reset window (3a), the measure window (3b) **and** each CNOT step
   (3c) is **not** double-counting: they are distinct physical waiting periods, and within a step (3c)
   only covers qubits with no gate that step.
4. **Per-step idling is modelled.** The entangling gates are scheduled into 4 parallel steps
   (`build_cnot_schedule`), and boundary qubits that rest during a step receive idle noise (3c) on top of
   the reset/measure windows. Because idle noise here is **biased**, its placement is not a neutral knob:
   it sets how much biased noise the data accumulates per round, which directly affects the high-bias XZZX
   threshold — so the per-step model matters for the CSS-vs-XZZX comparison rather than merely shifting
   both codes equally.



