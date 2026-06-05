# Comparison of Rotated Surface Codes and XZZX Codes under Biased Noise Models

This project simulates and compares the **rotated CSS surface code** against the
**XZZX surface code** under Z-biased Pauli noise, and estimates their error-correction
**thresholds** via finite-size scaling. The motivation is that the XZZX code — a
Hadamard-deformed variant of the surface code — is known to perform better than the
plain CSS surface code when the physical noise is biased toward dephasing (Z) errors,
as is common in many real qubit platforms.

The simulation is built on [Stim](https://github.com/quantumlib/Stim) for fast stabilizer
circuit sampling and [PyMatching](https://github.com/oscarhiggott/PyMatching) for
minimum-weight perfect matching decoding.

## Project Structure

```
css_xzzx_simulation/
├── README.md                      # This file
├── src/
│   ├── main.py                    # Entry point: sweep, decode, FSS fit, plot, diagrams
│   ├── code_builder.py            # Code definitions (rotated surface + XZZX)
│   ├── circuit_builder/           # Noise-model circuit builders (see its own README)
│   │   ├── README.md              # Detailed explanation of the three noise models
│   │   ├── shared.py              # BaseCircuitBuilder + shared helpers
│   │   ├── code_capacity.py       # Code capacity model
│   │   ├── phenomenological.py    # Phenomenological model
│   │   └── circuit_level.py       # Circuit-level model
│   ├── simulation.py              # Monte Carlo logical error rate estimation
│   └── threshold.py               # Finite-size scaling (FSS) threshold fit
└── results/                       # Generated plots and circuit diagrams (git-ignored)
```

## The Codes

Both codes are defined in [code_builder.py](./src/code_builder.py) and share the common
`QecCode` dataclass, which records the data/ancilla qubit layout on the lattice, the
stabilizer legs (ancilla → {data qubit → Pauli type}), and the logical X/Z operators.

### Rotated surface code (`build_rotated_surface_code`)

A distance-$d$ rotated surface code. Data qubits sit at odd-odd lattice coordinates
$(2c+1, 2r+1)$ and plaquette centres (ancillas) at even-even coordinates. Stabilizer
type is assigned by a checkerboard $(r+c)$ parity (X vs Z plaquettes). Weight-2 boundary
checks are kept only on the appropriate edges (Z on top/bottom, X on left/right) following
the standard rotated-code boundary convention. The logical X is a horizontal string and
the logical Z a vertical string of data qubits.

### XZZX code (`build_xzzx_code`)

The XZZX code is constructed as the **Hadamard-deformed** version of the rotated surface
code. A checkerboard subset of data qubits (`hset`) has H applied conjugation-wise, which
maps each affected leg via $X \leftrightarrow Z$, $Y \leftrightarrow Y$ (`deform`). As a
result, every stabilizer becomes a mixed-Pauli check of the form $X\,Z\,Z\,X$ around a
plaquette. The stabilizers and the logical operators are deformed consistently, and the
set `hset` is carried in `QecCode` so the circuit builders know each qubit's effective
basis. For the plain CSS code `hset` is empty, so all the deformation logic degrades to
the identity.

## Noise Model and Bias Parameterization

Biased Pauli noise is parameterized by the bias

$$
\begin{align*}
& \eta = \frac{p_Z}{p_X + p_Y} \\
& p_X = p_Y \\
& p = p_X + p_Y + p_Z \\ \\
\end{align*}
$$

so that

$$
\begin{align*}
& p_X = p_Y = \frac{p}{2 \cdot (1 + \eta)} \\
& p_Z = \frac{p \cdot \eta}{1 + \eta}
\end{align*}
$$

$\eta = 0.5$ is standard depolarizing noise; $\eta \rightarrow \infty$ is pure dephasing
(Z only). $p$ is the total single-qubit error probability. This mapping is implemented by
`biased_pauli_rates(p, eta)` in [circuit_builder/shared.py](./src/circuit_builder/shared.py)
(the $\eta = \infty$ case returns the pure-Z channel $(0, 0, p)$).

The circuit-level model is a **biased SD6** model: the standard depolarizing SD6 circuit model (single
parameter $p$, every operation noisy) with each channel made **Z-biased** by the factor $\eta$. It assumes
**bias-preserving two-qubit gates** (as realizable on e.g. cat qubits — Darmawan *et al.*,
[arXiv:2104.09539](https://arxiv.org/abs/2104.09539)), so **both** the CZ (Z-type checks) and CX (X-type
checks) gates get the **same biased correlated two-qubit channel**: the 15 two-qubit Paulis split into a
high-rate Z-subgroup $\{IZ, ZI, ZZ\}$ and 12 low-rate errors, with bias
$\eta = P(\{IZ,ZI,ZZ\})/P(\text{rest})$, generalizing the single-qubit $\eta = p_Z/(p_X+p_Y)$; as
$\eta\to\infty$ it concentrates uniformly on $\{IZ, ZI, ZZ\}$ (implemented by `biased_two_qubit_rates(p, eta)`).
This uniform biasing is what lets XZZX's advantage keep growing with $\eta$ (rather than saturating, as it
would under an unbiased/depolarizing CX on two-level qubits). See
[circuit_builder/README.md](./src/circuit_builder/README.md) for the full convention and limits.

### Noise models

Three noise models are provided as separate circuit builders, all extending a common
`BaseCircuitBuilder` and differing only in how noise is injected:

| Model | Data noise | Measurement / gate / idle noise | Rounds |
|---|---|---|---|
| **Code capacity** | one biased Pauli channel on data | perfect | effectively 1 |
| **Phenomenological** | bulk channel every round | measurement flips | $d$ (default) |
| **Circuit-level** | from operations only | reset + 2q gate + per-step idle + measurement | $d$ (default) |

For a full explanation of all three models, the shared building blocks, and the design
rationale behind them, see the dedicated
[circuit_builder/README.md](./src/circuit_builder/README.md).

## Memory Basis

The experiments run an X-basis memory (they track the logical X operator), which is
sensitive to the dominant Z errors at high bias. The XZZX code is deformed on a sublattice,
so initialization, readout, and the data-derived final detectors use a per-qubit effective
basis rather than a single global basis; for the plain surface code this reduces to a
uniform X basis.

## Syndrome Measurement

A single uniform gadget works for any X/Z/Y leg mix: the ancilla is prepared in $\ket{+}$,
a controlled-P (CX/CY/CZ, ancilla as control) is applied for each leg, and the ancilla is
measured in the X basis (phase kickback). The same code path therefore serves both the CSS
and the XZZX checks. Detectors compare each round to the previous one. Code-capacity and
phenomenological use a perfect reference round 0 to make the detectors deterministic; the
**circuit-level** model instead keeps round 0 noisy and closes both time boundaries with
explicit boundary detectors (initial + final) plus a reset and readout error, for a faithful
all-operations-noisy circuit (see [circuit_builder/README.md](./src/circuit_builder/README.md) §4).

For the **circuit-level** model the per-leg controlled gates are not applied ancilla-by-ancilla
but scheduled into **4 parallel time steps** (`build_cnot_schedule` in
[circuit_builder/shared.py](./src/circuit_builder/shared.py)). Each stabilizer's legs are
assigned to steps by their lattice offset so that no data qubit is engaged by two ancillas in
the same step, and so that the effective fault distance stays equal to $d$ (a distance-preserving
order in the spirit of Fowler 2012 / Tomita–Svore 2014; verified by exhaustive search with Stim's
`shortest_graphlike_error`). This time-step structure is what lets the idle noise be attached
per step to exactly the qubits that are resting — see
[circuit_builder/README.md](./src/circuit_builder/README.md) §4.

## Decoding and Logical Error Rate

For each configuration the built Stim circuit is converted to a detector error model
(`decompose_errors=True`) and decoded with PyMatching's minimum-weight perfect matching.
[simulation.py](./src/simulation.py) samples many shots, compares the decoded prediction
against the true observable, and reports the logical error rate

$$
p_L = \frac{\text{number of mismatched shots}}{\text{shots}}
$$

together with the Bernoulli standard deviation
$\sigma = \sqrt{p_L (1 - p_L) / \text{shots}}$.

## Threshold Estimation

The threshold is estimated using finite-size scaling (FSS) analysis of the logical error
rate data obtained from Monte Carlo simulations ([threshold.py](./src/threshold.py)).

We assume the standard scaling form

$$
p_L(p, d) = f\left((p - p_{\mathrm{th}}) \cdot d^{\frac{1}{\nu}}\right)
$$

where $p$ is the physical error rate, $d$ is the code distance, $p_{\mathrm{th}}$ is the
threshold, and $\nu$ is the critical exponent.

The scaling function $f(x)$ is approximated by a second-order polynomial,

$$
f(x) \approx a + b x + c x^2
$$

The parameters $(p_{\mathrm{th}}, \nu, a, b, c)$ are obtained by nonlinear least-squares
fitting (`scipy.optimize.curve_fit`), weighting each point by its standard deviation.

The reported threshold corresponds to the value of $p_{\mathrm{th}}$ that best fits the
scaling collapse over all simulated code distances and physical error rates.

## Running

From the `src/` directory:

```bash
python main.py [--outdir results] [--shots N] [--threshold] [--diagramonly]
```

- `--outdir` — output directory for plots and diagrams (default `results`).
- `--shots` — number of Monte Carlo shots per data point.
- `--threshold` — overlay the fitted CSS/XZZX thresholds as vertical lines on the plots.
- `--diagramonly` — skip the simulation and only regenerate the circuit diagrams.

It produces, in `--outdir`:

- `result_<eta>.png` — logical-vs-physical error rate curves (log-scaled y-axis) for both
  codes across the simulated distances, with optional threshold lines.
- `diagrams/<code>_detslice.svg` and `diagrams/<code>_timeline.svg` — Stim detector-slice
  and timeline visualizations of each code's circuit.

## Dependencies

- Python 3
- [stim](https://pypi.org/project/stim/)
- [pymatching](https://pypi.org/project/PyMatching/)
- numpy
- scipy
- matplotlib
- [rich](https://pypi.org/project/rich/) (progress bars)

```bash
pip install stim pymatching numpy scipy matplotlib rich
```
