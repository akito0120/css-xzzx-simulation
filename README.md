# Comparison of Rotated Surface Codes and XZZX Codes under Biased Noise Models


## Noise Model and Bias Parameterization

The noise model is the code capacity model, where independent biased Pauli noise is applied only to data qubits, and syndrome measurements are assumed to be perfect. In the simulation, one noisy round is compared against an ideal reference round for detector construction.

Biased Pauli noise is parameterized by the bias

```
eta = p_Z / (p_X + p_Y), with p_X = p_Y, p = p_X + p_Y + p_Z
```

so that

```
p_X = p_Y = p / (2 (1 + eta))
p_Z = p * eta / (1 + eta)
```

eta = 0.5 is standard depolarizing noise; eta -> inf is pure dephasing (Z only). p is the total single-qubit error probability.

## Memory basis

The experiments run an X-basis memory (they track the logical X operator), which is sensitive to the dominant Z errors at high bias. The XZZX code is deformed on a sublattice, so initialization, readout, and the data-derived final detectors use a per-qubit effective basis rather than a single global basis; for the plain surface code this reduces to a uniform X basis.

## Syndrome Measurement

A single uniform gadget works for any X/Z/Y leg mix: the ancilla is prepared in $\ket{+}$, a controlled-P (CX/CY/CZ, ancilla as control) is applied for each leg, and the ancilla is measured in the X basis (phase kickback). The same code path therefore serves both the CSS and the XZZX checks. Detectors compare each round to the previous one, with a perfect reference round 0 so that every detector is deterministic regardless of check type.

