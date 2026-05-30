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

## Memory Basis

The experiments run an X-basis memory (they track the logical X operator), which is sensitive to the dominant Z errors at high bias. The XZZX code is deformed on a sublattice, so initialization, readout, and the data-derived final detectors use a per-qubit effective basis rather than a single global basis; for the plain surface code this reduces to a uniform X basis.

## Syndrome Measurement

A single uniform gadget works for any X/Z/Y leg mix: the ancilla is prepared in $\ket{+}$, a controlled-P (CX/CY/CZ, ancilla as control) is applied for each leg, and the ancilla is measured in the X basis (phase kickback). The same code path therefore serves both the CSS and the XZZX checks. Detectors compare each round to the previous one, with a perfect reference round 0 so that every detector is deterministic regardless of check type.

## Threshold Estimation

The threshold is estimated using finite-size scaling (FSS) analysis of the logical error rate data obtained from Monte Carlo simulations.

We assume the standard scaling form

$$
p_L = f\!\left((p - p_{\mathrm{th}}) \cdot d^{\frac{1}{\nu}}\right)
$$

where $p$ is the physical error rate, $d$ is the code distance, $p_{\mathrm{th}}$ is the threshold, and $\nu$ is the critical exponent.

The scaling function $f(x)$ is approximated by a second-order polynomial,

$$
f(x) \approx a + b x + c x^2
$$

The parameters $(p_{\mathrm{th}}, \nu, a, b, c)$ are obtained by nonlinear least-squares fitting.

The reported threshold corresponds to the value of $p_{\mathrm{th}}$ that minimizes the least-squares error over all simulated code distances and physical error rates.

