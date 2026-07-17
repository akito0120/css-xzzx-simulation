# CSS vs. XZZX Surface Code Simulation

## Overview

This project compares the error-correction performance of the rotated CSS surface code and the XZZX surface code under Z-biased noise, where phase-flip errors dominate. The approach is a Monte Carlo stabilizer simulation built on Stim and Sinter. For each code, the simulator runs a memory experiment under circuit-level biased noise and measures the logical error rate while sweeping the noise bias, the code distance, and the physical error rate. Thresholds are then extracted from the resulting curves by finite-size-scaling analysis.

## Project Structure

```
css_xzzx_simulation/
├── src/
│   ├── main.py              # Entry point: runs the sweep and renders outputs
│   ├── config.py            # Swept parameters (bias, code types, distances, error-rate windows)
│   ├── code_builder.py      # Builds the rotated CSS and XZZX codes
│   ├── simulation.py        # Monte Carlo sweep, decoder selection, error-rate statistics
│   ├── threshold.py         # Finite-size-scaling threshold estimation
│   ├── uf_decoder.py        # Union-Find decoder wrapper
│   ├── visualization.py     # Figure and circuit-diagram rendering
│   └── circuit_builder/     # Noise-model circuit construction (code-capacity / phenomenological / circuit-level)
├── full_results/            # Simulation outputs: samples.csv, figures/, diagrams/
└── requirements.txt         # Pinned Python dependencies
```
