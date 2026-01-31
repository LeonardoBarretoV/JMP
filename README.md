# Debt-Financed Fiscal Stimulus, Heterogeneity, and Welfare

This repository contains the Julia code used to compute steady states, transition dynamics, and welfare effects in an open-economy heterogeneous-agent model with debt-financed fiscal stimulus.

Author: **Leonardo Barreto**

---

## Repository Structure

- `src/`  
  Core model definitions, calibration routines, household problem, equilibrium solvers, and transition algorithms.

- `scripts/`  
  Main execution file: `main_JMP.jl`.

- `mckay_process/`, `kmv_process/`  
  Income process calibration data.

- `matlab/`  
  MATLAB scripts used to generate figures from exported `.mat` files.

- `output/`  
  Generated results (figures, tables, `.mat`, `.jld2` files).  
  Contents are not version-controlled.

---

## Requirements

- [Julia](https://julialang.org/downloads/)  
- Tested on **Julia 1.10.4**

The code uses multi-threading (via `Polyester` and base Julia threading) for faster computation.

To let Julia automatically determine the number of threads:

```bash
julia -t auto
```
See [Julia multi-threading documentation](https://docs.julialang.org/en/v1/manual/multi-threading/#man-multithreading) for details.

---

## Setup

From the root of the repository, install all required packages:

```bash
julia --project -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"
```
This installs all dependencies specified in Project.toml.

---

## Running the Code

From the root directory of the repository, execute:

```bash
julia --project -t auto scripts/main_JMP.jl
```

The script performs the following steps:

1. Constructs and calibrates the model
2. Computes the initial steady state
3. Computes the final steady state
4. Solves the transition path
5. Generates summary tables
6. Exports time paths to .mat files
7. Saves model objects in output/

All generated files are written to the output/ directory.

---

## Expected Runtime

Expected runtime (after package installation):

11 minutes

Tested on:
- Windows 11 Pro
- Intel(R) Core(TM) i9-9880H CPU @ 2.30GHz
- 8 cores (16 logical processors)
- 32 GB RAM
