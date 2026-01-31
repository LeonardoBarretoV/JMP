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
