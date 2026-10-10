# anisotropic_pf_dcla_frc

A user-defined element (UEL) implementation in Abaqus for simulating fatigue damage in fiber-reinforced composites (FRCs) using an anisotropic phase-field model and a new degradation-informed constant-load-accumulation (D-CLA) scheme.

---

## 🔬 Overview

This repository contains the Fortran-based Abaqus UEL subroutine developed to model progressive damage in composite laminates under quasi-static and fatigue loading conditions. The implementation combines:

- An **anisotropic phase-field formulation** with modified driving force for fracture
- A **new degradation-informed constant-load-accumulation (D-CLA)** strategy for efficient fatigue simulation
- Coupling of damage and displacement fields in a monolithic framework using BFGS algorithm 

> The model is detailed in our article:  
> **Degradation-informed load accumulation framework for phase-field modeling of high-cycle fatigue**  
> _Composites Part A: Applied Science and Manufacturing_ (2025)  
> [doi:10.1016/j.compositesa.2025.109377](https://doi.org/10.1016/j.compositesa.2025.109377)

---

## 📁 Contents

- `src/anisotropic_pf_dcla_frc.for`: Fortran source code for the UEL
- `examples/anisotropic_pf_dcla_frc.inp`: Abaqus input file, single-edge notched tension (constant stiffness)
- `examples/Amplitude.txt`: load amplitude table used by the example
- `README.md`: Project documentation

Planned additions: single-edge notched tension (variable stiffness), double cantilever beam (DCB), open-hole tension (OHT), and Python post-processing scripts.

---

## 🛠️ Requirements

- **Abaqus 2020+** with UEL support
- **Fortran compiler** (e.g., Intel Fortran Compiler or Intel OneAPI Compiler)
- Optional: Python 3.x with `numpy`, `matplotlib` for postprocessing

---

## 🚀 Getting Started

1. Compile the UEL with your Abaqus installation:
   ```bash
   abaqus make library=src/anisotropic_pf_dcla_frc.for
2. Run an example job from the command line:
   ```bash
   abaqus job=<job_name> user=src/anisotropic_pf_dcla_frc.for input=examples/<example_file>.inp interactive
