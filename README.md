# anisotropic_pf_cla_frc

A user-defined element (UEL) implementation in Abaqus for simulating fatigue damage in fiber-reinforced composites (FRCs) using an anisotropic phase-field model and a new constant-load-accumulation (CLA) scheme.

---

## 🔬 Overview

This repository contains the Fortran-based Abaqus UEL subroutine developed to model progressive damage in composite laminates under quasi-static and fatigue loading conditions. The implementation combines:

- An **anisotropic phase-field formulation** for fracture
- A **cycle-level averaging (CLA)** strategy for efficient fatigue simulation
- Coupling of damage and displacement fields in a monolithic framework

> The model is detailed in our article:  
> **[Full Title of Your Paper]**  
> _Submitted to Journal_  
> [DOI/link to be added upon publication]

---

## 📁 Contents

- `uel_aniso_pf_cla.for`: Fortran source code for the UEL
- `input_examples/`: Sample Abaqus `.inp` files for:
  - Single-edge notched tension (Constant Stiffness)
  - Single-edge notched tension (Variable Stiffness)
  - Double cantilever beam (DCB)
  - Open-hole tension (OHT)
- `postprocessing/`: Python scripts to extract and visualize results
- `README.md`: Project documentation

---

## 🛠️ Requirements

- **Abaqus 2020+** with UEL support
- **Fortran compiler** (e.g., Intel Fortran Compiler or Intel OneAPI Compiler)
- Optional: Python 3.x with `numpy`, `matplotlib` for postprocessing

---

## 🚀 Getting Started

1. Compile the UEL with your Abaqus installation:
   ```bash
   abaqus make library=uel_aniso_pf_cla.for
2. Run an example job from the command line:
   ```bash
   abaqus job=<job_name> user=uel_aniso_pf_cla.for input=input_examples/<example_file>.inp interactive
