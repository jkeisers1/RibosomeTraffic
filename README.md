# RibosomeTraffic.jl 🧬

**A high-performance stochastic simulation engine for ribosome traffic on mRNA.**

`RibosomeTraffic.jl` implements a totally asymmetric simple exclusion process (TASEP) with extended particles to model protein translation dynamics. It uses a custom **Gillespie Direct Method** solver optimized for local dependency updates, making it significantly faster than generic solvers for 1D lattice problems.

## 🚀 Features

* **Biologically Realistic:** Models ribosomes as large particles (footprint $\ell = 10$ codons) with steric exclusion.
* **Stochastic Dynamics:** Fully rigorous continuous-time Monte Carlo (Gillespie algorithm).
* **Detailed Kinetics:**
    * Initiation ($\alpha$) and Termination ($\beta$).
    * Codon-specific elongation rates ($\omega_i$).
    * Ribosome pausing and rescue (switching between active and paused states).
* **High Performance:** Custom $O(1)$ rate update logic avoids scanning the lattice every step.

## 📦 Installation

To use this package, clone the repository and activate it in Julia:

```julia
using Pkg
Pkg.activate(".") # Activate the local environment
Pkg.instantiate() # Install dependencies (StaticArrays, Distributions, etc.)