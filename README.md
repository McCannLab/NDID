# NDID
Nutrient-Driven Instability at a Distance -- Research compendium

This repository includes all the code we used in "" (DOI:[Ecosystem Entanglement and the Propagation of Nutrient-Driven Instability](https://doi.org/10.1101/2020.04.20.050302)).

## Installation

This requires Julia (>=v.1.4.0) and the following packages:

- LinearAlgebra
- Statistics
- RecursiveArrayTools
- Parameters
- DifferentialEquations
- PyPlot


## Content

Briefly, in this repository:

- `src/ndid_core.jl`: includes all core functions, i.e. functions that encode the systems of differential equations, the function used for the bifurcation analysis, the figures
- `src/simulation_example.jl`: Example of a simulation
- `src/simulations_main.jl`: bifurcation analysis as presented in the main text
- `equilibrium_results.jl`:equilibrium computation