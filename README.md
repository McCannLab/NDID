# NDID: Nutrient-Driven Instability at a Distance -- Research compendium
[![Test](https://github.com/McCannLab/NDID/workflows/Test/badge.svg)](https://github.com/McCannLab/NDID/actions)

This repository includes all the code we used in *Landscape Modification and Nutrient-Driven Instability at a Distance* (DOI:XXX, see also the preprint [DOI: 10.1101/2020.04.20.050302 *](https://doi.org/10.1101/2020.04.20.050302)).


## Installation

We used [Julia](https://julialang.org/) v.1.5.0 and the following packages:

- [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
- [Statistics](https://docs.julialang.org/en/v1/stdlib/Statistics/)
- [NLsolve](https://github.com/JuliaNLSolvers/NLsolve.jl)
- [RecursiveArrayTools](https://github.com/SciML/RecursiveArrayTools.jl)
- [Parameters](https://github.com/mauro3/Parameters.jl)
- [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl)
- [PyPlot](https://github.com/JuliaPy/PyPlot.jl)

See `src/install_packages.jl` as well as the [Action](https://github.com/McCannLab/NDID/actions) tab for more details.


## Content

Briefly:

- `src/ndid_core.jl`: includes all core functions:
    - functions that encode the systems of differential equations,
    - bifurcation analysis
    - figures
- `src/simulation_example.jl`: includes one example of a simulation (i.e., solving the differential equations for a given time span) for both systems,
- `src/simulations_main.jl`: includes the bifurcation analyses as presented in the main text,
- `equilibrium_results.jl`: uses the NLsolve packages to compute the equilibrium for both systems studied.