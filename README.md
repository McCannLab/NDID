## NDID: Nutrient-Driven Instability at a Distance (research compendium)
[![Check](https://github.com/McCannLab/NDID/workflows/Check/badge.svg)](https://github.com/McCannLab/NDID/actions)
[![DOI](https://zenodo.org/badge/287106212.svg)](https://zenodo.org/badge/latestdoi/287106212)


This repository includes the code we used in *Landscape Modification and Nutrient- Driven Instability at a Distance*: 

- Ecology Letters paper (DOI: [10.1111/ele.13644](https://onlinelibrary.wiley.com/doi/abs/10.1111/ele.13644))
- Preprint [DOI: 10.1101/2020.04.20.050302 *](https://doi.org/10.1101/2020.04.20.050302)).


### Installation

We used [Julia](https://julialang.org/) v.1.5.1 with the following packages:

- [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
- [Statistics](https://docs.julialang.org/en/v1/stdlib/Statistics/)
- [NLsolve](https://github.com/JuliaNLSolvers/NLsolve.jl)
- [RecursiveArrayTools](https://github.com/SciML/RecursiveArrayTools.jl)
- [Parameters](https://github.com/mauro3/Parameters.jl)
- [DifferentialEquations](https://github.com/SciML/DifferentialEquations.jl)
- [PyPlot](https://github.com/JuliaPy/PyPlot.jl)

See `src/install_packages.jl` as well as the [Actions](https://github.com/McCannLab/NDID/actions) tab for more details.


### Content

Briefly, this repository includes:

- `src/ndid_core.jl`: includes all core functions:
    - functions that encode the systems of differential equations,
    - bifurcation analysis
    - figures
- `src/simulation_example.jl`: includes one example of a simulation (i.e., solving the differential equations for a given time span) for both systems,
- `src/simulations_main.jl`: includes the bifurcation analyses as presented in the main text,
- `src/equilibrium_results.jl`: uses the NLsolve packages to compute the equilibrium for both systems studied.
- `src/sup/`: includes all Julia files to reproduce supplementary figures. 


### Running the analysis 

For the main analysis, you can either run the code in `simulation_main.jl` line by line in an interactive session, or you can run it as a script: 

```sh
julia src/simulations_main.jl
```

See also [Actions](https://github.com/McCannLab/NDID/actions) tab. Note that all files in `scr/sup` can be run independently. 


## Troubleshooting

Please use the [issue feature](https://github.com/McCannLab/NDID/issues) if you encounter any difficulty while trying to reproduce the analysis. 