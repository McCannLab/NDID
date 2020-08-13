# Example of simulation
include("ndid_core.jl")

# System NCDRe
## set paramaters
param = EcoPar()
# NB this is valid for the two kind of ecosystem considered
## and  time span
tspan = (0.0, 30000.0)

# System NCDRe
## Nutrient (N), Detritus (D), Consumer (C), Re (Ressource edible)
# initial conditions
u01 = fill(0.03, 20)
# set the problem NDCReRle
prob1 = ODEProblem(eco_NDCRe!, u01, tspan, param)
# # slow hub node would be
# prob1 = ODEProblem(eco_NDCRe_sh!, u0, tspan, param)
## solve it
sol1 = solve(prob1, Rodas4P(), reltol = 1e-16)
# See results
show_results(sol1, 28000:30000, 4)
# Plot resutls
plot(sol1.t, sol1.u)

# System NDReRle
## Nutrient (N), Detritus (D), Consumer (C), Re (Ressource edible) and
## Rle (Ressource less edible)
param = EcoPar()
tspan = (0.0, 30000.0)
u02 = fill(0.03, 25)
# set the problem NDCReRle
prob2 = ODEProblem(eco_NDCReRle!, u02, tspan, param)
sol2 = solve(prob2, Rodas4P(), reltol = 1e-32)
# prob2 = ODEProblem(eco_NDCReRle_sh!, u0, tspan, param)
show_results(sol2, 28000:30000)
plot(sol2.t, sol2.u)

