# -------------- Core functions for Supplementary Information --------------

# Packages required
using LinearAlgebra, Parameters, DifferentialEquations, RecursiveArrayTools
using Statistics: mean, std
using PyPlot


function local_max(ts)
    lmaxs = Float64[]
    for i in 2:(length(ts) - 3)
        if ts[i - 1] < ts[i] && ts[i] > ts[i + 1]
            push!(lmaxs, ts[i])
            # @show ts[i]
        end
    end

    if length(lmaxs) == 0
        push!(lmaxs, ts[length(ts) - 3])
    end

    return lmaxs
end

function local_min(ts)
    lmins = Float64[]
    for i in 2:(length(ts) - 3)
        if ts[i - 1] > ts[i] && ts[i] < ts[i + 1]
            push!(lmins, ts[i])
            #@show ts[i]
        end
    end
    if length(lmins) == 0
        push!(lmins, ts[length(ts) - 3])
    end

    return lmins
end
