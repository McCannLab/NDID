include("ndid_core_supp.jl")

function eco_CR!(du, u, p, t)
    @unpack In1, In2, rn, a1, a11, a23, a2, a32, h32, b1, b11, b2, f, e1, e2, e3, dD, df, dc, d1, d11, d2 = p

    # I have started making sure things do not shoot through zeroes numerically, this is a numerical error but there
    # might be a better way of doing this and it is not working the way I want yet.

    du[1] = In1 - rn * u[1] -a1 * u[1] * u[2] / ( b1 + u[1]) + dD * u[4]
    du[2] = f * a1 * u[1] * u[2] / (b1 + u[1]) - a32 * u[2] * u[3] / (1. + a32 * h32 * u[2]) - (d1 + e1) * u[2]
    du[3] = f * a32 * u[2] * u[3] / (1. + a32 * h32 * u[2]) - (d2 + e2) * u[3] - dc * u[3]
    du[4] = d1 * u[2] + d2 * u[3] -  (dD + e3) * u[4]

    du[5] = In1 - rn * u[5] -a1 * u[5] * u[6] / ( b1 + u[5]) + dD * u[8]
    du[6] = f * a1 * u[5] * u[6] / (b1 + u[5]) - a32 * u[6] * u[7] / (1. + a32 * h32 * u[6]) - (d1 + e1) * u[6]
    du[7] = f * a32 * u[6] * u[7] / (1. + a32 * h32 * u[6]) - (d2 + e2) * u[7] - dc * u[7]
    du[8] = d1 * u[6] + d2 * u[7] -  (dD + e3) * u[8]

    du[9] = In1 - rn * u[9] -a1 * u[9] * u[10] / ( b1 + u[9]) + dD * u[12]
    du[10] = f * a1 * u[9] * u[10] / (b1 + u[9]) - a32 * u[10] * u[11] / (1. + a32 * h32 * u[10]) - (d1 + e1) * u[10]
    du[11] = f * a32 * u[10] * u[11] / (1. + a32 * h32 * u[10]) - (d2 + e2) * u[11] - dc * u[11]
    du[12] = d1 * u[10] + d2 * u[11] -  (dD + e3) * u[12]

   # added /3 to dilute consumers in stream 
    du[13] = In1 - rn * u[13] -a1 * u[13] * u[14] / ( b1 + u[13]) + dD * u[16]
    du[14] = f * a1 * u[13] * u[14] / (b1 + u[13]) - a32 * u[14] * u[15] / (1. + a32 * h32 * u[14]) - (d1 + e1) * u[14]
    du[15] = f * a32 * u[14] * u[15] / (1. + a32 * h32 * u[14]) - (d2 + e2) * u[15] + dc * (u[3] + u[7] + u[11])/3. - dc * u[15]
    du[16] = d1 * u[14] + d2 * u[15] -  (dD + e3) * u[16]

    du[17] = In1 - rn * u[17] -a1 * u[17] * u[18] / ( b1 + u[17]) + dD * u[20]
    du[18] = f * a1 * u[17] * u[18] / (b1 + u[17]) - a32 * u[18] * u[19] / (1. + a32 * h32 * u[18]) - (d1 + e1) * u[18]
    du[19] = f * a32 * u[18] * u[19] / (1. + a32 * h32 * u[18]) - (d2 + e2) * u[19] + dc * u[15]
    du[20] = d1 * u[18] + d2 * u[19] -  (dD + e3) * u[20]

    return du
end

@with_kw mutable struct EcoPar
    In1 = 0.25
    In2 = 0.0005
    rn = 0.05
    f = 0.8
    e1 = 0.05
    e2 = 0.05
    e3 = 0.05
    a1 = 0.6
    a11 = 0.4
    a23 = .07
    a2 = 0.25
    h32=.005
    b1 = 0.040
    b11 = 0.047
    b2 = 0.040
    a32=.20
    dD = 0.1
    dfd = 0.0
    dcd = 0.0
    d1 = 0.10
    d11 = 0.11
    noise = 0.0
    df2 = 0.00
    df = 0.0
# dc not working yet
    dc = 0.06
    d2 = 0.001
end


# Setup a default set
param = EcoPar()

# Example of simulation 
u0 = fill(0.3, 20)
tspan = (0.0, 10000.0)
prob = ODEProblem(eco_CR!, u0, tspan, param)
sol = solve(prob, KenCarp5(), reltol = 1e-16)
plot(sol.t, sol.u)

#C's
@show "nutrients"
@show sol[1, end]
@show sol[5, end]
@show sol[9, end]
@show sol[13, end]
@show sol[17, end]

#Rs
@show "resources"
@show sol[2, end]
@show sol[6, end]
@show sol[10, end]
@show sol[14, end]
@show sol[18, end]

#Cs
@show "consumers"
@show sol[3, end]
@show sol[7, end]
@show sol[11, end]
@show sol[15, end]
@show sol[19, end]

#Ds
@show "detritus"
@show sol[4, end]
@show sol[8, end]
@show sol[12, end]
@show sol[16, end]
@show sol[20, end]

plot(sol.t, sol[17, :])



function bif_analysis(vals)
    tspan = (0.0, 10000.0)
    u0 = fill(0.03, 20)

    # make a local set of parameters
    pbif = EcoPar()

    apoints3 = []
    apoints14 = []
    apoints5 = []
    apoints2 = []
    apoints18 = []

    acv3 = []
    acv23 = []

    for dc in vals
        pbif.dc = dc
        prob = ODEProblem(eco_CR!, u0, tspan, pbif)
        sol = solve(prob, KenCarp5(), reltol = 1e-16)
        #@show sol.t
        maxes3 = local_max(sol(5000:1:10000)[3, :])
        mins3 = local_min(sol(5000:1:10000)[3, :])

        maxes14 = local_max(sol(5000:1:10000)[14, :])
        mins14 = local_min(sol(5000:1:10000)[14, :])

        maxes2 = local_max(sol(5000:1:10000)[2, :])
        mins2 = local_min(sol(5000:1:10000)[2, :])

        maxes18 = local_max(sol(5000:1:10000)[18, :])
        mins18 = local_min(sol(5000:1:10000)[18, :])

        maxes4 = local_max(sol(5000:1:10000)[4, :])
        mins4 = local_min(sol(5000:1:10000)[4, :])

        maxes5 = local_max(sol(5000:1:10000)[5, :])
        mins5 = local_min(sol(5000:1:10000)[5, :])

        for m in maxes3
            push!(apoints3, [dc, m])
        end
        for l in mins3
            push!(apoints3, [dc, l])
        end

        for m in maxes14
            push!(apoints14, [dc, m])
        end
        for l in mins14
            push!(apoints14, [dc, l])
        end

        for m in maxes18
            push!(apoints18, [dc, m])
        end
        for l in mins18
            push!(apoints18, [dc, l])
        end

        for m in maxes2
            push!(apoints2, [dc, m])
        end
        for l in mins2
            push!(apoints2, [dc, l])
        end

        for m in maxes5
        push!(apoints5, [dc, m])
        end
        for l in mins5
        push!(apoints5, [dc, l])
        end

        @show dc
    end

    return VectorOfArray(apoints2),VectorOfArray(apoints14),VectorOfArray(apoints18)
end


# Fig.S5 B i 

bif = bif_analysis(0.0:0.001:.10)


let
    figure()

    subplot(131)
    plot(bif[1][1, :], bif[1][2, :], "ko",markersize = 1.4,color = "black")
    xlabel("df", fontname = "Times New Roman", fontsize = 12)
    ylabel("Max & Min R")
    xlim(0.0,0.10)
    ylim(0.0,0.8)

    subplot(132)
    plot(bif[2][1, :], bif[2][2, :], "ko",markersize = 1.4)
    xlabel("df", fontname = "Times New Roman", fontsize = 12)
    ylabel("Max & Min R")
    xlim(0.0,0.10)
    ylim(0.0,0.8)

    subplot(133)
    plot(bif[3][1, :], bif[3][2, :], "ko",markersize = 1.4)
    xlabel("df", fontname = "Times New Roman", fontsize = 12)
    ylabel("Max & Min R")
    xlim(0.0,0.1)
    ylim(0.0,0.8)

    tight_layout()

end
