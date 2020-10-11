include("ndid_core_supp.jl")

function eco_CR!(du, u, p, t)
    @unpack In1, In2, rn, a1, a11, a23, a2, a42, a43, h42, h43, b1, b11, b2, f, e1, e2, e3, dD, df, dc, d1, d11, d2 = p

    du[1] = In1 - rn * u[1] -a1 * u[1] * u[2] / ( b1 + u[1]) - a11 * u[1] * u[3] / ( b11 + u[1]) + dD * u[5]
    du[2] = f * a1 * u[1] * u[2] / (b1 + u[1]) - a42 * u[2] * u[4] / (1. + a42*h42*u[2]+a43*h43*u[3]) - (d1 + e1) * u[2]
    du[3] = f * a11 * u[1] * u[3] / (b11 + u[1]) -  a43 * u[3] * u[4] / (1. + a42*h42*u[2]+a43*h43*u[3]) - (d11 + e1) * u[3]
    du[4] = f * a42 * u[2] * u[4] / (1. + a42*h42*u[2]+a43*h43*u[3]) +  f * a43 * u[3] * u[4] / (1. + a42*h42*u[2]+a43*h43*u[3]) - (d2 + e2) * u[4] - df * u[4]
    du[5] = d1 * u[2] + d11 * u[3] + d2 * u[4] -  (dD + e3) * u[5]

    du[6] = In1 - rn * u[6] -a1 * u[6] * u[7] / ( b1 + u[6]) - a11 * u[6] * u[8] / ( b11 + u[6]) + dD * u[10]
    du[7] = f * a1 * u[6] * u[7] / (b1 + u[6]) - a42 * u[7] * u[9] / (1. + a42*h42*u[7]+a43*h43*u[8]) - (d1 + e1) * u[7]
    du[8] = f * a11 * u[6] * u[8] / (b11 + u[6]) - a43 * u[8] * u[9] / (1. + a42*h42*u[7]+a43*h43*u[8]) - (d11 + e1) * u[8]
    du[9] = f * a42 * u[7] * u[9] / (1. + a42*h42*u[7]+a43*h43*u[8]) + f * a43 * u[8] * u[9] / (1. + a42*h42*u[7]+a43*h43*u[8]) - (d2 + e2) * u[9] - df * u[9]
    du[10] = d1 * u[7] + d11 * u[8] + d2 * u[9] -  (dD + e3) * u[10]

    du[11] = In1 - rn * u[11] -a1 * u[11] * u[12] / ( b1 + u[11]) - a11 * u[11] * u[13] / ( b11 + u[11]) + dD * u[15]
    du[12] = f * a1 * u[11] * u[12] / (b1 + u[11]) - a42 * u[12] * u[14] / (1. + a42*h42*u[12]+a43*h43*u[13]) - (d1 + e1) * u[12]
    du[13] = f * a11 * u[11] * u[13] / (b11 + u[11])  - a43 * u[13] * u[14] / (1. + a42*h42*u[12]+a43*h43*u[13]) - (d11 + e1) * u[13]
    du[14] = f * a42 * u[12] * u[14] / (1. + a42*h42*u[12]+a43*h43*u[13]) + f * a43 * u[13] * u[14] / (1. + a42*h42*u[12]+a43*h43*u[13]) - (d2 + e2) * u[14] - df * u[14]
    du[15] = d1 * u[12] + d11 * u[13] + d2 * u[14] -  (dD + e3) * u[15]

# division by 3 to acount for 3fold increase in water volume
    du[16] = In1 - rn * u[16] -a1 * u[16] * u[17] / ( b1 + u[16]) - a11 * u[16] * u[18] / ( b11 + u[16]) + dD * u[20] 
    du[17] = f * a1 * u[16] * u[17] / (b1 + u[16]) - a42 * u[17] * u[19] / (1. + a42*h42*u[17]+a43*h43*u[18]) - (d1 + e1) * u[17]
    du[18] = f * a11 * u[16] * u[18] / (b11 + u[16]) - a43 * u[18] * u[19] / (1. + a42*h42*u[17]+a43*h43*u[18]) - (d11 + e1) * u[18]
    du[19] = f * a42 * u[17] * u[19] / (1. + a42*h42*u[17]+a43*h43*u[18]) + f * a43 * u[18] * u[19] / (1. + a42*h42*u[17]+a43*h43*u[18]) - (d2 + e2) * u[19] + df * (u[4] + u[9] + u[14])/3. - df * u[19]
    du[20] = d1 * u[17] + d11 * u[18] + d2 * u[19] -  (dD + e3) * u[20]

    du[21] = In1 - rn * u[21] -a1 * u[21] * u[22] / ( b1 + u[21]) - a11 * u[21] * u[23] / ( b11 + u[21]) + dD * u[25]
    du[22] = f * a1 * u[21] * u[22] / (b1 + u[21]) - a42 * u[22] * u[24] / (1. + a42*h42*u[22]+a43*h43*u[23]) - (d1 + e1) * u[22]
    du[23] = f * a11 * u[21] * u[23] / (b11 + u[21]) - a43 * u[23] * u[24] / (1. + a42*h42*u[22]+a43*h43*u[23]) - (d11 + e1) * u[23]
    du[24] = f * a42 * u[22] * u[24] / (1. + a42*h42*u[22]+a43*h43*u[23]) + f * a43 * u[23] * u[24] / (1. + a42*h42*u[22]+a43*h43*u[23]) - (d2 + e2) * u[24] + df*u[19]
    du[25] = d1 * u[22] + d11 * u[23] + d2 * u[24] -  (dD + e3) * u[25]

    return du
end


@with_kw mutable struct EcoPar
    In1 = 0.15
    In2 = 0.0005
    rn = 0.05
    f = 0.8
    e1 = 0.05
    e2 = 0.05
    e3 = 0.05
    a1 = 0.6
    a11 = 0.55
    a23 = 0.07
    a2 = 0.25
    h42=.0005
    h43=.005
    b1 = 0.040
    b11 = 0.043
    b2 = 0.040
    a42 = 0.350
    a43 = 0.05
    dD = 0.1
    dfd = 0.0
    dcd = 0.0
    d1 = 0.10
    d11 = 0.11
    df2 = 0.001
    df = 0.2
    dc = 0.01
    d2 = 0.001
end

# Setup a default set
param = EcoPar()


# Example of simulation
u0 = fill(0.3, 25)
tspan = (0.0, 10000.0)
prob = ODEProblem(eco_CR!, u0, tspan, param)
sol = solve(prob, KenCarp5(), reltol = 1e-16)
plot(sol.t, sol.u)

#C's
@show "nutrients"
@show sol[1, end]
@show sol[6, end]
@show sol[11, end]
@show sol[16, end]
@show sol[21, end]

#Rs
@show "resources"
@show sol[2, end]
@show sol[7, end]
@show sol[12, end]
@show sol[17, end]
@show sol[22, end]

#RInedibles
@show "resources inedible"
@show sol[3, end]
@show sol[8, end]
@show sol[13, end]
@show sol[18, end]
@show sol[23, end]
#Cs
@show "consumers"
@show sol[4, end]
@show sol[9, end]
@show sol[14, end]
@show sol[19, end]
@show sol[24, end]

#Ds
@show "detritus"
@show sol[5, end]
@show sol[10, end]
@show sol[15, end]
@show sol[20, end]
@show sol[25, end]

plot(sol.t, sol[23, :])


function bif_analysis(vals)
    tspan = (0.0, 10000.0)
    u0 = fill(0.03, 25)
    # inedibles non-zero

    pbif = EcoPar()

    apoints3 = []
    apoints23 = []
    apoints4 = []
    apoints5 = []
    apoints24 = []
    apoints25 = []
    apoints2 = []
    apoints17 = []
    apoints22 = []
    acv3 = []
    acv23 = []

    for df in vals
        pbif.df = df
        prob = ODEProblem(eco_CR!, u0, tspan, pbif)
        sol = solve(prob, KenCarp5(), reltol = 1e-16)
        #@show sol.t
        maxes3 = local_max(sol(5000:1:10000)[3, :])
        mins3 = local_min(sol(5000:1:10000)[3, :])

        maxes23 = local_max(sol(5000:1:10000)[23, :])
        mins23 = local_min(sol(5000:1:10000)[23, :])

        maxes2 = local_max(sol(5000:1:10000)[2, :])
        mins2 = local_min(sol(5000:1:10000)[2, :])

        maxes17 = local_max(sol(5000:1:10000)[17, :])
        mins17 = local_min(sol(5000:1:10000)[17, :])

        maxes22 = local_max(sol(5000:1:10000)[22, :])
        mins22 = local_min(sol(5000:1:10000)[22, :])

        maxes4 = local_max(sol(5000:1:10000)[4, :])
        mins4 = local_min(sol(5000:1:10000)[4, :])

        maxes5 = local_max(sol(5000:1:10000)[5, :])
        mins5 = local_min(sol(5000:1:10000)[5, :])

        maxes24 = local_max(sol(5000:1:10000)[24, :])
        mins24 = local_min(sol(5000:1:10000)[24, :])

        maxes25 = local_max(sol(5000:1:10000)[25, :])
        mins25 = local_min(sol(5000:1:10000)[25, :])

        SD23 = std(sol(5000:1:10000)[23, :])
        MN23 = mean(sol(5000:1:10000)[23, :])
        CV23 = SD23 / MN23

        SD3 = std(sol(5000:1:10000)[3, :])
        MN3 = mean(sol(5000:1:10000)[3, :])
        CV3 = SD3 / MN3

        for m in maxes3
            push!(apoints3, [df, m])
        end
        for l in mins3
            push!(apoints3, [df, l])
        end

        for m in maxes23
            push!(apoints23, [df, m])
        end
        for l in mins23
            push!(apoints23, [df, l])
        end

        for m in maxes22
            push!(apoints22, [df, m])
        end
        for l in mins22
            push!(apoints22, [df, l])
        end

        for m in maxes17
            push!(apoints17, [df, m])
        end
        for l in mins17
            push!(apoints17, [df, l])
        end

        for m in maxes2
            push!(apoints2, [df, m])
        end
        for l in mins2
            push!(apoints2, [df, l])
        end

        for m in maxes4
            push!(apoints4, [df, m])
        end
        for l in mins4
            push!(apoints4, [df, l])
        end

        for m in maxes24
            push!(apoints24, [df, m])
        end
        for l in mins24
            push!(apoints24, [df, l])
        end

        for m in maxes25
            push!(apoints25, [df, m])
        end
        for l in mins25
            push!(apoints25, [df, l])
        end

        for m in maxes5
        push!(apoints5, [df, m])
        end
        for l in mins5
        push!(apoints5, [df, l])
        end

        for m in CV23
            push!(acv23, [df, m])
        end
        for m in CV3
            push!(acv3, [df, m])
        end

        @show df
    end

    return VectorOfArray(apoints2),VectorOfArray(apoints17),VectorOfArray(apoints22)
end


# Fig.S.5.B.ii (model 2)

bif = bif_analysis(0.0:0.001:.10)

let
    figure()

    subplot(131)
    plot(bif[1][1, :], bif[1][2, :], "ko", markersize = 1.4,color = "black")
    xlabel("df", fontname = "Times New Roman", fontsize = 12)
    ylabel("Max & Min R")
    xlim(0.0,.05)
    ylim(0.0,0.2)

    subplot(132)
    plot(bif[2][1, :], bif[2][2, :], "ko", markersize = 1.4)
    xlabel("df", fontname = "Times New Roman", fontsize = 12)
    ylabel("Max & Min R")
    xlim(0.0,.05)
    ylim(0.0,0.2)

    subplot(133)
    plot(bif[3][1, :], bif[3][2, :], "ko", markersize = 1.4)
    xlabel("df", fontname = "Times New Roman", fontsize = 12)
    ylabel("Max & Min R")
    xlim(0.0,.05)
    ylim(0.0,0.2)

    tight_layout()

end
