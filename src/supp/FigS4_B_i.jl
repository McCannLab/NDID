# Notes: if C capable of hitting Rined pretty good can do two things it seems, allow Hopf as one might expect -- parallel strong channels or
# two) allow a feddback up the edible to suppress the inedible, exclusion of the iinedible in effect
# SEEMS TO FLIP BETWEEN Red ALWAYS WINS OR Rined ALWAYS WINS WHICH IS ODD -- this occurs when the trde-off defied
# SEEMS a numerical problem with Red overshooting to extinction followed by C even at df=0.0
include("ndid_core_supp.jl")

function eco_CR!(du, u, p, t)
    @unpack In1, In2, rn, a1, a11, a23, a2, a42, a43, h42, h43, b1, b11, b2, f, e1, e2, e3, dD, df, dc, d1, d11, d2, z1, z2 = p

    du[1] = In1 - rn * u[1] -a1 * u[1] * u[2] / ( b1 + u[1]) - a11 * u[1] * u[3] / ( b11 + u[1]) + dD * u[5] - df * u[1]
    du[2] = f * a1 * u[1] * u[2] / (b1 + u[1]) - a42 * u[2]/u[4]^z2 * u[4] / (1. + a42*h42*u[2]/u[4]^z2+a43*h43*u[3]/u[4]^z2) - (d1 + e1) * u[2]
    du[3] = f * a11 * u[1] * u[3] / (b11 + u[1]) -  a43 * u[3]/u[4]^z2 * u[4] / (1. + a42*h42*u[2]/u[4]^z2+a43*h43*u[3]/u[4]^z2) - (d11 + e1) * u[3]
    du[4] = f * a42 * u[2]/u[4]^z2 * u[4] / (1. + a42*h42*u[2]/u[4]^z2+a43*h43*u[3]/u[4]^z2) +  f * a43 * u[3]/u[4]^z2 * u[4] / (1. + a42*h42*u[2]/u[4]^z2+a43*h43*u[3]/u[4]^z2) - (d2 + e2) * u[4]^z1
    du[5] = d1 * u[2] + d11 * u[3] + d2 * u[4]^z1 -  (dD + e3) * u[5]

    du[6] = In1 - rn * u[6] -a1 * u[6] * u[7] / ( b1 + u[6]) - a11 * u[6] * u[8] / ( b11 + u[6]) + dD * u[10] - df * u[6]
    du[7] = f * a1 * u[6] * u[7] / (b1 + u[6]) - a42 * u[7]/u[9]^z2 * u[9] / (1. + a42*h42*u[7]/u[9]^z2+a43*h43*u[8]/u[9]^z2) - (d1 + e1) * u[7]
    du[8] = f * a11 * u[6] * u[8] / (b11 + u[6]) - a43 * u[8]/u[9]^z2 * u[9] / (1. + a42*h42*u[7]/u[9]^z2+a43*h43*u[8]/u[9]^z2) - (d11 + e1) * u[8]
    du[9] = f * a42 * u[7]/u[9]^z2 * u[9] / (1. + a42*h42*u[7]/u[9]^z2+a43*h43*u[8]/u[9]^z2) + f * a43 * u[8]/u[9]^z2 * u[9] / (1. + a42*h42*u[7]/u[9]^z2+a43*h43*u[8]/u[9]^z2) - (d2 + e2) * u[9]^z1
    du[10] = d1 * u[7] + d11 * u[8] + d2 * u[9]^z1 -  (dD + e3) * u[10]

    du[11] = In1 - rn * u[11] -a1 * u[11] * u[12] / ( b1 + u[11]) - a11 * u[11] * u[13] / ( b11 + u[11]) + dD * u[15] - df * u[11]
    du[12] = f * a1 * u[11] * u[12] / (b1 + u[11]) - a42 * u[12]/u[14]^z2 * u[14] / (1. + a42*h42*u[12]/u[14]^z2+a43*h43*u[13]/u[14]^z2) - (d1 + e1) * u[12]
    du[13] = f * a11 * u[11] * u[13] / (b11 + u[11])  - a43 * u[13]/u[14]^z2 * u[14] / (1. + a42*h42*u[12]/u[14]^z2+a43*h43*u[13]/u[14]^z2) - (d11 + e1) * u[13]
    du[14] = f * a42 * u[12]/u[14]^z2 * u[14] / (1. + a42*h42*u[12]/u[14]^z2+a43*h43*u[13]/u[14]^z2) + f * a43 * u[13]/u[14]^z2 * u[14] / (1. + a42*h42*u[12]/u[14]^z2+a43*h43*u[13]/u[14]^z2) - (d2 + e2) * u[14]^z1
    du[15] = d1 * u[12] + d11 * u[13] + d2 * u[14]^z1 -  (dD + e3) * u[15]

# 3 fold increase in water so divide by 3 to keep density
    du[16] = In1 - rn * u[16] -a1 * u[16] * u[17] / ( b1 + u[16]) - a11 * u[16] * u[18] / ( b11 + u[16]) + dD * u[20] + df * (u[1] + u[6] + u[11])/3. - df * u[16]
    du[17] = f * a1 * u[16] * u[17] / (b1 + u[16]) - a42 * u[17]/u[19]^z2 * u[19] / (1. + a42*h42*u[17]/u[19]^z2+a43*h43*u[18]/u[19]^z2) - (d1 + e1) * u[17]
    du[18] = f * a11 * u[16] * u[18] / (b11 + u[16]) - a43 * u[18]/u[19]^z2 * u[19] / (1. + a42*h42*u[17]/u[19]^z2+a43*h43*u[18]/u[19]^z2) - (d11 + e1) * u[18]
    du[19] = f * a42 * u[17]/u[19]^z2 * u[19] / (1. + a42*h42*u[17]/u[19]^z2+a43*h43*u[18]/u[19]^z2) + f * a43 * u[18]/u[19]^z2 * u[19] / (1. + a42*h42*u[17]/u[19]^z2+a43*h43*u[18]/u[19]^z2) - (d2 + e2) * u[19]^z1
    du[20] = d1 * u[17] + d11 * u[18] + d2 * u[19]^z1 -  (dD + e3) * u[20]

    du[21] = In1 - rn * u[21] -a1 * u[21] * u[22] / ( b1 + u[21]) - a11 * u[21] * u[23] / ( b11 + u[21]) + dD * u[25] + df * u[16]
    du[22] = f * a1 * u[21] * u[22] / (b1 + u[21]) - a42 * u[22]/u[24]^z2 * u[24] / (1. + a42*h42*u[22]/u[24]^z2+a43*h43*u[23]/u[24]^z2) - (d1 + e1) * u[22]
    du[23] = f * a11 * u[21] * u[23] / (b11 + u[21]) - a43 * u[23]/u[24]^z2 * u[24] / (1. + a42*h42*u[22]/u[24]^z2+a43*h43*u[23]/u[24]^z2) - (d11 + e1) * u[23]
    du[24] = f * a42 * u[22]/u[24]^z2 * u[24] / (1. + a42*h42*u[22]/u[24]^z2+a43*h43*u[23]/u[24]^z2) + f * a43 * u[23]/u[24]^z2 * u[24] / (1. + a42*h42*u[22]/u[24]^z2+a43*h43*u[23]/u[24]^z2) - (d2 + e2) * u[24]^z1
    du[25] = d1 * u[22] + d11 * u[23] + d2 * u[24]^z1 -  (dD + e3) * u[25]

    return du
end



@with_kw mutable struct EcoPar
    In1 = 0.075
    In2 = 0.005
    rn = 0.05
    f = 0.8
    e1 = 0.05
    e2 = 0.005
    e3 = 0.05
    a1 = 0.6
    a11 = 0.42
    a23 = .07
    a2 = 0.25
    h42=.005
    h43=.005
    b1 = 0.040
    b11 = 0.05
    b2 = 0.040
    a42=.15
    a43=.03
    dD = 0.1
    dfd = 0.0
    dcd = 0.0
    d1 = 0.10
    d11 = 0.11
    noise = 0.0
    df2 = 0.001
    df = .10
    dc = 0.0
    d2 = 0.02
    z1=1.0
    z2=.0
end
# Setup a default set
param = EcoPar()

u0 = fill(0.1, 25)
u0[3]=0.01
u0[8]=0.01
u0[13]=0.01
u0[18]=0.01
u0[23]=0.01

tspan = (0.0, 10000.0)

prob = ODEProblem(eco_CR!, u0, tspan, param)
sol = solve(prob, KenCarp5(), reltol = 1e-16)
plot(sol.t, sol.u)

# Note here that self-regulation is a powerful stabilzing force as is well known in this simple model
# BUT it nonetheless drives the proliferation of indedibles and detritus which in lakes can drive
# destabilziation through another mechanism not in this odel, loss of oxygen due to proliferation of bacteria
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

plot(sol.t, sol[21, :])


function bif_analysis(vals)
    tspan = (0.0, 10000.0)
    u0 = fill(0.05, 25)
    u0[3]=0.01
    u0[8]=0.01
    u0[13]=0.01
    u0[18]=0.01
    u0[23]=0.01

    # make a local set of parameters
    pbif = EcoPar()

    apoints3 = []
    apoints23 = []
    apoints4 = []
    apoints5 = []
    apoints24 = []
    apoints25 = []
    apoints2 = []
    apoints17 = []
    apoints18 = []
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

        maxes18 = local_max(sol(5000:1:10000)[18, :])
        mins18 = local_min(sol(5000:1:10000)[18, :])

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

                for m in maxes18
                    push!(apoints18, [df, m])
                end
                for l in mins18
                    push!(apoints18, [df, l])
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


# Fig.S4B.i
bif = bif_analysis(0.0:0.01:1.0)

let
    figure()
    subplot(131)
    plot(bif[1][1, :], bif[1][2, :], "ko",markersize = 1.4,color = "black")
    xlabel("df", fontname = "Times New Roman", fontsize = 12)
    ylabel("Max & Min R")
    xlim(0.0,1.0)
    ylim(0.0,1.0)

    subplot(132)
    plot(bif[2][1, :], bif[2][2, :], "ko",markersize = 1.4)
    xlabel("df", fontname = "Times New Roman", fontsize = 12)
    ylabel("Max & Min R")
    xlim(0.0,1.0)
    ylim(0.0,1.0)

    subplot(133)
    plot(bif[3][1, :], bif[3][2, :], "ko",markersize = 1.4)

    xlabel("df", fontname = "Times New Roman", fontsize = 12)
    ylabel("Max & Min R")
    xlim(0.0,1.0)
    ylim(0.0,1.0)

    tight_layout()

end
