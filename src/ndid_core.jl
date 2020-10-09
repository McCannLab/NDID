# -------------- Core functions --------------

# Packages required
using LinearAlgebra, Parameters, DifferentialEquations, RecursiveArrayTools
using Statistics: mean, std
using PyPlot


# Reference set of parameters
@with_kw mutable struct EcoPar
    In1 = 0.3
    In2 = 0.25  # used in slow node
    rn  = 0.05
    f   = 0.8
    e1  = 0.05
    e2  = 0.05
    e3  = 0.05
    a1  = 0.6
    a11 = 0.4
    a23 = 0.07
    a2  = 0.25
    h42 = 0.005   # => treated as h32 in NCDRe
    h43 = 0.005
    b1  = 0.040
    b11 = 0.047
    b2  = 0.040
    a42 = 0.250   # => treated as a32 in NCDRe
    a43 = 0.09
    dD  = 0.1
    d1  = 0.10
    d11 = 0.11
    df2 = 0.001   # used in slow node
    df  = 0.05
    d2  = 0.001
    hub_slow = false
end


# I- ODE - Edible resource

## interaction + diffusion 
function eco_NDCRe!(du, u, p, t)
    @unpack In1, df = p
    #
    j = 1
    for i in 1:5
        ind = j:(j + 3)
        du[ind] = eco_NDCRe_unit!(du[ind], u[ind], p, t)
        j += 4
    end
    # diffusion and input
    du[1] += In1 - df * u[1]
    du[5] += In1 - df * u[5]
    du[9] += In1 - df * u[9]
    du[13] += In1 + df * (u[1] + u[5] + u[9])/3. - df * u[13]
    du[17] += In1 + df * u[13]

    return du
end

## interaction + diffusion - slow hub node (sh)
function eco_NDCRe_sh!(du, u, p, t)
    @unpack In1, In2, df, df2 = p
    #
    j = 1
    for i in 1:5
        ind = j:(j + 3)
        du[ind] = eco_NDCRe_unit!(du[ind], u[ind], p, t)
        j += 4
    end
    # diffusion and input
    du[1] += In1 - df * u[1]
    du[5] += In1 - df * u[5]
    du[9] += In1 - df * u[9]
    du[13] += In2 + df * (u[1] + u[5] + u[9])/3. - df2 * u[13]
    du[17] += In1 + df2 * u[13]

    return du
end

## core ODE 
function eco_NDCRe_unit!(du, u, p, t)
    # u[1]: Phosphorus; u[2]: Resource; u[3]: Consumer; u[4]: Detritus
    @unpack rn, a1, a23, a2, a42, h42, b1, f, e1, e2, e3, dD, d1, d2 = p

    tmp1 = a1 * u[1] * u[2] / (b1 + u[1])
    tmp24 = a42 * u[2] * u[3] / (1. + a42 * h42 * u[2])
    du[1] = - rn * u[1] - tmp1 + dD * u[4]
    du[2] = f * tmp1 - tmp24 - (d1 + e1) * u[2]
    du[3] = f * tmp24 - (d2 + e2) * u[3]
    du[4] = d1 * u[2] + d2 * u[3] - (dD + e3) * u[4] # + d11 * u[3]

    return du
end



# II- ODE - edible and less edible resource

## interaction + diffusion
function eco_NDCReRle!(du, u, p, t)
    @unpack In1, df= p;
    #
    j = 1
    for i in 1:5
        ind = j:(j + 4)
        du[ind] = eco_NDCReRle_unit!(du[ind], u[ind], p, t)
        j += 5
    end
    # diffusion and input
    du[1] += In1 - df * u[1]
    du[6] += In1 - df * u[6]
    du[11] += In1 - df * u[11]
    du[16] += In1 + df * (u[1] + u[6] + u[11])/3. - df * u[16]
    du[21] += In1 + df * u[16]

    return du
end

## interaction + diffusion - slow hub node (sh)
function eco_NDCReRle_sh!(du, u, p, t)
    @unpack In1, In2, df, df2 = p;
    #
    j = 1
    for i in 1:5
        ind = j:(j + 4)
        du[ind] = eco_NDCReRle_unit!(du[ind], u[ind], p, t)
        j += 5
    end
    # diffusion and input
    du[1] += In1 - df * u[1]
    du[6] += In1 - df * u[6]
    du[11] += In1 - df * u[11]
    du[16] += In2 + df * (u[1] + u[6] + u[11])/3. - df2 * u[16]
    du[21] += In1 + df2 * u[16]

    return du
end

## core ODE 
function eco_NDCReRle_unit!(du, u, p, t)
    # u[1]: Phosphorus; u[2]: Resource; u[3]: Resource less edible;
    # u[4]: Consumer; u[5]: Detritus
    @unpack rn, a1, a11, a23, a2, a42, a43, h42, h43, b1, b11, f, e1, e2, e3, dD, d1, d11, d2 = p

    tmp1 = a1 * u[1] * u[2] / (b1 + u[1])
    tmp11 = a11 * u[1] * u[3] / (b11 + u[1])
    K = 1. + a42 * h42 * u[2] + a43 * h43 * u[3]
    tmp24 = a42 * u[2] * u[4] / K
    tmp34 = a43 * u[3] * u[4] / K
    du[1] = -rn * u[1] - tmp1 - tmp11 + dD * u[5]
    du[2] = f * tmp1 - tmp24 - (d1 + e1) * u[2]
    du[3] = f * tmp11 - tmp34 - (d11 + e1) * u[3]
    du[4] = f * tmp24 + f * tmp34 - (d2 + e2) * u[4]
    du[5] = d1 * u[2] + d11 * u[3] + d2 * u[4] - (dD + e3) * u[5]

    return du
end




# III- Bifurcation analysis

# we use bif_analysis2 cause bif_analysis was used in a previous version 
# see SI where it is still used.

function bif_analysis2(vals, FUN, u0, tspan, ind, pbif = EcoPar())
    nsp = size(u0, 1)
    mat_max = zeros(Float64, nsp, size(vals, 1))
    mat_min = zeros(Float64, nsp, size(vals, 1))
    mat_mean = zeros(Float64, nsp, size(vals, 1))
    mat_cv  = zeros(Float64, nsp, size(vals, 1))

    for (i, df) in enumerate(vals)
        printstyled("rep #", i, "  df=", df, "\n", color = :blue)
        pbif.df = df
        prob = ODEProblem(FUN, u0, tspan, pbif)
        # sol = solve(prob, Trapezoid(), reltol = 1e-16)
        sol = solve(prob, Rodas4P(), reltol = 1e-16)
        mat_max[:, i] = minimum(sol(ind), dims = 2)#[local_max(sol(ind)[k, :])[1] for k in 1:nsp]
        mat_min[:, i] = maximum(sol(ind), dims = 2)#[local_min(sol(ind)[k, :])[1] for k in 1:nsp]
        tmp = mean(sol(ind), dims = 2)
        tmp2 = std(sol(ind), dims = 2)
        mat_mean[:, i] = mean(sol(ind), dims = 2)
        mat_cv[:, i] = tmp2 ./ tmp
    end

    return (vals, mat_min, mat_max, mat_mean, mat_cv)
end



# IV- Figures 

function fig_bif(bif, ind, file, xlab, ylab, ymin, ymax, col = "black", lw = 1.8)
    figure(figsize = (2.8, 9))

    subplot(311)
    plot(bif[1], bif[2][ind[1], :], color = col, linestyle = "-", linewidth = lw)
    plot(bif[1], bif[3][ind[1], :], color = col, linestyle = "-", linewidth = lw)
    xlabel(xlab, fontname = "Times New Roman", fontsize = 12)
    # ylabel(ylab)
    xlim(0., 1.0)
    ylim(ymin, ymax)

    subplot(312)
    plot(bif[1], bif[2][ind[2], :], color = col, linestyle = "-", linewidth = lw)
    plot(bif[1], bif[3][ind[2], :], color = col, linestyle = "-", linewidth = lw)
    xlabel(xlab, fontname = "Times New Roman", fontsize = 12)
    # ylabel(ylab)
    xlim(0., 1.0)
    ylim(ymin, ymax)

    subplot(313)
    plot(bif[1], bif[2][ind[3], :], color = col, linestyle = "-", linewidth = lw)
    plot(bif[1], bif[3][ind[3], :], color = col, linestyle = "-", linewidth = lw)
    xlabel(xlab, fontname = "Times New Roman", fontsize = 12)
    # ylabel(ylab)
    xlim(0., 1.0)
    ylim(ymin, ymax)

    tight_layout()
    savefig(file, dpi = 300)
    plt.close()

end

function barplotmean(bif, ind, file, ymin, ymax, col = "black")
    figure(figsize = (2.6, 3.8))
    bar(1:3, bif[4][ind, 81][1,:], color = col, align = "center")
    ylim(ymin, ymax)
    tight_layout()
    savefig(file, dpi = 300)
    plt.close()
end

function fig_bif2(bif, ind, ind2, col, file, xlab, ylab, ymin, ymax)
    figure(figsize = (2.8, 9))

    subplot(311)
    plot(bif[1], bif[4][ind[1], :], color = col[1], linestyle = "-")
    plot(bif[1], bif[4][ind2[1], :], color = col[2], linestyle = "-")
    xlim(0., 1.0)
    ylim(ymin, ymax)

    subplot(312)
    plot(bif[1], bif[4][ind[2], :], color = col[1], linestyle = "-")
    plot(bif[1], bif[4][ind2[2], :], color = col[2], linestyle = "-")
    xlim(0., 1.0)
    ylim(ymin, ymax)

    subplot(313)
    plot(bif[1], bif[4][ind[3], :], color = col[1], linestyle = "-")
    plot(bif[1], bif[4][ind2[3], :], color = col[2], linestyle = "-")
    xlabel(xlab, fontname = "Times New Roman", fontsize = 12)
    xlim(0., 1.0)
    ylim(ymin, ymax)

    tight_layout()
    savefig(file, dpi = 300)
    plt.close()

end


function fig_bif3(bif, ind, col, file, xlab, ylab, ymin, ymax)
    figure(figsize = (2.8, 9))

    subplot(311)
    for i in 1:3
        plot(bif[1], bif[2][ind[i][1], :], color = col[i], linestyle = "-")
        plot(bif[1], bif[3][ind[i][1], :], color = col[i], linestyle = "-")
    end
    xlim(0., 1.0)
    ylim(ymin, ymax)

    subplot(312)
    for i in 1:3
        plot(bif[1], bif[2][ind[i][2], :], color = col[i], linestyle = "-")
        plot(bif[1], bif[3][ind[i][2], :], color = col[i], linestyle = "-")
    end
    xlim(0., 1.0)
    ylim(ymin, ymax)

    subplot(313)
    for i in 1:3
        plot(bif[1], bif[2][ind[i][3], :], color = col[i], linestyle = "-")
        plot(bif[1], bif[3][ind[i][3], :], color = col[i], linestyle = "-")
    end
    xlim(0., 1.0)
    ylim(ymin, ymax)

    tight_layout()
    savefig(file, dpi = 300)
    plt.close()

end


# V- Print results

function show_results(sol, ind, nsp = 5, nnd = 5)

    nm = ["Nutrients" "Resources" "Resources inedible" "Consumers" "Detritus"]

    for i in 1:nsp
        printstyled(nm[i], "\n", color = :blue)
        [@show mean(sol(ind)[k, :]) for k in i:nsp:((nnd - 1)*nsp + i)]
    end

    return mean(sol(ind), dims = 2)
end