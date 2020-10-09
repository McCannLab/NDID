# -------------- Solve equation at equilibrium  --------------

using NLsolve

# I- NCDRe

## Parameter values
In1 = .3
rn  = 0.05
a1  = 0.6
b1  = 0.04
dD  = 0.1
df  = 0.05
f   = 0.8
a32 = 0.25
h32 = 0.005
d1  = 0.1
d2  = 0.001
e1  = 0.05
e2  = 0.05
e3  = 0.05

function f!(F, u)
    F[1] = In1 - rn * u[1] - a1 * u[1] * u[2] / ( b1 + u[1]) + dD * u[4] - df * u[1]
    F[2] = f * a1 * u[1] * u[2] / (b1 + u[1]) - a32 * u[2] * u[3] / (1. + a32*h32*u[2]) - (d1 + e1) * u[2]
    F[3] = f * a32 * u[2] * u[3] / (1. + a32*h32*u[2]) - (d2 + e2) * u[3]
    F[4] = d1 * u[2] + d2 * u[3] - (dD + e3) * u[4]
end

initU2 = (d2 + e2)/(f * a32 - a32 * h32 * (d2 + e2))

# function, jacobian (autodiff) and initial guess
# equilibrium for initial node i
soli = nlsolve(f!, [1.1; initU2; .6; .4], autodiff = :forward)
# check against results in example

# Hub node
function fh!(Fh, uh)
    Fh[1] = In1 - rn * uh[1] -a1 * uh[1] * uh[2] / ( b1 + uh[1]) + dD*uh[4] + df * soli.zero[1] - df * uh[1]
    Fh[2] = f * a1 * uh[1] * uh[2] / (b1 + uh[1]) - a32 * uh[2] * uh[3] / (1. + a32 * h32 * uh[2]) - (d1 + e1) * uh[2]
    Fh[3] = f * a32 * uh[2] * uh[3] / (1. + a32 * h32 * uh[2]) - (d2 + e2) * uh[3]
    Fh[4] = d1 * uh[2] + d2 * uh[3] - (dD + e3) * uh[4]
end

# function, jacobian (autodiff) and initial guess
# equilibrium for initial node i
solh = nlsolve(fh!, [1.1; initU2; .6; .4], autodiff = :forward)


function ft!(Ft, ut)
    Ft[1] = In1 - rn * ut[1] -a1 * ut[1] * ut[2] / ( b1 + ut[1]) + dD*ut[4] + df * solh.zero[1]
    Ft[2] = f * a1 * ut[1] * ut[2] / (b1 + ut[1]) - a32 * ut[2] * ut[3] / (1. + a32 * h32 * ut[2]) - (d1 + e1) * ut[2]
    Ft[3] = f * a32 * ut[2] * ut[3] / (1. + a32 * h32 * ut[2]) - (d2 + e2) * ut[3]
    Ft[4] = d1 * ut[2] + d2 * ut[3] - (dD + e3) * ut[4]
end

# function, jacobian (autodiff) and initial guess
# equilibrium for initial node i
solt = nlsolve(ft!, [ 1.1; initU2; .6; .4], autodiff = :forward)





# II- NCDReRle

## Additional parameters (remaining are set as above)

a11 = 0.35
b11 = 0.05
a42 = 0.25
h42 = 0.005
h43 = 0.005
a43 = 0.09
d11 = 0.11

function f!(F, u)
    F[1] = In1 - rn * u[1] -a1 * u[1] * u[2] / ( b1 + u[1]) - a11 * u[1] * u[3] / ( b11 + u[1]) + dD * u[5] - df * u[1]
    F[2] = f * a1 * u[1]*u[2] / (b1 + u[1]) - a42 * u[2] * u[4] / (1. + a42*h42*u[2]+a43*h43*u[3]) - (d1 + e1)*u[2]
    F[3] = f * a11 * u[1]*u[3] / (b11 + u[1]) - a43 * u[3] * u[4] / (1. + a42*h42*u[2]+a43*h43*u[3]) - (d11 + e1)*u[3]
    F[4] = f * a42 * u[2] * u[4] / (1. + a42*h42*u[2]+a43*h43*u[3]) +  f * a43 * u[3] * u[4] / (1. + a42*h42*u[2]+a43*h43*u[3]) - (d2 + e2) * u[4]
    F[5] = d1 * u[2] + d11 * u[3] + d2 * u[4] -  (dD + e3) * u[5]
end

#a43=0
initU1 = b1 * (d1 + e1) / (f * a1 -  (d1 + e1))
initU2 = (d2 + e2)/(f * a42 -  a42 * h42 * (d2 + e2))

#i nvasion criteria for inedible with N-values
Nstar = b11 * (d11 + e1)/(f * a11 - (d11+e1))
    if Nstar < 0
        Nstar = 1000000.
    end

if initU1 > Nstar
    initU1 = Nstar
end

# function, jacobian (autodiff) and initial guess
# equilibrium for initial node i
soli = nlsolve(f!, [ initU1; .01; .01;.1;.3], autodiff = :forward,ftol=1.0e-16,iterations = 1000000)
# mcpsolve(f!, [-1.e-16,0.0,0.0,0.0,0.0], [100,100,100,100,100], [ initU1; .01; .01;.1;.3], autodiff = :forward, ftol = 1.0e-16) # returns 6.25915e8


# not hub node
function fh!(Fh, uh)
    Fh[1] = In1 - rn * uh[1] -a1 * uh[1] * uh[2] / ( b1 + uh[1]) - a11 * uh[1] * uh[3] / ( b11 + uh[1]) + dD * uh[5] + df * soli.zero[1] - df * uh[1]
    Fh[2] = f * a1 * uh[1] * uh[2] / (b1 + uh[1]) - a42 * uh[2] * uh[4] / (1. + a42*h42*uh[2]+a43*h43*uh[3]) - (d1 + e1) * uh[2]
    Fh[3] = f * a11 * uh[1] * uh[3] / (b11 + uh[1]) - a43 * uh[3] * uh[4] / (1. + a42*h42*uh[2]+a43*h43*uh[3]) - (d11 + e1) * uh[3]
    Fh[4] = f * a42 * uh[2] * uh[4] / (1. + a42*h42*uh[2]+a43*h43*uh[3]) + f * a43 * uh[3] * uh[4] / (1. + a42*h42*uh[2]+a43*h43*uh[3]) - (d2 + e2) * uh[4]
    Fh[5] = d1 * uh[2] + d11 * uh[3] + d2 * uh[4] -  (dD + e3) * uh[5]
end

# function, jacobian (autodiff) and initial guess
# equilibrium for initial node i
solh = nlsolve(fh!, [ initU1; .01; .01;.1;.3], autodiff = :forward,ftol=1.0e-16,iterations=1000000)
NLsolve.mcpsolve(fh!, [0.,0.,0.,0.0,0.0], [Inf,Inf,Inf,Inf,Inf], [ initU1; .01; .01;.1;.3], reformulation = :smooth, autodiff = :forward,ftol=1.0e-16) # returns -1.48382
NLsolve.mcpsolve(fh!, [0.,0.,0.,0.0,0.0], [Inf,Inf,Inf,Inf,Inf], [ initU1; .01; .01;.1;.3], reformulation = :minmax, autodiff = :forward,ftol=1.0e-16) # returns -1.0
NLsolve.mcpsolve(fh!, [0.,0.,0.,0.0,0.0], [Inf,Inf,Inf,Inf,Inf], [ initU1; .01; .01;.1;.3], reformulation = :smooth, autodiff = :forward,ftol=1.0e-16, method = :newton) # returns 6.25915e8


function ft!(Ft, ut)
    Ft[1] = In1 - rn * ut[1] -a1 * ut[1] * ut[2] / ( b1 + ut[1]) - a11 * ut[1] * ut[3] / ( b11 + ut[1]) + dD * ut[5] + df * solh.zero[1]
    Ft[2] = f * a1 * ut[1] * ut[2] / (b1 + ut[1]) - a42 * ut[2] * ut[4] / (1. + a42*h42*ut[2]+a43*h43*ut[3]) - (d1 + e1) * ut[2]
    Ft[3] = f * a11 * ut[1] * ut[3] / (b11 + ut[1]) - a43 * ut[3] * ut[4] / (1. + a42*h42*ut[2]+a43*h43*ut[3]) - (d11 + e1) * ut[3]
    Ft[4] = f * a42 * ut[2] * ut[4] / (1. + a42*h42*ut[2]+a43*h43*ut[3]) + f * a43 * ut[3] * ut[4] / (1. + a42*h42*ut[2]+a43*h43*ut[3]) - (d2 + e2) * ut[4]
    Ft[5] = d1 * ut[2] + d11 * ut[3] + d2 * ut[4] -  (dD + e3) * ut[5]
end


# function, jacobian (autodiff) and initial guess
# equilibrium for initial node i
solt = nlsolve(ft!, [initU1; .01; .01; .1;.3], autodiff = :forward, ftol = 1.0e-16, iterations = 1e6)
