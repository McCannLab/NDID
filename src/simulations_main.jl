# -------------- Main results  --------------

include("ndid_core.jl")
colC = "black"
colD = "#c7254e"
colN = "#ffdd55"
colRe = "#3fb3b2"
colRle = "#8edb7c"

# Bifurcation analysis (increasing df from 0 to 1)

## A- All nodes are fast

### i- NDCRe
param = EcoPar()

### results for all components
bife = bif_analysis2(0.0:0.01:1.0, eco_NDCRe!, fill(0.03, 20), (0.0, 30000.0),
    25000:30000, param)
### bifurcation analysis for the all components
fig_bif(bife, [1 13 17], "figRe_N_30.svg", "df", "N (max/min)", 0, 12.5)
fig_bif(bife, [2 14 18], "figRe_R_30.svg", "df", "Re (max/min)", 0, 2, colRe, 2)
fig_bif(bife, [3 15 19], "figRe_C_30.svg", "df", "C (max/min)", 0, 4)
fig_bif(bife, [4 16 20], "figRe_D_30.svg", "df", "D (max/min)", 0, 2)
### combining 2 results (not used in the last version of the article)
fig_bif2(bife, [1 13 17], [3 15 19], [colN colC], "figRe_p1.svg", "", "",  0, 12.5)
fig_bif2(bife, [1 13 17], [4 16 20], [colN colD], "figRe_p2.svg", "", "",  0, 10.5)
### actual results
fig_bif(bife, [2 14 18], "figRe_Re.svg", "", "", 0, 2, colRe, 3)
barplotmean(bife, [1 13 17], "hist_bife_N.svg", 0, 11, colN)
barplotmean(bife, [4 16 20], "hist_bife_D.svg", 0.170, 0.180, colD)



### ii- NDCReRle
param = EcoPar()
param.a42 = 0.15
param.a43 = 0.05
# NB
# a43=.15; a42=.05 -- loss of edible R
# a43=.15; a42=0.01 -- no loss but edible weakly declines while Ri rises -- this is universal
# a42= .15; a42=0; Redible flatlines and Ri rises

### bifurcation analysis for the all components
bifle = bif_analysis2(0.0:0.01:1.0, eco_NDCReRle!, fill(0.03, 25), (0.0, 30000.0), 25000:30000, param)
### results for all components
fig_bif(bifle, [1 16 21], "figRle_N.svg", "df", "N max/min", 0, 1.5)
fig_bif(bifle, [2 17 22], "figRle_Re.svg", "df", "Re max/min", 0, .6, colRe, 2)
fig_bif(bifle, [3 18 23], "figRle_Rle.svg", "df", "Rle max/min", 0, 1.5)
fig_bif(bifle, [4 19 24], "figRle_C.svg", "df", "C max/min", 0, 2)
fig_bif(bifle, [5 20 25], "figRle_D.svg", "df", "D max/min", 0, 1.5)
### combining 2 or 3 results (not used in the last version of the article)
fig_bif2(bifle, [1 16 21], [5 20 25], [colN colD], "figRle_p2.svg", "", "",  0, 1)
fig_bif3(bifle, ([2 14 18], [4 16 20], [3 18 23]), [colRe colD colRle], "figRle_p2.svg", "", "Max/min",  0, 1.5)
### actual results
fig_bif(bifle, [2 17 22], "figRle_Re.svg", "", "", 0, .6, colRe, 3)
barplotmean(bifle, [1 16 21], "hist_bifle_N.svg", 0, .15, colN)
barplotmean(bifle, [5 20 25], "hist_bifle_D.svg", 0, 1, colD)





## B- All nodes are fast but the hub node that is slow

### i- NDCRe

param = EcoPar()

### bifurcation analysis for the all components
bife_sh = bif_analysis2(0.0:0.01:1.0, eco_NDCRe_sh!, fill(0.03, 20),
    (0.0, 30000.0), 25000:30000, param)

### results for all components
fig_bif(bife_sh, [1 13 17], "figRe_N_sh.svg", "df", "N (max/min)", 0, 20)
fig_bif(bife_sh, [2 14 18], "figRe_R_sh.svg", "df", "Re (max/min)", 0, .8)
fig_bif(bife_sh, [3 15 19], "figRe_C_sh.svg", "df", "C (max/min)", 0, 10)
fig_bif(bife_sh, [4 16 20], "figRe_D_sh.svg", "df", "D (max/min)", 0, 1)

### combining more than 2 or results
fig_bif2(bife_sh, [1 13 17], [3 15 19], [colN colD], "figRe_sh_p1.svg", "", "Max/min",  0, 9.2)
fig_bif2(bife_sh, [2 14 18], [4 16 20], [colRe colD], "figRe_sh_p2.svg", "", "Max/min",  0, .8)
fig_bif2(bife_sh, [1 13 17], [4 16 20], [colN colD], "figRe_sh_p2.svg", "", "",  0, 8.5)
### actual results
fig_bif(bife_sh, [2 14 18], "figRe_R_sh.svg", "", "", 0, .5, colRe, 3)
barplotmean(bife_sh, [1 13 17], "hist_bife_sh_N.svg", 0, 6, colN)
barplotmean(bife_sh, [4 16 20], "hist_bife_sh_D.svg", 0.17, 0.18, colD)



### ii- NDCReRle

param = EcoPar()
param.a42 = 0.15
param.a43 = 0.05

bifle_sh = bif_analysis2(0.0:0.01:1.0, eco_NDCReRle_sh!, fill(0.03, 25), (0.0, 30000.0), 25000:30000, param)

fig_bif(bifle_sh, [1 16 21], "figRle_N_sh.svg", "df", "N (max/min)", 0, 1)
fig_bif(bifle_sh, [2 17 22], "figRle_Re_sh.svg", "df", "Re (max/min)", 0, .8)
fig_bif(bifle_sh, [3 18 23], "figRle_Rle_sh.svg", "df", "Rle (max/min)", 0, 1)
fig_bif(bifle_sh, [4 19 24], "figRle_C_sh.svg", "df", "C (max/min)", 0, 2)
fig_bif(bifle_sh, [5 20 25], "figRle_D_sh.svg", "df", "D (max/min)", 0, 1)

### combining more than 2 or results
fig_bif2(bifle_sh, [1 16 21], [3 15 19], [colN colC], "figRle_sh_p1.svg", "df", "Max/min",  0, 1.45)
fig_bif3(bifle_sh, ([2 14 18], [4 16 20], [3 18 23]), [colRe colD colRle], "figRle_sh_p2.svg", "df", "Max/min",  0, 1.45)

### actual results
fig_bif(bifle_sh, [2 17 22], "figRle_sh_Re.svg", "", "", 0, .6, colRe, 3)
barplotmean(bifle_sh, [1 16 21], "hist_bifle_sh_N.svg", 0, .17, colN)
barplotmean(bifle_sh, [5 20 25], "hist_bifle_sh_D.svg", 0, 1, colD)
