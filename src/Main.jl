
include("PDEFourierAnalysis.jl")

import .PDEFourierAnalysis

using Plots
import GR

plotly()

i = 1
m = 10
dt = 1
dtau = dt/10

ERK = PDEFourierAnalysis.ERKTemporal("SSPRK_3")
DIRK = PDEFourierAnalysis.SDIRKTemporal("SDIRK_33(0.435866)", 1)
Dual = PDEFourierAnalysis.DualTemporal(DIRK, ERK, dt, dtau, m)

p = 5
h = 1
alpha = 1.
FR = PDEFourierAnalysis.AdvFRSpatial(p, h, alpha)

P = PDEFourierAnalysis.amplification_factor_poly(ERK)

System = PDEFourierAnalysis.FullScheme(FR, ERK)
println(PDEFourierAnalysis.CFL_limit(System, n_k=1000, t_tol=1e-7, r_tol=1e-7))
#C, k = PDEFourierAnalysis.mod_wavenumber(FR)

#n_lx = 1000
#n_ly = 1001

#lx = range(start=-30, stop=0.5, length=n_lx)
#ly = range(start=-40, stop=40, length=n_ly)
#Lambda = lx'.*ones(n_ly) + ones(n_lx)'.*(1im.*ly)

#Sigma = PDEFourierAnalysis.amplification_factor.(Lambda, Ref(Dual))
#Sigma_abs = abs.(Sigma)
#Sigma_abs[Sigma_abs.>1] .= 1

#fig1 = heatmap(lx, ly, Sigma_abs, title=string("Dirk Stage ", i))
#p = heatmap(real.(Lambda), imag.(Lambda), abs.(Sigma))
#display(fig1)
