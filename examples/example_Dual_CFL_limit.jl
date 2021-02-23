
include("PDEFourierAnalysis.jl")

import .PDEFourierAnalysis

dirk_stage = 1  # this can be changed later
pseudo_step = 2 # pseudo step number
dt = 1          # physical time step
dτ = dt/10      # pseudo time step

ERK = PDEFourierAnalysis.ERKTemporal("SSPRK_3")
DIRK = PDEFourierAnalysis.SDIRKTemporal("SDIRK_33(0.435866)", dirk_stage)
Dual = PDEFourierAnalysis.DualTemporal(DIRK, ERK, dt, dτ, pseudo_step)

p = 3
h = 1
alpha = 1.
FR = PDEFourierAnalysis.AdvFRSpatial(p, h, alpha)

System = PDEFourierAnalysis.FullScheme(FR, Dual)

println(PDEFourierAnalysis.CFL_limit(System))
