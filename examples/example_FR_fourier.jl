include("../src/PDEFourierAnalysis.jl")

import .PDEFourierAnalysis

p = 3
h = 1
alpha = 1.

FR = PDEFourierAnalysis.AdvFRSpatial(p, h, alpha)

C, k = PDEFourierAnalysis.mod_wavenumber(FR, nk=4)

println(C[1], " ", k[1])
println(C[2], " ", k[2])