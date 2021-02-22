module PDEFourierAnalysis

import LinearAlgebra; const LA = LinearAlgebra
import Polynomials; const Poly = Polynomials
import GaussQuadrature
import SpecialPolynomials; const SPoly = SpecialPolynomials

include("spatial_base.jl")
include("temporal_base.jl")

include("spatial_fr_adv.jl")

include("rk_tableaux.jl")

include("temporal_rk.jl")
include("temporal_erk.jl")
include("temporal_sdirk.jl")

include("temporal_dual.jl")

end
