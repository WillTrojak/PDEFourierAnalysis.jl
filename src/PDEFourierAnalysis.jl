module PDEFourierAnalysis

using Combinatorics: permutations
using FastGaussQuadrature: gausslegendre
using LinearAlgebra
using Memoize
using Polynomials
using SpecialPolynomials

include("utils_polynomials.jl")
include("utils_matrix.jl")

include("spatial_base.jl")
include("spatial_fr.jl")
include("spatial_fr_adv.jl")
include("spatial_fr_advdiff.jl")
include("spatial_fr_advhypediff.jl")

include("rk_tableaux.jl")

include("temporal_base.jl")
include("temporal_rk.jl")
include("temporal_erk.jl")
include("temporal_sdirk.jl")
include("temporal_dual.jl")

include("scheme_base.jl")
include("scheme_full.jl")

end # module
