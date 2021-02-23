abstract type AbstractFRSpatial <: AbstractSpatial end

@memoize Dmatrix(X::T) where {T<:AbstractFRSpatial} = [legendre_d(j-1, i-1) for i=1:X.p+1, j=1:X.p+1]

@memoize P2Xmatrix(X::T) where {T<:AbstractFRSpatial} = [SpecialPolynomial.basis.(Legendre, i)(X.xs[j]) for i=0:X.p, j=1:X.p+1]
@memoize X2Pmatrix(X::T) where {T<:AbstractFRSpatial} = inv(P2Xmatrix(X))

@memoize function correction_modes(X::T) where {T<:AbstractFRSpatial}
    corr_l_modes = zeros(eltype(X.xs_ref), X.p + 2)
    corr_r_modes = zeros(eltype(X.xs_ref), X.p + 2)
    corr_l_modes[X.p + 1] =  ((-1)^X.p)/2
    corr_l_modes[X.p + 2] = -((-1)^X.p)/2
    corr_r_modes[X.p + 1] = 1/2
    corr_r_modes[X.p + 2] = 1/2

    h_l_leg = Legendre(corr_l_modes)
    h_r_leg = Legendre(corr_r_modes)

    g_l_mono = derivative(convert(Polynomial, h_l_leg))
    g_r_mono = derivative(convert(Polynomial, h_r_leg))

    g_l = convert(Legendre, g_l_mono).coeffs
    g_r = convert(Legendre, g_r_mono).coeffs
    return g_l, g_r
end
