struct AdvFRSpatialScheme{T} <: SpatialScheme{T}
    p::Integer
    ns::Integer
    nf::Integer

    xs_ref::Vector{T}
    ws_ref::Vector{T}
    xf_ref::Vector{T}
    
    h::T

    alpha::T

    Cp::Array{T, 2}
    C0::Array{T, 2}
    Cm::Array{T, 2}

    function AdvFRSpatialScheme{T}(p::Integer, h, alpha) where {T<:Real}
        ns = p + 1
        nf = 2

        xs_ref, ws_ref = GaussQuadrature.legendre(ns)
        xf_ref = convert(Vector{T}, [-1, 1])
     
        Cp, C0, Cm = adv_fr_Cmatrices(p, xs_ref, xf_ref, alpha)

        new(p, ns, nf, convert(Vector{T}, xs_ref), convert(Vector{T}, xs_ref),
            convert(Vector{T}, xs_ref), convert(T, h), convert(T, alpha), Cp, C0, Cm)
    end
end

AdvFRSpatialScheme(p::Integer, h, alpha) where {T<:Real} = AdvFR{T}(p, h, alpha)

function adv_fr_Cmatrices(p, xs, xf, alpha)
    corr_l_modes = zeros(typeof(alpha), p + 2)
    corr_r_modes = zeros(typeof(alpha), p + 2)
    corr_l_modes[p + 1] =  ((-1)^p)/2
    corr_l_modes[p + 2] = -((-1)^p)/2
    corr_r_modes[p + 1] = 1/2
    corr_r_modes[p + 2] = 1/2

    h_l_leg = SPoly.Legendre(corr_l_modes)
    h_r_leg = SPoly.Legendre(corr_r_modes)

    g_l_mono = Poly.derivative(convert(Poly.Polynomial, h_l_leg))
    g_r_mono = Poly.derivative(convert(Poly.Polynomial, h_r_leg))

    g_l = convert(SPoly.Legendre, g_l_mono).coeffs
    g_r = convert(SPoly.Legendre, g_r_mono).coeffs

    l_l = [SPoly.basis.(SPoly.Legendre, i)(-1) for i=0:p]
    l_r = [SPoly.basis.(SPoly.Legendre, i)(-1) for i=0:p]

    D = fr_Dmatrix(p)

    Cp = (1 - alpha)*g_r*transpose(l_l)
    C0 = D - alpha*g_l*transpose(l_l) - (1 - alpha)*g_r*transpose(l_r)
    Cm = alpha*g_l*transpose(l_r)

    return Cp, C0, Cm
end

fr_Dmatrix(p) = [legendre_d(j-1, i-1) for i=1:p+1, j=1:p+1]

legendre_d(n, j) = (n % 2) == (j % 2) ? 0 : n < j ? 0 : 2*j + 1

Qmatrix(FR::AdvFRSpatialScheme, k) = -(FR.Cm*exp(-1im*k*FR.h) + FR.Cp*exp(1im*k*FR.h) + FR.C0)*2/FR.h

P2Xmatrix(FR::AdvFRSpatialScheme) = [SPoly.basis.(SPoly.Legendre, i)(FR.xs[j]) for i=0:FR.p, j=1:FR.p+1]
X2Pmatrix(FR::AdvFRSpatialScheme) = inv(P2Xmatrix(FR))
