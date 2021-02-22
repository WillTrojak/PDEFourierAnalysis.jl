struct AdvFRSpatial <: FRSpatial
    p::Integer
    ns::Integer
    nf::Integer

    xs_ref::Vector{AbstractFloat}
    ws_ref::Vector{AbstractFloat}
    xf_ref::Vector{AbstractFloat}

    h::Real

    alpha::Real

    Cp::Array{AbstractFloat, 2}
    C0::Array{AbstractFloat, 2}
    Cm::Array{AbstractFloat, 2}

    function AdvFRSpatial(p::Integer, h, alpha)
        ns = p + 1
        nf = 2

        xs_ref, ws_ref = FGQ.gausslegendre(ns)
        xf_ref = convert(typeof(xs_ref), [-1, 1])

        new(p, ns, nf, xs_ref, ws_ref, xf_ref, h, alpha)
    end
end

@memoize function Cmatrices(X::T) where {T<:AdvFRSpatial}
    corr_l_modes = zeros(typeof(X.alpha), X.p + 2)
    corr_r_modes = zeros(typeof(X.alpha), X.p + 2)
    corr_l_modes[X.p + 1] =  ((-1)^X.p)/2
    corr_l_modes[X.p + 2] = -((-1)^X.p)/2
    corr_r_modes[X.p + 1] = 1/2
    corr_r_modes[X.p + 2] = 1/2

    h_l_leg = SPoly.Legendre(corr_l_modes)
    h_r_leg = SPoly.Legendre(corr_r_modes)

    g_l_mono = Poly.derivative(convert(Poly.Polynomial, h_l_leg))
    g_r_mono = Poly.derivative(convert(Poly.Polynomial, h_r_leg))

    g_l = convert(SPoly.Legendre, g_l_mono).coeffs
    g_r = convert(SPoly.Legendre, g_r_mono).coeffs

    l_l = [SPoly.basis.(SPoly.Legendre, i)(-1) for i=0:X.p]
    l_r = [SPoly.basis.(SPoly.Legendre, i)(-1) for i=0:X.p]

    D = Dmatrix(X)

    Cp = (1 - X.alpha)*g_r*transpose(l_l)
    C0 = D - X.alpha*g_l*transpose(l_l) - (1 - X.alpha)*g_r*transpose(l_r)
    Cm = X.alpha*g_l*transpose(l_r)

    return Cp, C0, Cm
end

function Qmatrix(FR::AdvFRSpatial, k)
    FR.Cp, FR.C0, FR.Cm = adv_fr_Cmatrices(FR)
    -(FR.Cm*exp(-1im*k*FR.h) + FR.Cp*exp(1im*k*FR.h) + FR.C0)*2/FR.h
end
