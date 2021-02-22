mutable struct AdvDiffFRSpatial <: AdvFRSpatial
    p::Integer
    ns::Integer
    nf::Integer

    xs_ref::Vector{AbstractFloat}
    ws_ref::Vector{AbstractFloat}
    xf_ref::Vector{AbstractFloat}

    h::Real

    alpha::Real
    nu::Real
    kappa::Real

    Cp::Array{AbstractFloat, 2}
    C0::Array{AbstractFloat, 2}
    Cm::Array{AbstractFloat, 2}

    Bp2::Array{AbstractFloat, 2}
    Bp1::Array{AbstractFloat, 2}
    B0::Array{AbstractFloat, 2}
    Bm1::Array{AbstractFloat, 2}
    Bm2::Array{AbstractFloat, 2}

    function AdvDiffFRSpatial(p::Integer, h, alpha, nu, kappa)
        ns = p + 1
        nf = 2

        xs_ref, ws_ref = gausslegendre(ns)
        xf_ref = convert(typeof(xs_ref), [-1, 1])

        new(p, ns, nf, xs_ref, ws_ref, xf_ref, h, alpha, nu, kappa)
    end
end

@memoize function Bmatrix(X::T, kappa) where {T<:AdvDiffFRScheme}
    Cp, C0, Cm = Cmatrices(X, kappa)
    return Cp*Cp, C0*Cp + Cp*C0, Cm*Cp + C0*C0 + Cp*Cm, Cm*C0 + C0*Cm, Cm*Cm
end

function Qmatrix(FR::AdvDiffFRSpatial, k)
    FR.Cp, FR.C0, FR.Cm = Cmatrices(FR, FR.alpha)
    FR.Bp2, FR.Bp1, FR.B0, FR.Bm1, FR.Bm2 = Bmatrix(FR, FR.kappa)
    Qa = -(FR.Cm*exp(-1im*k*FR.h) + FR.Cp*exp(1im*k*FR.h) + FR.C0)*2/FR.h
    Qd = 4(FR.Bm2*exp(-im*k*2h) + FR.Bm1*exp(-im*k*h)+ FR.Bp2*exp(im*k*2h) + FR.Bp1*exp(im*k*h) + B0)/(FR.h*FR.h)
    2Qa - 4*FR.nu*Qd
end
