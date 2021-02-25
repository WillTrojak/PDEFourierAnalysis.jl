abstract type AbstractAdvFRSpatial <: AbstractFRSpatial end

mutable struct AdvFRSpatial <: AbstractAdvFRSpatial
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

        xs_ref, ws_ref = gausslegendre(ns)
        xf_ref = convert(typeof(xs_ref), [-1, 1])

        new(p, ns, nf, xs_ref, ws_ref, xf_ref, h, alpha)
    end
end


@memoize function Cmatrices(X::T, alpha) where {T<:AbstractAdvFRSpatial}
    g_l, g_r = correction_modes(X)

    l_l = [SpecialPolynomials.basis.(Legendre, i)(-1) for i=0:X.p]
    l_r = [SpecialPolynomials.basis.(Legendre, i)( 1) for i=0:X.p]

    D = Dmatrix(X)

    Cp = (1 - alpha)*g_r*transpose(l_l)
    C0 = D - alpha*g_l*transpose(l_l) - (1 - alpha)*g_r*transpose(l_r)
    Cm = alpha*g_l*transpose(l_r)

    return Cp, C0, Cm
end

function Qmatrix(FR::AdvFRSpatial, k)
    FR.Cp, FR.C0, FR.Cm = Cmatrices(FR, FR.alpha)
    -(FR.Cm*exp(-1im*k*FR.h) + FR.Cp*exp(1im*k*FR.h) + FR.C0)*2/FR.h
end

function mod_wavenumber(FR::AdvFRSpatial; nk=100, k_min=1e-6, k_max=2*pi)
    k = LinRange(k_min*(FR.p + 1), k_max*(FR.p + 1), nk)
    C = zeros(ComplexF64, FR.p + 1, nk)
    for i=1:nk
        C[:, i] = k[i].*eigvals(-Qmatrix(FR, k[i])/(im*k[i]))
    end
    organise_rows!(C)
    return C, k
end

t_derivative_switch(FR::AdvFRSpatial) = UniformScaling(1)
