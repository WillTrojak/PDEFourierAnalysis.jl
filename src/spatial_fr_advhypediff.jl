abstract type AbstractAdvHypeDiffFRSpatial <: AbstractAdvFRSpatial end

mutable struct AdvHypeDiffFRSpatial <: AbstractAdvFRSpatial
    p::Integer
    ns::Integer
    nf::Integer

    xs_ref::Vector{AbstractFloat}
    ws_ref::Vector{AbstractFloat}
    xf_ref::Vector{AbstractFloat}

    h::Real

    alpha::Real
    nu::Real
    T_r::Real
    kappa::Real

    Cp_a::Array{AbstractFloat, 2}
    C0_a::Array{AbstractFloat, 2}
    Cm_a::Array{AbstractFloat, 2}
    Cp_d::Array{AbstractFloat, 2}
    C0_d::Array{AbstractFloat, 2}
    Cm_d::Array{AbstractFloat, 2}

    function AdvHypeDiffFRSpatial(p::Integer, h, alpha, nu, T_r, kappa)
        ns = p + 1
        nf = 2

        xs_ref, ws_ref = gausslegendre(ns)
        xf_ref = convert(typeof(xs_ref), [-1, 1])

        new(p, ns, nf, xs_ref, ws_ref, xf_ref, h, alpha, nu, T_r, kappa)
    end
end

function Qmatrix(FR::AdvHypeDiffFRSpatial, k)
    FR.Cp_a, FR.C0_a, FR.Cm_a = Cmatrix(FR, FR.alpha)
    FR.Cp_d, FR.C0_d, FR.Cm_d = Cmatrix(FR, FR.kappa)
    Q = zeros(2FR.ns, 2FR.ns)
    Q[1:FR.ns, 1:FR.ns] += -(FR.Cm_a*exp(-1im*k*FR.h) + FR.Cp_a*exp(1im*k*FR.h) + FR.C0_a)*2/FR.h
    Q[1:FR.ns, FR.ns+1:2FR.ns] += FR.nu*(FR.Cm_a*exp(-1im*k*FR.h) + FR.Cp_a*exp(1im*k*FR.h) + FR.C0_a)*2/FR.h
    Q[FR.ns+1:2FR.ns, 1:FR.ns] += (FR.Cm_a*exp(-1im*k*FR.h) + FR.Cp_a*exp(1im*k*FR.h) + FR.C0_a)*2/(FR.h*FR.T_r)
    Q[FR.ns+1:2FR.ns, FR.ns+1:2FR.ns] += UniformScaling(1/FR.T_r)
    return Q
end

function t_derivative_switch(FR::AdvHypeDiffFRSpatial)
    A = zeros(2FR.ns,2FR.ns) + UniformScaling(1)
    A[FR.ns+1:2FR.ns, FR.ns+1:2FR.ns] -= UniformScaling(1)
    return A
end
